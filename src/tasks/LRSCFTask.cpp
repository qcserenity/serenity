/**
 * @file   LRSCFTask.cpp
 *
 * @date   Aug 17, 2016
 * @author M. Boeckers, N. Niemeyer, J. Toelle
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "tasks/LRSCFTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h" //keep engines alive.
#include "postHF/LRSCF/Analysis/DeltaSpinSquared.h"
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
#include "postHF/LRSCF/Analysis/ExcitationSpectrum.h"
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"
#include "postHF/LRSCF/Analysis/LRSCFPopulationAnalysis.h"
#include "postHF/LRSCF/Analysis/ResponseProperties.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Sigmavectors/RICC2/XWFController.h"
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "postHF/LRSCF/Tools/CouplingConstruction.h"
#include "postHF/LRSCF/Tools/EigenvalueSolver.h"
#include "postHF/LRSCF/Tools/LRSCFRestart.h"
#include "postHF/LRSCF/Tools/LRSCFSetup.h"
#include "postHF/LRSCF/Tools/NonlinearEigenvalueSolver.h"
#include "postHF/LRSCF/Tools/ResponseLambda.h"
#include "postHF/LRSCF/Tools/ResponseSolver.h"
#include "postHF/LRSCF/Tools/SigmaCalculator.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LRSCFTask<SCFMode>::LRSCFTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                              const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _act(activeSystems), _env(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
void LRSCFTask<SCFMode>::run() {
  std::cout << std::scientific << std::setprecision(3);

  printSectionTitle("LRSCF");

  // set type of calculation
  Options::LRSCF_TYPE type = Options::LRSCF_TYPE::ISOLATED;
  if (_act.size() == 1 && _env.empty()) {
    type = Options::LRSCF_TYPE::ISOLATED;
  }
  else if (_act.size() == 1 && _env.size() >= 1) {
    type = Options::LRSCF_TYPE::UNCOUPLED;
  }
  else if (_act.size() > 1) {
    type = Options::LRSCF_TYPE::COUPLED;
  }

  // creates lrscf controller that contains all data needed for the LRSCF calculation
  for (const auto& system : _act) {
    _lrscf.push_back(std::make_shared<LRSCFController<SCFMode>>(system, settings));
  }
  // if special coupling pattern is requested
  if (!settings.couplingPattern.empty()) {
    const unsigned int numberOfRows = std::sqrt(settings.couplingPattern.size());
    Eigen::MatrixXi couplingPattern(numberOfRows, numberOfRows);
    couplingPattern.setZero();
    for (unsigned int i = 0; i < numberOfRows; i++) {
      if (i < _act.size()) {
        _referenceLoadingType.push_back(Options::LRSCF_TYPE::COUPLED);
      }
      else {
        _referenceLoadingType.push_back(settings.loadType);
      }
      for (unsigned int j = 0; j < numberOfRows; j++) {
        couplingPattern(i, j) = settings.couplingPattern[(i * numberOfRows) + j];
      }
    }
    _couplingPattern = couplingPattern;
    printf("  ----------------------------------------------------------------\n");
    printf("                 Special coupling pattern is used!                \n");
    printf("  ----------------------------------------------------------------\n");
    std::cout << _couplingPattern << std::endl;
    // adds additional lrscfController if one system is used coupled and uncoupled
    if ((unsigned int)_couplingPattern.rows() > _act.size()) {
      for (unsigned int iRow = _act.size(); iRow < _couplingPattern.rows(); iRow++) {
        // the systemcontroller which occurs two times
        unsigned int index = _couplingPattern(iRow, iRow) - 1;
        _lrscf.push_back(std::make_shared<LRSCFController<SCFMode>>(_act[index], settings));
      }
    }
  }
  // regular coupled calculation
  else if (type == Options::LRSCF_TYPE::COUPLED && settings.couplingPattern.empty()) {
    Eigen::MatrixXi couplingPattern(_act.size(), _act.size());
    couplingPattern.setZero();
    for (unsigned int i = 0; i < _act.size(); i++) {
      couplingPattern(i, i) = i + 1;
      _referenceLoadingType.push_back(settings.loadType);
    }
    _couplingPattern = couplingPattern;
  }
  // sets the environment systems, important for EOSigmavector with more than two subsystems
  for (auto lrscf : _lrscf) {
    lrscf->setEnvSystems(_env);
  }

  LRSCFSetup<SCFMode>::setupLRSCFController(settings, _act, _env, _lrscf);

  Eigen::VectorXd diagonal = LRSCFSetup<SCFMode>::getDiagonal(_lrscf);
  unsigned nDimension = diagonal.size();
  if (nDimension == 0) {
    throw SerenityError("Orbital-transition space of zero is impossible.");
  }
  settings.nEigen = std::min(settings.nEigen, nDimension);

  // sanity check for gauge-origin input
  if (settings.gaugeOrigin.size() != 3) {
    throw SerenityError("Gauge-origin expects three cartesian coordinates!");
  }

  // information about the TDDFT calculation
  LRSCFSetup<SCFMode>::printInfo(_lrscf, settings, _env, type);

  // dipole integrals (electric (length/velocity) and magnetic)
  auto gaugeOrigin = LRSCFSetup<SCFMode>::getGaugeOrigin(settings, _act, _env);
  auto dipoles = std::make_shared<DipoleIntegrals<SCFMode>>(_lrscf, gaugeOrigin);

  // prepare solutionvectors, eigenvectors and eigenvalues
  std::shared_ptr<std::vector<Eigen::MatrixXd>> solutionvectors = nullptr;
  std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors = nullptr;
  std::shared_ptr<std::vector<Eigen::MatrixXd>> transDensityMatrices = nullptr;
  Eigen::VectorXd eigenvalues;

  // initialize restart system
  LRSCFRestart<SCFMode> restart(_lrscf, settings, type);

  // prepare solutions vector in case of a FDEc calculation
  if (type == Options::LRSCF_TYPE::COUPLED && settings.nEigen > 0 && !settings.restart) {
    eigenvectors = LRSCFSetup<SCFMode>::setupFDEcTransformation(settings, _couplingPattern, _referenceLoadingType,
                                                                _lrscf, _act, nDimension);
    // set number of the eigenvalues to be determined in the coupled step
    settings.nEigen = (*eigenvectors)[0].cols();
  }
  else if (settings.restart) {
    eigenvectors = restart.fetchEigenpairs(eigenvalues);
  }

  // Turn off grid output globally.
  bool oldIOGridInfo = iOOptions.printGridInfo;
  iOOptions.printGridInfo = (GLOBAL_PRINT_LEVEL != Options::GLOBAL_PRINT_LEVELS::DEBUGGING) ? false : true;

  // Check if local orbitals are used.
  bool fockDiagonal = true;
  for (auto& lrscf : _lrscf) {
    fockDiagonal = fockDiagonal && lrscf->isMOFockMatrixDiagonal();
  }
  if (!fockDiagonal) {
    WarningTracker::printWarning(" Fock matrix is not diagonal, please make sure this is intended.", true);
  }

  // Setup lambda functions for response matrix sigmavectors.
  ResponseLambda<SCFMode> lambda(_act, _env, _lrscf, diagonal, settings);
  lambda.setupTDDFTLambdas();
  SigmaCalculator sigmaCalculator = lambda.setEigensolverMode(fockDiagonal);

  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  if (lambda.usesLRExchange()) {
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 4);
  }

  // Eigenvalue Solver.
  if (settings.nEigen > 0) {
    // Adjust maximum subspace dimension if it was untouched in the input.
    if (settings.maxSubspaceDimension == 1e9) {
      settings.maxSubspaceDimension = 25 * settings.nEigen;
    }

    // Set initial subspace size.
    unsigned initialSubspace = std::min(settings.nEigen * 2, settings.nEigen + 8);
    initialSubspace = (type == Options::LRSCF_TYPE::COUPLED) ? settings.nEigen : initialSubspace;
    initialSubspace = std::min(initialSubspace, nDimension);

    if (settings.method == Options::LR_METHOD::TDA || settings.method == Options::LR_METHOD::TDDFT) {
      if (type == Options::LRSCF_TYPE::COUPLED && settings.partialResponseConstruction && !settings.fullFDEc) {
        CouplingConstruction<SCFMode>::solve(_lrscf, settings, _referenceLoadingType, sigmaCalculator, eigenvectors, eigenvalues);
      }
      else {
        // If this is a coupled calculation and no full FDEc was requested, only perform a single iteration.
        EigenvalueSolver eigensolver(settings.saveResponseMatrix, nDimension, settings.nEigen, diagonal, settings.conv,
                                     (type == Options::LRSCF_TYPE::COUPLED && !settings.fullFDEc) ? 1 : settings.maxCycles,
                                     settings.maxSubspaceDimension, initialSubspace, settings.algorithm,
                                     sigmaCalculator, eigenvectors, restart.getWriteToDisk());

        eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(eigensolver.getEigenvectors());
        eigenvalues = eigensolver.getEigenvalues();

        // Switch to finer grid if necessary and perform one final iteration.
        if (lambda.usesKernel() && settings.grid.accuracy != settings.grid.smallGridAccuracy) {
          printSmallCaption("Switching from small to default grid");
          printf("  Residual norms might not match the desired convergence anymore.\n");
          printf("  Grid accuracy:    %1i -> %1i\n\n", settings.grid.smallGridAccuracy, settings.grid.accuracy);
          lambda.setupKernel(Options::GRID_PURPOSES::DEFAULT);
          // At this point, there is X and Y stored in (*eigenvectors)[0] and [1], respectively.
          if (settings.method == Options::LR_METHOD::TDDFT) {
            (*eigenvectors)[0] += (*eigenvectors)[1];
            (*eigenvectors)[1] = ((*eigenvectors)[0] - 2 * (*eigenvectors)[1]).eval();
          }
          eigensolver = EigenvalueSolver(settings.saveResponseMatrix, nDimension, settings.nEigen, diagonal,
                                         settings.conv, 1, settings.maxSubspaceDimension, settings.nEigen,
                                         settings.algorithm, sigmaCalculator, eigenvectors, restart.getWriteToDisk());
          eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(eigensolver.getEigenvectors());
          eigenvalues = eigensolver.getEigenvalues();
        }
      }

      // For (TDA-)TDDFT transition density matrices: (X+Y) (X-Y).
      transDensityMatrices = std::make_shared<std::vector<Eigen::MatrixXd>>(2);
      (*transDensityMatrices)[0] = (*eigenvectors)[0];
      (*transDensityMatrices)[1] = (*eigenvectors)[0];
      if (settings.method == Options::LR_METHOD::TDDFT) {
        (*transDensityMatrices)[0] += (*eigenvectors)[1];
        (*transDensityMatrices)[1] -= (*eigenvectors)[1];
      }

      // If double hybrid is used, calculate CIS(D) correction.
      if (lambda.usesDoubleHybrid()) {
        lambda.setupCC2Lambdas();
        Eigen::VectorXd cisd = eigenvalues;

        NonlinearEigenvalueSolver nlEigenSolver(settings.nEigen, diagonal, settings.conv, settings.preopt,
                                                settings.diis, settings.diisStore, settings.maxCycles,
                                                Options::LR_METHOD::CISD, (*transDensityMatrices)[0], cisd,
                                                lambda.getRightXWFSigma(), restart.getWriteToDisk());
        nlEigenSolver.solve();

        // Add correction to excitation energies.
        eigenvalues += lambda.getDHRatio() * cisd;

        // Get rid of CC2 and RIIntegral controller.
        for (auto& lrscf : _lrscf) {
          lrscf->finalizeXWFController();
          lrscf->finalizeRIIntegrals(LIBINT_OPERATOR::coulomb);
        }
      } /* Double hybrid correction */
    }
    else {
      lambda.setupCC2Lambdas();
      if (!eigenvectors || settings.method == Options::LR_METHOD::CISD) {
        // RI-CIS calculation.
        printSubSectionTitle("RI-CIS Calculation");
        double cisthresh = (settings.method == Options::LR_METHOD::CISD) ? settings.conv : settings.preopt;
        EigenvalueSolver eigensolver(false, nDimension, settings.nEigen, diagonal, cisthresh, settings.maxCycles,
                                     settings.maxSubspaceDimension, initialSubspace,
                                     Options::RESPONSE_ALGORITHM::SYMMETRIC, sigmaCalculator);

        eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(eigensolver.getEigenvectors());
        eigenvalues = eigensolver.getEigenvalues();
      }
      /* RI-CIS calculation. */

      Eigen::VectorXd nleigenvalues(settings.nEigen);
      printSubSectionTitle("Right Eigenvector Calculation");
      NonlinearEigenvalueSolver nlEigenSolver(settings.nEigen, diagonal, settings.conv, settings.preopt, settings.diis,
                                              settings.diisStore, settings.maxCycles, settings.method, (*eigenvectors)[0],
                                              nleigenvalues, lambda.getRightXWFSigma(), restart.getWriteToDisk());
      nlEigenSolver.solve();

      if (settings.method == Options::LR_METHOD::CISD) {
        // In the CISD case, this is just the correction.
        eigenvalues += nleigenvalues;
        transDensityMatrices = std::make_shared<std::vector<Eigen::MatrixXd>>(2, (*eigenvectors)[0]);
      }
      else {
        eigenvalues = nleigenvalues;
        if (settings.ccprops) {
          Eigen::VectorXd rightEigenvalues = eigenvalues;
          if (settings.method != Options::LR_METHOD::ADC2) {
            eigenvectors->push_back((*eigenvectors)[0]);
            printSubSectionTitle("Left Eigenvector Calculation");
            NonlinearEigenvalueSolver nlEigenSolver(settings.nEigen, diagonal, settings.conv, settings.preopt,
                                                    settings.diis, settings.diisStore, settings.maxCycles, settings.method,
                                                    (*eigenvectors)[1], eigenvalues, lambda.getLeftXWFSigma());
            nlEigenSolver.solve();
            _lrscf[0]->getXWFController()->normalizeEigenvectors(*eigenvectors, eigenvalues);
          }
          else {
            _lrscf[0]->getXWFController()->normalizeEigenvectors(*eigenvectors, eigenvalues);
          }

          if ((eigenvalues - rightEigenvalues).cwiseAbs().maxCoeff() > settings.conv) {
            printf("  Careful: left and right eigenvalues are not converged to %5.2e Hartree!\n", settings.conv);
          }

          eigenvectors->push_back(Eigen::MatrixXd::Zero(nDimension, settings.nEigen));
          if (settings.method == Options::LR_METHOD::CC2) {
            printSubSectionTitle("Ground-state Langrange Multiplier");
            _lrscf[0]->getXWFController()->calcGroundStateLagrangeMultiplier();
            printSubSectionTitle("Transition-moment Langrange Multiplier");
            _lrscf[0]->getXWFController()->calcTransitionMomentLagrangeMultiplier(*eigenvectors, eigenvalues);
          }

          transDensityMatrices = std::make_shared<std::vector<Eigen::MatrixXd>>(2);
          _lrscf[0]->getXWFController()->calcDensityMatrices(*eigenvectors, eigenvalues, *transDensityMatrices);
        }
        else {
          _lrscf[0]->getXWFController()->normalizeEigenvectors(*eigenvectors, eigenvalues);
        }
      }
    }
    restart.storeConvergedSolution(eigenvectors, eigenvalues);
  } /* Eigenvalue calculation. */

  // Response Solver.
  if (!settings.frequencies.empty() || !settings.frequencyRange.empty()) {
    // frequencyRange input.
    if (!settings.frequencyRange.empty()) {
      // Sanity check.
      bool freqsGood = settings.frequencyRange.size() == 3;
      freqsGood = freqsGood && settings.frequencyRange[0] < settings.frequencyRange[1];
      freqsGood = freqsGood && (settings.frequencyRange[0] + settings.frequencyRange[1]) > settings.frequencyRange[2];
      if (!freqsGood) {
        throw SerenityError("LRSCFTask: frequencyRange input error!");
      }

      for (double freq = settings.frequencyRange[0]; freq <= settings.frequencyRange[1]; freq += settings.frequencyRange[2]) {
        settings.frequencies.push_back(freq);
      }
    }

    // It's better to work with a.u. internally.
    for (auto& freq : settings.frequencies) {
      freq *= EV_TO_HARTREE;
    }
    settings.damping *= EV_TO_HARTREE;

    // Append static limit if velocity gauge.
    if (settings.gauge == Options::GAUGE::VELOCITY) {
      settings.frequencies.push_back(0);
    }

    settings.maxSubspaceDimension = settings.frequencies.size() * 100;

    // Just use normal grid for property calculation.
    if (lambda.usesKernel() && settings.grid.accuracy != settings.grid.smallGridAccuracy && settings.nEigen < 1) {
      printf("  Will use default grid accuracy for response property calculation.\n\n");
      lambda.setupKernel(Options::GRID_PURPOSES::DEFAULT);
    }

    // Determine rhs of the linear equation system.
    std::vector<Eigen::MatrixXd> rhs(settings.damping != 0 ? 4 : 2, Eigen::MatrixXd::Zero(nDimension, 3));
    if (settings.gauge == Options::GAUGE::LENGTH) {
      rhs[0] = 2 * (*dipoles->getLengths());
    }
    else {
      rhs[1] = 2 * (*dipoles->getVelocities());
    }

    ResponseSolver responseSolver(diagonal, settings.conv, settings.maxCycles, settings.maxSubspaceDimension,
                                  settings.frequencies, settings.damping, rhs, lambda.getRPASigma(), nullptr,
                                  [](std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&) {});

    // Get converged solution vectors for analysis.
    solutionvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(responseSolver.getEigenvectors());
    Eigen::VectorXd freq(settings.frequencies.size());
    for (unsigned afreq = 0; afreq < settings.frequencies.size(); afreq++) {
      freq[afreq] = settings.frequencies[afreq];
    }
    restart.storeConvergedResponse(solutionvectors, freq);
  } /* Response calculation. */

  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  if (lambda.usesLRExchange()) {
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 4);
  }

  // Analysis: dominant contributions, excitation spectrum, CD Spectrum, response properties.
  if (settings.analysis) {
    if (eigenvectors) {
      _excitations = Eigen::MatrixXd::Zero(settings.nEigen, 6);
      // Prints dominant contributions etc.
      LRSCFAnalysis<SCFMode>::printDominantContributions(_lrscf, (*eigenvectors), eigenvalues, settings.dominantThresh);
      // Print delta s^2 values in unrestricted case.
      if (settings.method == Options::LR_METHOD::TDA || settings.method == Options::LR_METHOD::TDDFT) {
        if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
          DeltaSpinSquared<SCFMode> spinSquared(_lrscf, settings.nEigen, settings.method, type);
          spinSquared.print();
        }
      }
      if (transDensityMatrices) {
        // Prints excitation and CD spectrum (oscillator and rotatory strengths, respectively).
        std::string name = (type != Options::LRSCF_TYPE::COUPLED)
                               ? _lrscf[0]->getSys()->getSystemPath() + _lrscf[0]->getSys()->getSystemName()
                               : "";
        ExcitationSpectrum<SCFMode>::printSpectrum(settings.method, dipoles, (*transDensityMatrices), eigenvalues,
                                                   _excitations, name);
      }
      // Calculate Löwdin charges and store on disk.
      if (settings.transitionCharges) {
        if (settings.method == Options::LR_METHOD::CC2) {
          WarningTracker::printWarning("CC2 transition densities might not be well defined.", true);
        }
        LRSCFPopulationAnalysis<SCFMode>::calculateTransitionCharges(_lrscf, (*transDensityMatrices));
      }
    }
    if (solutionvectors) {
      // Calculate specific rotation factor (important for ORD).
      double molWeight = 0.0;
      for (const auto& sys : _act) {
        auto atoms = sys->getGeometry()->getAtoms();
        for (const auto& atom : atoms) {
          if (!atom->isDummy()) {
            molWeight += atom->getAtomType()->getMass();
          }
        }
      }
      // See Grimme Chem. Phys. Lett. 361 (2002) 321-328.
      double rotFactor = 1.343e-4 / molWeight * HARTREE_TO_OOCM * HARTREE_TO_OOCM;

      // Print properties.
      ResponseProperties<SCFMode>::printProperties(dipoles, (*solutionvectors), settings.frequencies, rotFactor,
                                                   settings.damping, settings.gauge, _properties);
    }
  }
  iOOptions.printGridInfo = oldIOGridInfo;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<LRSCFController<SCFMode>>> LRSCFTask<SCFMode>::getLRSCFControllers() {
  return _lrscf;
}

template<Options::SCF_MODES SCFMode>
const Eigen::MatrixXd& LRSCFTask<SCFMode>::getTransitions() {
  return _excitations;
}

template<Options::SCF_MODES SCFMode>
const std::vector<std::tuple<double, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d>>&
LRSCFTask<SCFMode>::getProperties() {
  return _properties;
}

template class LRSCFTask<Options::SCF_MODES::RESTRICTED>;
template class LRSCFTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
