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
#include "geometry/gradients/TDDFTGradientCalculator.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h" //keep engines alive.
#include "io/FormattedOutputStream.h"
#include "parameters/Constants.h" // EV_TO_HARTREE
#include "postHF/LRSCF/Analysis/DeltaSpinSquared.h"
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
#include "postHF/LRSCF/Analysis/ExcitationSpectrum.h"
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"
#include "postHF/LRSCF/Analysis/LRSCFPopulationAnalysis.h"
#include "postHF/LRSCF/Analysis/ResponseProperties.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/RICC2/CC2Controller.h"
#include "postHF/LRSCF/RICC2/CC2HelperFunctions.h"
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "postHF/LRSCF/Tools/CouplingConstruction.h"
#include "postHF/LRSCF/Tools/EigenvalueSolver.h"
#include "postHF/LRSCF/Tools/LRSCFRestart.h"
#include "postHF/LRSCF/Tools/LRSCFSetup.h"
#include "postHF/LRSCF/Tools/NonlinearEigenvalueSolver.h"
#include "postHF/LRSCF/Tools/NonlinearResponseSolver.h"
#include "postHF/LRSCF/Tools/ResponseLambda.h"
#include "postHF/LRSCF/Tools/ResponseSolver.h"
#include "postHF/LRSCF/Tools/SigmaCalculator.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

void LRSCFTaskSettings::printSettings(std::basic_string<char> filename) {
  std::string field;
  std::string value;
  std::ofstream ofs;
  ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
  ofs << "#=======================================================" << std::endl;
  ofs << "# NOTE:" << std::endl;
  ofs << "# This file contains the list of LRSCFTaskSettings " << std::endl;
  ofs << "#  that were at some point used for this system." << std::endl;
  ofs << "# To reuse this in an input, note that the active and " << std::endl;
  ofs << "#  passive systems need to be added to this manually." << std::endl;
  ofs << "#=======================================================" << std::endl;
  ofs << "+task lrscf" << std::endl;
  print_visitor visitor(field, value, ofs);
  visit_each((*this), visitor);
  ofs << "+emb" << std::endl;
  visit_each(this->embedding, visitor);
  ofs << "-emb" << std::endl;
  ofs << "+grid" << std::endl;
  visit_each(this->grid, visitor);
  ofs << "-grid" << std::endl;
  ofs << "+customFunc" << std::endl;
  visit_each(this->customFunc, visitor);
  ofs << "-customFunc" << std::endl;
  ofs << "-task" << std::endl;
  ofs.close();
}

template<Options::SCF_MODES SCFMode>
LRSCFTask<SCFMode>::LRSCFTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                              const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _act(activeSystems), _env(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
void LRSCFTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _act);
  printSectionTitle("LRSCF");

  // Save settings.
  LRSCFTaskSettings oldSettings = settings;

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

  bool isNotCC2 = (settings.method == Options::LR_METHOD::TDA || settings.method == Options::LR_METHOD::TDDFT);

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

  // Sanity check for gauge-origin input.
  if (settings.gaugeOrigin.size() != 3) {
    throw SerenityError("Gauge-origin expects three cartesian coordinates!");
  }

  // Approximate Coulomb interaction printout.
  LRSCFSetup<SCFMode>::printApproximateCoulombInfo(_lrscf, settings);

  // Stability analysis printout.
  if (settings.scfstab != Options::STABILITY_ANALYSIS::NONE) {
    LRSCFSetup<SCFMode>::prepareStabilityAnalysis(_lrscf, settings);
  }

  // Information about the TDDFT calculation.
  LRSCFSetup<SCFMode>::printInfo(_lrscf, settings, _env, type);

  Eigen::VectorXd diagonal = LRSCFSetup<SCFMode>::getDiagonal(_lrscf);
  unsigned nDimension = diagonal.size();
  if (nDimension == 0) {
    throw SerenityError("Orbital-transition space of zero is impossible.");
  }
  settings.nEigen = std::min(settings.nEigen, nDimension);

  // Dipole integrals (electric (length/velocity) and magnetic).
  auto gaugeOrigin = LRSCFSetup<SCFMode>::getGaugeOrigin(settings, _act, _env);
  auto dipoles = std::make_shared<DipoleIntegrals<SCFMode>>(_lrscf, gaugeOrigin);

  // Calculate molecular weight (important for ORD) and nuclear part of dipole moment.
  double molWeight = 0;
  Eigen::Vector3d nucDipoleMoment = Eigen::Vector3d::Zero();
  LRSCFSetup<SCFMode>::calculateMolecularWeightandNuclearDipoleMoment(_act, _env, gaugeOrigin, molWeight, nucDipoleMoment);

  // Prepare solutionvectors, eigenvectors and eigenvalues.
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::Zero(settings.nEigen);
  std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors = nullptr;
  std::shared_ptr<std::vector<Eigen::MatrixXd>> transitiondensities = nullptr;

  std::shared_ptr<std::vector<Eigen::MatrixXd>> solutionvectors = nullptr;
  std::shared_ptr<std::vector<Eigen::MatrixXd>> perturbeddensities = nullptr;
  std::vector<Eigen::Matrix3d> Fdipdip, Fdipmag;

  // Initialize restart system.
  LRSCFRestart<SCFMode> restart(_lrscf, settings, type);

  // Prepare solutions vector in case of a FDEc calculation.
  if (type == Options::LRSCF_TYPE::COUPLED && settings.nEigen > 0 && !settings.restart) {
    eigenvectors = LRSCFSetup<SCFMode>::setupFDEcTransformation(eigenvalues, settings, _couplingPattern,
                                                                _referenceLoadingType, _lrscf, _act, nDimension);
    // Set number of the eigenvalues to be determined in the coupled step.
    settings.nEigen = (*eigenvectors)[0].cols();
  }
  else if (settings.restart) {
    eigenvectors = restart.fetchEigenpairs(eigenvalues);
  }

  // Turn off grid output globally.
  bool oldIOGridInfo = iOOptions.printGridInfo;
  iOOptions.printGridInfo = (GLOBAL_PRINT_LEVEL != Options::GLOBAL_PRINT_LEVELS::DEBUGGING) ? false : true;

  // Setup lambda functions for response matrix sigmavectors.
  ResponseLambda<SCFMode> lambda(_act, _env, _lrscf, diagonal, settings);
  lambda.setupTDDFTLambdas();
  if (!isNotCC2) {
    lambda.setupCC2Lambdas();
  }
  SigmaCalculator sigmaCalculator =
      (settings.method == Options::LR_METHOD::TDDFT) ? lambda.getRPASigma() : lambda.getTDASigma();

  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  if (lambda.usesLRExchange()) {
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 4);
  }

  // Set initial subspace size.
  unsigned initialSubspace = std::min(settings.nEigen * 2, settings.nEigen + 8);
  initialSubspace = (type == Options::LRSCF_TYPE::COUPLED) ? settings.nEigen : initialSubspace;
  initialSubspace = std::min(initialSubspace, nDimension);

  settings.partialResponseConstruction =
      type == Options::LRSCF_TYPE::COUPLED && settings.partialResponseConstruction && !settings.fullFDEc;

  if (!settings.frequencies.empty() || !settings.frequencyRange.empty()) {
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
  }

  // TDDFT.
  if (isNotCC2) {
    std::unique_ptr<EigenvalueSolver> eigensolver;
    std::unique_ptr<ResponseSolver> responsesolver;

    // Eigenvalue Solver.
    if (settings.nEigen > 0) {
      if (settings.partialResponseConstruction) {
        CouplingConstruction<SCFMode>::solve(_lrscf, settings, _referenceLoadingType, sigmaCalculator, eigenvectors, eigenvalues);
      }
      else {
        eigensolver = std::make_unique<EigenvalueSolver>(
            nDimension, settings.nEigen, diagonal, settings.conv,
            (type == Options::LRSCF_TYPE::COUPLED && !settings.fullFDEc) ? 1 : settings.maxCycles, settings.maxSubspaceDimension,
            initialSubspace, settings.method, sigmaCalculator, eigenvectors, restart.getWriteToDisk());
        eigensolver->solve();
      }
    }

    // Response Solver.
    if (!settings.frequencies.empty()) {
      // Determine rhs of the linear equation system.
      std::vector<Eigen::MatrixXd> rhs(settings.damping != 0 ? 4 : 2, Eigen::MatrixXd::Zero(nDimension, 3));
      if (settings.gauge == Options::GAUGE::LENGTH) {
        rhs[0] = 2 * (*dipoles->getLengths());
      }
      else {
        rhs[1] = 2 * (*dipoles->getVelocities());
      }

      responsesolver =
          std::make_unique<ResponseSolver>(diagonal, settings.conv, settings.maxCycles, settings.maxSubspaceDimension,
                                           settings.frequencies, settings.damping, rhs, lambda.getRPASigma());
      responsesolver->solve();
    }

    // Switch to finer grid if necessary and perform one final iteration.
    if (lambda.usesKernel() && settings.grid.accuracy != settings.grid.smallGridAccuracy) {
      printSmallCaption("Switching from small to default grid");
      printf("  Residual norms might not match the desired convergence anymore.\n");
      printf("  Grid accuracy:    %1i -> %1i\n\n", settings.grid.smallGridAccuracy, settings.grid.accuracy);
      lambda.setupKernel(Options::GRID_PURPOSES::DEFAULT);

      // Eigenvalue Solver.
      if (settings.nEigen > 0 && !settings.partialResponseConstruction) {
        printBigCaption("Eigenvalue Solver");
        eigensolver->iterate();
        eigensolver->postProcessing();
      }

      // Response Solver.
      if (!settings.frequencies.empty()) {
        printBigCaption("Response Solver");
        responsesolver->iterate();
      }
    }

    if (settings.nEigen > 0) {
      if (!settings.partialResponseConstruction) {
        eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(eigensolver->getEigenvectors());
        eigenvalues = eigensolver->getEigenvalues();
      }

      // For TDA: (X+Y) == (X-Y).
      if (settings.method == Options::LR_METHOD::TDA) {
        if (eigenvectors->size() == 1) {
          eigenvectors->push_back((*eigenvectors)[0]);
        }
        else {
          (*eigenvectors)[1] = (*eigenvectors)[0];
        }
      }

      // For (TDA-)TDDFT transition density matrices: (X+Y), (X-Y).
      if (settings.scfstab == Options::STABILITY_ANALYSIS::NONE) {
        transitiondensities = eigenvectors;
      }

      // If double hybrid is used, calculate CIS(D) correction.
      if (lambda.usesDoubleHybrid()) {
        lambda.setupCC2Lambdas();
        Eigen::VectorXd cisd = eigenvalues;

        NonlinearEigenvalueSolver nlEigenSolver(settings.nEigen, diagonal, settings.conv, settings.preopt, settings.diis,
                                                settings.diisStore, settings.maxCycles, false, Options::LR_METHOD::CISD,
                                                (*eigenvectors)[0], cisd, lambda.getRightCC2Sigma());
        nlEigenSolver.solve();

        // Add correction to excitation energies.
        eigenvalues += lambda.getDHRatio() * cisd;

        // Get rid of CC2 and RIIntegral controller.
        for (auto& lrscf : _lrscf) {
          lrscf->finalizeCC2Controller();
          lrscf->finalizeRIIntegrals(LIBINT_OPERATOR::coulomb);
        }
      } /* Double hybrid correction */
    }

    if (!settings.frequencies.empty()) {
      solutionvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(responsesolver->getEigenvectors());
      perturbeddensities = solutionvectors;
    }
  }

  // CC2.
  else {
    // Eigenvalue Solver.
    if (settings.nEigen > 0) {
      if (!eigenvectors || settings.method == Options::LR_METHOD::CISD) {
        CC2HelperFunctions<SCFMode>::calculateCISEigenvectors(settings, eigenvectors, eigenvalues, nDimension, diagonal,
                                                              initialSubspace, sigmaCalculator);
      }

      Eigen::VectorXd initEigenvalues = eigenvalues;
      CC2HelperFunctions<SCFMode>::calculateRightEigenvectors(settings, eigenvectors, eigenvalues, diagonal,
                                                              lambda.getRightCC2Sigma(), type, restart.getWriteToDisk());

      if (settings.method == Options::LR_METHOD::CISD) {
        eigenvalues += initEigenvalues;
        eigenvectors->push_back((*eigenvectors)[0]);
        transitiondensities = std::make_shared<std::vector<Eigen::MatrixXd>>(2, (*eigenvectors)[0]);
      }
      else {
        if (eigenvectors->size() == 1) {
          eigenvectors->push_back((*eigenvectors)[0]);
        }

        if (!(settings.cctrdens || settings.ccexdens) || settings.triplet || settings.method == Options::LR_METHOD::ADC2) {
          // In an FDEc calculation, make sure that [1] are not still the FDEu vectors in the ADC(2) case.
          (*eigenvectors)[1] = (*eigenvectors)[0];
          CC2HelperFunctions<SCFMode>::normalizeRightEigenvectors(_lrscf, settings, eigenvectors, eigenvalues);
        }
        else {
          // Trick: in an FDEc calculation, initEigenvalues holds the FDEu eigenvalues, this is why
          // we plug them in here to be consistent with the right FDEc eigenvector procedure.
          if (type == Options::LRSCF_TYPE::COUPLED) {
            eigenvalues = initEigenvalues;
          }
          CC2HelperFunctions<SCFMode>::calculateLeftEigenvectors(settings, eigenvectors, eigenvalues, diagonal,
                                                                 lambda.getLeftCC2Sigma(), type);
          CC2HelperFunctions<SCFMode>::normalizeBothEigenvectors(_lrscf, settings, eigenvectors, eigenvalues);
        }

        if ((settings.cctrdens || settings.ccexdens) && !settings.triplet) {
          CC2HelperFunctions<SCFMode>::prepareVectors(_lrscf, settings, eigenvectors, transitiondensities);
          if (settings.cctrdens || settings.ccexdens) {
            if (settings.method == Options::LR_METHOD::CC2) {
              CC2HelperFunctions<SCFMode>::calculateStateMultipliers(_lrscf, settings, eigenvectors, eigenvalues,
                                                                     lambda.getLeftCC2Sigma(), diagonal);
            }
            CC2HelperFunctions<SCFMode>::calculateStateDensities(_lrscf, transitiondensities, eigenvectors, eigenvalues);
          }
          if (settings.cctrdens) {
            if (settings.method == Options::LR_METHOD::CC2) {
              CC2HelperFunctions<SCFMode>::calculateTransitionMultipliers(_lrscf, settings, eigenvectors, eigenvalues,
                                                                          lambda.getLeftCC2Sigma(), diagonal);
            }
            CC2HelperFunctions<SCFMode>::calculateTransitionDensities(_lrscf, transitiondensities, eigenvectors, eigenvalues);
          }
        }
      }
    }
    // AR: the transitiondensities-vector contains three matrices: first the right transition density, then the left
    // transition density, and as the third element the state densities. in each matrix, the columns represent the
    // different states (CC2 case: zeroth column of the third matrix is the ground state density), and there are (n_MO *
    // n_MO) rows (see prepareVectors function of CC2HelperFunctions)

    // Response Solver.
    if (!settings.frequencies.empty()) {
      if (settings.method != Options::LR_METHOD::CC2) {
        WarningTracker::printWarning(
            "Response Properties not properly defined and implemented for second-order methods other than CC2.", true);
      }

      unsigned nFreqs = settings.frequencies.size();

      // Check if there is a static frequency anywhere for index book-keeping.
      for (unsigned iFreq = 0; iFreq < settings.frequencies.size(); ++iFreq) {
        if (settings.frequencies[iFreq] == 0) {
          // Move static frequency to back.
          auto it = settings.frequencies.begin() + iFreq;
          std::rotate(it, it + 1, settings.frequencies.end());
        }
      }

      // Create vector containing the negative frequencies.
      auto frequencies = settings.frequencies;
      for (unsigned iFreq = 0; iFreq < nFreqs; ++iFreq) {
        // Don't append zero frequency.
        if (settings.frequencies[iFreq]) {
          frequencies.push_back(-settings.frequencies[iFreq]);
        }
      }

      // Ground state multipliers and densities not yet calculated? Do so now.
      if (settings.nEigen == 0 || !(settings.cctrdens || settings.ccexdens)) {
        CC2HelperFunctions<SCFMode>::prepareVectors(_lrscf, settings, eigenvectors, transitiondensities);
        if (settings.method == Options::LR_METHOD::CC2) {
          CC2HelperFunctions<SCFMode>::calculateStateMultipliers(_lrscf, settings, eigenvectors, eigenvalues,
                                                                 lambda.getLeftCC2Sigma(), diagonal);
          CC2HelperFunctions<SCFMode>::calculateStateDensities(_lrscf, transitiondensities, eigenvectors, eigenvalues);
        }
      }

      Eigen::MatrixXd dips(dipoles->getLengths()->rows(), 6);
      dips.leftCols(3) = (settings.gauge == Options::GAUGE::LENGTH) ? (*dipoles->getLengths()) : (*dipoles->getVelocities());
      dips.rightCols(3) = (*dipoles->getMagnetics());
      CC2HelperFunctions<SCFMode>::calculatePerturbedAmplitudes(_lrscf, settings, solutionvectors, frequencies,
                                                                lambda.getRightCC2Sigma(), diagonal, dips);
      CC2HelperFunctions<SCFMode>::calculatePerturbedDensities(_lrscf, settings, solutionvectors, perturbeddensities,
                                                               frequencies, dips, Fdipdip, Fdipmag);

      // Do not assume that something useful can be found under the eigenvectors pointer.
      if (settings.nEigen == 0) {
        eigenvectors = nullptr;
      }
    }
  }

  if (eigenvectors) {
    restart.storeConvergedSolution(eigenvectors, eigenvalues);
  }
  if (solutionvectors && isNotCC2) {
    restart.storeConvergedResponse(solutionvectors,
                                   Eigen::Map<Eigen::VectorXd>(settings.frequencies.data(), settings.frequencies.size()));
  }

  if (settings.excGradList.size()) {
    if (!((settings.method == Options::LR_METHOD::TDA) || (settings.method == Options::LR_METHOD::TDDFT)) ||
        (type != Options::LRSCF_TYPE::ISOLATED)) {
      throw SerenityError(
          "Excited-state gradients are only implemented for isolated TDA/TDDFT/CIS/TDHF calculations so far!");
    }
    for (unsigned iExc : settings.excGradList) {
      if (!iExc)
        throw SerenityError(
            "Excgradlist in the LRSCFTask is 1-based, i.e. the lowest-lying excited-state is number 1, not 0.");
    }
    TDDFTGradientCalculator<SCFMode> tddftgrads(_lrscf[0], settings.hypthresh);
    tddftgrads.calculateGradients();
  }

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
      _excitations.col(0) = eigenvalues;
      // Prints dominant contributions etc.
      LRSCFAnalysis<SCFMode>::printDominantContributions(_lrscf, (*eigenvectors), eigenvalues, settings.dominantThresh);
      // Print delta s^2 values in unrestricted case.
      if (isNotCC2 && SCFMode == Options::SCF_MODES::UNRESTRICTED) {
        DeltaSpinSquared<SCFMode> spinSquared(_lrscf, settings.nEigen, settings.method, type);
        spinSquared.print();
      }
      // No excitation spectrum for triplet excitations or stability analyses.
      if (transitiondensities && !(settings.triplet || settings.scfstab != Options::STABILITY_ANALYSIS::NONE)) {
        // Prints excitation and CD spectrum (oscillator and rotatory strengths, respectively).
        std::string name = (type != Options::LRSCF_TYPE::COUPLED)
                               ? _lrscf[0]->getSys()->getSystemPath() + _lrscf[0]->getSys()->getSystemName()
                               : "";
        if (isNotCC2 || (!isNotCC2 && settings.cctrdens) || settings.method == Options::LR_METHOD::CISD) {
          ExcitationSpectrum<SCFMode>::printTransitionMoments(settings.method, dipoles, (*transitiondensities),
                                                              eigenvalues, _excitations, name);
          // _lrscf only contains more than one LRSCFController in the coupled case
          if (isNotCC2 && (type == Options::LRSCF_TYPE::ISOLATED))
            ExcitationSpectrum<SCFMode>::printStateMoments(settings.method, dipoles, *_lrscf[0]->getUnrelaxedDiffDensities(),
                                                           *_lrscf[0]->getRelaxedDiffDensities(), eigenvalues,
                                                           nucDipoleMoment, _lrscf[0]->getNOccupied());
        }
        if (!isNotCC2 && settings.ccexdens && settings.method != Options::LR_METHOD::CISD) {
          ExcitationSpectrum<SCFMode>::printStateMoments(settings.method, dipoles, (*transitiondensities), eigenvalues,
                                                         _excitations, nucDipoleMoment);
        }
      }
      // Calculate Löwdin charges and store on disk.
      if (settings.transitionCharges) {
        if (settings.method == Options::LR_METHOD::CC2) {
          WarningTracker::printWarning("CC2 transition densities might not be well defined.", true);
        }
        LRSCFPopulationAnalysis<SCFMode>::calculateTransitionCharges(_lrscf, (*transitiondensities));
      }
    }
    if (solutionvectors) {
      // See Grimme Chem. Phys. Lett. 361 (2002) 321-328.
      double rotFactor = 1.343e-4 / molWeight * HARTREE_TO_OOCM * HARTREE_TO_OOCM;

      // Print properties.
      ResponseProperties<SCFMode>::printProperties(isNotCC2, dipoles, (*perturbeddensities), settings.frequencies, rotFactor,
                                                   settings.damping, settings.gauge, Fdipdip, Fdipmag, _properties);
    }
  }
  iOOptions.printGridInfo = oldIOGridInfo;

  // In case of stability anaysis: Rotate orbitals according to found instabilities.
  if (settings.stabroot != 0) {
    for (auto& lrscf : _lrscf) {
      lrscf->rotateOrbitalsSCFInstability();
    }
  }
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
