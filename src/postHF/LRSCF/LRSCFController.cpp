/**
 * @file LRSCFController.cpp
 *
 * @date Dec 06, 2018
 * @author Michael Boeckers, Niklas Niemeyer, Johannes Toelle
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
#include "postHF/LRSCF/LRSCFController.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Atom.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutputStream.h"
#include "io/HDF5.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "misc/WarningTracker.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/RICC2/ADC2Controller.h"
#include "postHF/LRSCF/RICC2/CC2Controller.h"
#include "postHF/LRSCF/Tools/Besley.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "postHF/MBPT/MBPT.h"
#include "potentials/bundles/PotentialBundle.h"
#include "settings/LRSCFOptions.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <numeric>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LRSCFController<SCFMode>::LRSCFController(std::shared_ptr<SystemController> system, const LRSCFTaskSettings& settings)
  : _system(system),
    _settings(settings),
    _nOcc(_system->getNOccupiedOrbitals<SCFMode>()),
    _nVirt(_system->getNVirtualOrbitalsTruncated<SCFMode>()),
    _coefficients(_system->getActiveOrbitalController<SCFMode>()->getCoefficients()),
    _particleCoefficients(_system->getActiveOrbitalController<SCFMode>()->getCoefficients()),
    _holeCoefficients(_system->getActiveOrbitalController<SCFMode>()->getCoefficients()),
    _orbitalEnergies(_system->getActiveOrbitalController<SCFMode>()->getEigenvalues()),
    _excitationVectors(nullptr),
    _CC2Controller(nullptr),
    _riints(nullptr),
    _riErfints(nullptr),
    _screening(nullptr),
    _envTransformation(nullptr),
    _inverseM(nullptr),
    _inverseErfM(nullptr) {
  Settings sysSettings = _system->getSettings();
  if (sysSettings.load != "") {
    _system->template getElectronicStructure<SCFMode>()->toHDF5(sysSettings.path + sysSettings.name, sysSettings.identifier);
  }
  // AR: the member variable _settings is declared as const in the LRSCFController and therefore cannot be given to
  // Options::resolve because resolve goes both ways and takes references
  LRSCFTaskSettings settingsCopyToPrint = _settings;
  settingsCopyToPrint.printSettings(this->_system->getSystemPath() + this->_system->getSystemName() + "_lrscf.settings");
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<std::vector<Eigen::MatrixXd>> LRSCFController<SCFMode>::getExcitationVectors(Options::LRSCF_TYPE type) {
  if (!_excitationVectors || _type != type) {
    this->loadFromH5(type);
  }
  return _excitationVectors;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Eigen::VectorXd> LRSCFController<SCFMode>::getExcitationEnergies(Options::LRSCF_TYPE type) {
  if (!_excitationEnergies || _type != type) {
    this->loadFromH5(type);
  }
  return _excitationEnergies;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::loadFromH5(Options::LRSCF_TYPE type) {
  auto settings = _system->getSettings();

  std::string tName;
  if (type == Options::LRSCF_TYPE::ISOLATED) {
    tName += "iso";
  }
  else if (type == Options::LRSCF_TYPE::UNCOUPLED) {
    tName += "fdeu";
  }
  else {
    tName += "fdec";
  }

  printSmallCaption("Loading " + tName + "-eigenpairs for: " + settings.name);

  std::string path = (settings.load == "") ? settings.path : settings.load;
  std::string fName = settings.name + "_lrscf.";

  std::vector<Eigen::MatrixXd> XY(2);
  Eigen::VectorXd eigenvalues;

  auto loadEigenpairs = [&](std::string input) {
    OutputControl::n.printf("\n   $  %-20s\n\n", input.c_str());

    HDF5::H5File file(input.c_str(), H5F_ACC_RDONLY);
    HDF5::dataset_exists(file, "RIGHT");
    HDF5::dataset_exists(file, "LEFT");
    HDF5::dataset_exists(file, "EIGENVALUES");
    HDF5::load(file, "RIGHT", XY[0]);
    HDF5::load(file, "LEFT", XY[1]);
    HDF5::load(file, "EIGENVALUES", eigenvalues);
    file.close();

    unsigned nDimI = 0;
    for_spin(_nOcc, _nVirt) {
      nDimI += _nOcc_spin * _nVirt_spin;
    };
    if (XY[0].rows() != nDimI || XY[1].rows() != nDimI) {
      throw SerenityError("The dimension of your loaded eigenpairs does not match with this response problem.");
    }
    if (eigenvalues.rows() != XY[0].cols()) {
      throw SerenityError("The number of loaded eigenvectors and eigenvalues does not match.");
    }

    OutputControl::n.printf("  Found %3i eigenpairs.\n\n\n", (int)eigenvalues.rows());

    _excitationEnergies = std::make_shared<Eigen::VectorXd>(eigenvalues);
    _excitationVectors = std::make_shared<std::vector<Eigen::MatrixXd>>(XY);
    _type = type;
  };

  fName += tName + (SCFMode == RESTRICTED ? ".res.h5" : ".unres.h5");

  try {
    loadEigenpairs(path + fName);
  }
  catch (...) {
    if (settings.load != "") {
      OutputControl::n.printf("  Could not find any. Instead\n");
      try {
        loadEigenpairs(settings.path + fName);
      }
      catch (...) {
        throw SerenityError("Cannot find eigenpairs from h5: " + (settings.path + fName));
      }
    }
    else {
      throw SerenityError("Cannot find eigenpairs from h5: " + (path + fName));
    }
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setSolution(std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors,
                                           std::shared_ptr<Eigen::VectorXd> eigenvalues, Options::LRSCF_TYPE type) {
  _excitationVectors = eigenvectors;
  _excitationEnergies = eigenvalues;
  _type = type;

  std::string fName = this->getSys()->getSystemPath() + this->getSys()->getSystemName() + "_lrscf.";
  if (type == Options::LRSCF_TYPE::ISOLATED) {
    fName += "iso.";
  }
  else if (type == Options::LRSCF_TYPE::UNCOUPLED) {
    fName += "fdeu.";
  }
  else {
    fName += "fdec.";
  }
  fName += (SCFMode == RESTRICTED) ? "res.h5" : "unres.h5";

  HDF5::H5File file(fName.c_str(), H5F_ACC_TRUNC);
  HDF5::save_scalar_attribute(file, "ID", _system->getSystemIdentifier());
  HDF5::save(file, "RIGHT", (*eigenvectors)[0]);
  HDF5::save(file, "LEFT", (*eigenvectors)[1]);
  HDF5::save(file, "EIGENVALUES", (*eigenvalues));
  file.close();

  for (unsigned i = 0; i < (*eigenvectors)[0].cols(); i++) {
    MatrixInBasis<SCFMode> T(this->getBasisController());

    unsigned iaStart = 0;
    for_spin(T, _nOcc, _nVirt) {
      Eigen::Map<Eigen::MatrixXd> right =
          Eigen::Map<Eigen::MatrixXd>((*eigenvectors)[0].col(i).data() + iaStart, _nVirt_spin, _nOcc_spin);
      Eigen::Map<Eigen::MatrixXd> left =
          Eigen::Map<Eigen::MatrixXd>((*eigenvectors)[1].col(i).data() + iaStart, _nVirt_spin, _nOcc_spin);
      T_spin.bottomRightCorner(_nVirt_spin, _nVirt_spin) = 0.5 * (right * right.transpose() + left * left.transpose());
      T_spin.topLeftCorner(_nOcc_spin, _nOcc_spin) = -0.5 * (right.transpose() * right + left.transpose() * left);
      iaStart += _nOcc_spin * _nVirt_spin;
    };
    _unrelaxedDiffDensities->push_back(T);
  }
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, unsigned> LRSCFController<SCFMode>::getNOccupied() {
  return _nOcc;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, unsigned> LRSCFController<SCFMode>::getNVirtual() {
  return _nVirt;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setNOccupied(SpinPolarizedData<SCFMode, unsigned> nOcc) {
  _nOcc = nOcc;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setNVirtual(SpinPolarizedData<SCFMode, unsigned> nVirt) {
  _nVirt = nVirt;
}

template<Options::SCF_MODES SCFMode>
CoefficientMatrix<SCFMode>& LRSCFController<SCFMode>::getCoefficients() {
  return _coefficients;
}

template<Options::SCF_MODES SCFMode>
CoefficientMatrix<SCFMode>& LRSCFController<SCFMode>::getParticleCoefficients() {
  if (_CC2Controller) {
    auto& P = _CC2Controller->_P;
    for_spin(P, _particleCoefficients) {
      _particleCoefficients_spin = P_spin;
    };
  }
  else {
    _particleCoefficients = _coefficients;
  }
  return _particleCoefficients;
}

template<Options::SCF_MODES SCFMode>
CoefficientMatrix<SCFMode>& LRSCFController<SCFMode>::getHoleCoefficients() {
  if (_CC2Controller) {
    auto& H = _CC2Controller->_H;
    for_spin(H, _holeCoefficients) {
      _holeCoefficients_spin = H_spin;
    };
  }
  else {
    _holeCoefficients = _coefficients;
  }
  return _holeCoefficients;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setCoefficients(CoefficientMatrix<SCFMode> coeff) {
  _coefficients = coeff;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<BasisController> LRSCFController<SCFMode>::getBasisController(Options::BASIS_PURPOSES basisPurpose) {
  return _system->getBasisController(basisPurpose);
}

template<>
std::shared_ptr<SPMatrix<RESTRICTED>> LRSCFController<RESTRICTED>::getMOFockMatrix() {
  if (!_system->template getElectronicStructure<RESTRICTED>()->checkFock()) {
    throw SerenityError("LRSCF: No Fock matrix present in your system.");
  }
  auto fock = _system->template getElectronicStructure<RESTRICTED>()->getFockMatrix();
  if (!_settings.rpaScreening) {
    unsigned nb = _nOcc + _nVirt;
    fock = (_coefficients.leftCols(nb).transpose() * fock * _coefficients.leftCols(nb)).eval();
  }
  else {
    Eigen::MatrixXd temp = _orbitalEnergies.segment(0, _nOcc + _nVirt).asDiagonal();
    fock = temp;
  }
  _fock = std::make_shared<SPMatrix<RESTRICTED>>(fock);
  return _fock;
}

template<>
std::shared_ptr<SPMatrix<UNRESTRICTED>> LRSCFController<UNRESTRICTED>::getMOFockMatrix() {
  if (!_system->template getElectronicStructure<UNRESTRICTED>()->checkFock()) {
    throw SerenityError("lrscf: no fock matrix present in your system.");
  }
  auto fock = _system->template getElectronicStructure<UNRESTRICTED>()->getFockMatrix();
  if (!_settings.rpaScreening) {
    if (_settings.scfstab != Options::STABILITY_ANALYSIS::SPINFLIP) {
      for_spin(_coefficients, fock, _nOcc, _nVirt) {
        unsigned nb = _nOcc_spin + _nVirt_spin;
        fock_spin = (_coefficients_spin.leftCols(nb).transpose() * fock_spin * _coefficients_spin.leftCols(nb)).eval();
      };
    }
    else {
      // For Spin-Flip TDDFT, the reference orbitals are occupied (spin excess) and the respective other spin as the
      // virtuals. Therefore, only one orbital transition space is needed which is the alpha part of everything here. To
      // have a correct MO representation of the Fock matrix, the constructed alpha part here must contain the alpha
      // occ-occ block and the beta virt-virt block. This needs to be taken account of. Check out the
      // setupSpinFlipReference() member function of this class to understand this better.
      unsigned nmos = _nOcc.alpha + _nVirt.alpha;
      Eigen::MatrixXd F = Eigen::MatrixXd::Zero(nmos, nmos);
      F.topLeftCorner(_nOcc.alpha, _nOcc.alpha) = _coefficients.alpha.middleCols(0, _nOcc.alpha).transpose() *
                                                  (_system->getSpin() > 0 ? fock.alpha : fock.beta) *
                                                  _coefficients.alpha.middleCols(0, _nOcc.alpha);
      F.bottomRightCorner(_nVirt.alpha, _nVirt.alpha) =
          _coefficients.alpha.middleCols(_nOcc.alpha, _nVirt.alpha).transpose() *
          (_system->getSpin() > 0 ? fock.beta : fock.alpha) * _coefficients.alpha.middleCols(_nOcc.alpha, _nVirt.alpha);
      fock.alpha = F;
      fock.beta = F;
    }
  }
  else {
    for_spin(fock, _orbitalEnergies, _nOcc, _nVirt) {
      Eigen::MatrixXd temp = _orbitalEnergies_spin.segment(0, _nOcc_spin + _nVirt_spin).asDiagonal();
      fock_spin = temp;
    };
  }
  _fock = std::make_shared<SPMatrix<UNRESTRICTED>>(fock);
  return _fock;
}

template<Options::SCF_MODES SCFMode>
bool LRSCFController<SCFMode>::isMOFockMatrixDiagonal() {
  SPMatrix<SCFMode> fock = *this->getMOFockMatrix();
  bool isDiagonal = true;
  for_spin(fock) {
    fock_spin.diagonal() -= fock_spin.diagonal();
    isDiagonal = isDiagonal && std::fabs(fock_spin.sum()) < 1e-6;
    if (std::fabs(fock_spin.sum()) > 1e-6)
      OutputControl::nOut << " Absolute largest Fock matrix off-Diagonal element  "
                          << fock_spin.array().abs().matrix().maxCoeff() << std::endl;
  };

  return isDiagonal;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd> LRSCFController<SCFMode>::getEigenvalues() {
  return _orbitalEnergies;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setEigenvalues(SpinPolarizedData<SCFMode, Eigen::VectorXd> eigenvalues) {
  _orbitalEnergies = eigenvalues;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<GridController> LRSCFController<SCFMode>::getGridController() {
  return _system->getGridController();
}

template<Options::SCF_MODES SCFMode>
const Settings& LRSCFController<SCFMode>::getSysSettings() {
  return _system->getSettings();
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<SystemController> LRSCFController<SCFMode>::getSys() {
  return _system;
}

template<Options::SCF_MODES SCFMode>
const LRSCFTaskSettings& LRSCFController<SCFMode>::getLRSCFSettings() {
  return _settings;
}

template<Options::SCF_MODES SCFMode>
Options::LR_METHOD LRSCFController<SCFMode>::getResponseMethod() {
  return _settings.method;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setEnvSystems(std::vector<std::shared_ptr<SystemController>> envSystems) {
  _envSystems = envSystems;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<SystemController>> LRSCFController<SCFMode>::getEnvSystems() {
  return _envSystems;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::initializeCC2Controller() {
  if (this->getResponseMethod() == Options::LR_METHOD::ADC2) {
    _CC2Controller = std::make_shared<ADC2Controller<SCFMode>>(this->shared_from_this());
  }
  else {
    _CC2Controller = std::make_shared<CC2Controller<SCFMode>>(this->shared_from_this());
  }
  _CC2Controller->initialize();
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<CC2Controller<SCFMode>> LRSCFController<SCFMode>::getCC2Controller() {
  return _CC2Controller;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::finalizeCC2Controller() {
  _CC2Controller = nullptr;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::initializeRIIntegrals(LIBINT_OPERATOR op, double mu, bool calcJia) {
  if (op == LIBINT_OPERATOR::coulomb) {
    _riints = std::make_shared<RIIntegrals<SCFMode>>(this->shared_from_this(), op, mu, calcJia);
  }
  else if (op == LIBINT_OPERATOR::erf_coulomb) {
    _riErfints = std::make_shared<RIIntegrals<SCFMode>>(this->shared_from_this(), op, mu, calcJia);
  }
  else {
    throw SerenityError("This operator for RI integrals is not yet supported.");
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::finalizeRIIntegrals(LIBINT_OPERATOR op) {
  if (op == LIBINT_OPERATOR::coulomb) {
    _riints = nullptr;
  }
  else if (op == LIBINT_OPERATOR::erf_coulomb) {
    _riErfints = nullptr;
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<RIIntegrals<SCFMode>> LRSCFController<SCFMode>::getRIIntegrals(LIBINT_OPERATOR op) {
  if (op == LIBINT_OPERATOR::coulomb) {
    return _riints;
  }
  else if (op == LIBINT_OPERATOR::erf_coulomb) {
    return _riErfints;
  }
  else {
    throw SerenityError("This operator for RI integrals is not yet supported.");
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::calculateScreening(const Eigen::VectorXd& eia) {
  if (_settings.nafThresh != 0) {
    throw SerenityError(" NAF not supported with environmetnal screening!");
  }
  if (_riints != nullptr) {
    printBigCaption("RPA Screening");
    auto& jia = *(_riints->getJiaPtr());
    unsigned nAux = _riints->getNTransformedAuxBasisFunctions();
    double prefactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
    Eigen::MatrixXd unit = Eigen::MatrixXd::Identity(nAux, nAux);
    Eigen::MatrixXd piPQ = unit;
    Eigen::VectorXd chi_temp = prefactor * (-2.0 / eia.array());

    if (_envSystems.size() == 0) {
      unsigned spincounter = 0;
      for_spin(jia) {
        piPQ -= jia_spin.transpose() * chi_temp.segment(spincounter, jia_spin.rows()).asDiagonal() * jia_spin;
        spincounter += jia_spin.rows();
      };
      _screening = std::make_shared<Eigen::MatrixXd>(piPQ.householderQr().solve(unit));
    }
    if (_envSystems.size() > 0) {
      if (_settings.nafThresh != 0) {
        throw SerenityError(" NAF not supported with environmetnal screening!");
      }
      GWTaskSettings gwSettings;
      gwSettings.environmentScreening = false;
      gwSettings.integrationPoints = 0;
      // sets the geometry of the active subsystem to evaluate integrals in the correct aux basis
      _riints->setGeo(_system->getGeometry());
      auto mbpt = std::make_shared<MBPT<SCFMode>>(this->shared_from_this(), gwSettings, _envSystems, _riints, 0, 0);
      auto envResponse = mbpt->environmentRespose();
      auto auxtrafo = _riints->getAuxTrafoPtr();
      auto metric = _riints->getAuxMetricPtr();
      auto metricsqrt = mSqrt_Sym(*metric);
      Eigen::MatrixXd coulAuxbasis = (*metric) + metricsqrt * envResponse * metricsqrt;

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(coulAuxbasis);
      Eigen::VectorXd eValues = eigenvalueSolver.eigenvalues();
      Eigen::MatrixXd eigenVectors = eigenvalueSolver.eigenvectors();

      unsigned int spincounter = 0;
      for_spin(jia) {
        auto piPQenv = (jia_spin * (*auxtrafo) * eigenVectors).transpose() *
                       chi_temp.segment(spincounter, jia_spin.rows()).asDiagonal() * (jia_spin * (*auxtrafo) * eigenVectors);
        piPQ = unit - piPQenv * eValues.asDiagonal();
        spincounter += jia_spin.rows();
      };
      Eigen::MatrixXd screen = piPQ.householderQr().solve(unit);
      screen =
          ((*auxtrafo) * eigenVectors * eValues.asDiagonal() * screen * eigenVectors.transpose() * (*auxtrafo).transpose())
              .eval();
      _riints = nullptr;
      _screening = std::make_shared<Eigen::MatrixXd>(screen);
    }
    OutputControl::n.printf(" .. done.\n\n");
  }
  else {
    throw SerenityError("No RI integrals for screening initialized!");
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Eigen::MatrixXd> LRSCFController<SCFMode>::getScreeningAuxMatrix() {
  return _screening;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Eigen::MatrixXd> LRSCFController<SCFMode>::getEnvTrafo() {
  return _envTransformation;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setInverseMetric(std::shared_ptr<Eigen::MatrixXd> metric) {
  _inverseM = metric;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Eigen::MatrixXd> LRSCFController<SCFMode>::getInverseMetric() {
  return _inverseM;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setInverseErfMetric(std::shared_ptr<Eigen::MatrixXd> metric) {
  _inverseErfM = metric;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Eigen::MatrixXd> LRSCFController<SCFMode>::getInverseErfMetric() {
  return _inverseErfM;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::editReference(SpinPolarizedData<SCFMode, std::vector<unsigned>> indexWhiteList) {
  unsigned iSpin = 0;
  for_spin(_coefficients, _orbitalEnergies, indexWhiteList, _nOcc, _nVirt) {
    // Update orbitals
    auto oldCoefficients = _coefficients_spin;
    auto oldOrbitalEnergies = _orbitalEnergies_spin;
    unsigned oldnOcc = _nOcc_spin;

    _coefficients_spin.setZero();
    _orbitalEnergies_spin.setZero();
    _nOcc_spin = 0;
    _nVirt_spin = 0;

    for (unsigned iMO = 0; iMO < indexWhiteList_spin.size(); ++iMO) {
      _coefficients_spin.col(iMO) = oldCoefficients.col(indexWhiteList_spin[iMO]);
      _orbitalEnergies_spin(iMO) = oldOrbitalEnergies(indexWhiteList_spin[iMO]);
      if (indexWhiteList_spin[iMO] < oldnOcc) {
        _nOcc_spin += 1;
      }
      else {
        _nVirt_spin += 1;
      }
    } /* Update orbitals */
    std::string system = _system->getSettings().name;
    // Print info
    if (SCFMode == Options::SCF_MODES::RESTRICTED) {
      OutputControl::nOut << " System: " << system << " \n";
      OutputControl::n.printf(" New reference orbitals (occ | virt): \n");
    }
    else {
      OutputControl::nOut << " System: " << system << " \n";
      OutputControl::n.printf("%s New reference orbitals (occ | virt): \n", (iSpin == 0) ? "Alpha" : "Beta");
    }

    bool printedOccVirtSep = false;
    for (unsigned iMO = 0; iMO < indexWhiteList_spin.size(); ++iMO) {
      if (indexWhiteList_spin[iMO] >= oldnOcc && !printedOccVirtSep) {
        OutputControl::n.printf("   |");
        printedOccVirtSep = true;
      }
      OutputControl::n.printf("%4i", indexWhiteList_spin[iMO] + 1);
      if ((iMO + 1) % 10 == 0) {
        OutputControl::n.printf("\n");
      }
    }
    OutputControl::n.printf("\n");
    iSpin += 1;
  };
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::applyFrozenCore() {
  if (_settings.coreOnly) {
    throw SerenityError("LRSCFTask: frozenCore and coreOnly are mutually exclusive!");
  }

  unsigned nCore = 0;
  for (const auto& atom : _system->getGeometry()->getAtoms()) {
    nCore += atom->getNCoreElectrons();
  }

  nCore /= 2;

  bool isAlpha = true;
  for_spin(_nOcc, _nVirt, _coefficients, _orbitalEnergies) {
    if (isAlpha) {
      OutputControl::nOut << "  Freezing alpha core orbitals: " << nCore << ",   " << 100.0 * nCore / _nOcc_spin
                          << " %." << std::endl;
    }
    else {
      OutputControl::nOut << "  Freezing beta core orbitals: " << nCore << ",   " << 100.0 * nCore / _nOcc_spin << " %."
                          << std::endl;
    }
    _coefficients_spin.middleCols(0, _nOcc_spin - nCore + _nVirt_spin) =
        _coefficients_spin.middleCols(nCore, _nOcc_spin - nCore + _nVirt_spin).eval();
    _orbitalEnergies_spin.middleRows(0, _nOcc_spin - nCore + _nVirt_spin) =
        _orbitalEnergies_spin.middleRows(nCore, _nOcc_spin - nCore + _nVirt_spin).eval();

    _coefficients_spin.middleCols(_nOcc_spin - nCore + _nVirt_spin, nCore).setZero();
    _orbitalEnergies_spin.middleRows(_nOcc_spin - nCore + _nVirt_spin, nCore).setZero();

    _nOcc_spin = _nOcc_spin - nCore;
    isAlpha = false;
  };
  OutputControl::nOut << std::endl;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::applyFrozenVirtual() {
  if (_settings.frozenVirtual < 0) {
    throw SerenityError("LRSCFTask: frozenVirtual keyword must be positive!");
  }
  OutputControl::nOut << "  Freezing all virtual orbitals beyond e(HOMO)+" << _settings.frozenVirtual << " eV." << std::endl;

  bool isAlpha = true;
  for_spin(_nOcc, _nVirt, _coefficients, _orbitalEnergies) {
    double e_homo = _orbitalEnergies_spin[_nOcc_spin - 1];
    for (unsigned a = 0; a < _nVirt_spin; ++a) {
      double e_virt = _orbitalEnergies_spin[_nOcc_spin + a];
      if (e_virt > e_homo + _settings.frozenVirtual * EV_TO_HARTREE) {
        if (isAlpha) {
          OutputControl::nOut << "  Freezing alpha virtual orbitals: " << _nVirt_spin - a << ",   "
                              << 100.0 - 100.0 * a / _nVirt_spin << " %." << std::endl;
        }
        else {
          OutputControl::nOut << "  Freezing beta virtual orbitals: " << _nVirt_spin - a << ",   "
                              << 100.0 - 100.0 * a / _nVirt_spin << " %." << std::endl;
        }
        _coefficients_spin.middleCols(_nOcc_spin + a, _nVirt_spin - a).setZero();
        _orbitalEnergies_spin.middleRows(_nOcc_spin + a, _nVirt_spin - a).setZero();
        _nVirt_spin = a;
        break;
      }
    }
    isAlpha = false;
  };
  OutputControl::nOut << std::endl;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::applyCoreOnly() {
  if (_settings.frozenCore) {
    throw SerenityError("LRSCFTask: frozenCore and coreOnly are mutually exclusive!");
  }

  unsigned nCore = 0;
  for (const auto& atom : _system->getGeometry()->getAtoms()) {
    nCore += atom->getNCoreElectrons();
  }

  nCore /= 2;
  OutputControl::nOut << "  Freezing all but core orbitals:  " << nCore << std::endl;

  for_spin(_nOcc, _nVirt, _coefficients, _orbitalEnergies) {
    _coefficients_spin.middleCols(nCore, _nVirt_spin) = _coefficients_spin.middleCols(_nOcc_spin, _nVirt_spin).eval();
    _orbitalEnergies_spin.middleRows(nCore, _nVirt_spin) = _orbitalEnergies_spin.middleRows(_nOcc_spin, _nVirt_spin).eval();

    _coefficients_spin.middleCols(nCore + _nVirt_spin, _nOcc_spin - nCore).setZero();
    _orbitalEnergies_spin.middleRows(nCore + _nVirt_spin, _nOcc_spin - nCore).setZero();

    _nOcc_spin = nCore;
  };
  OutputControl::nOut << std::endl;
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::setupSpinFlipReference() {
  if (SCFMode == Options::SCF_MODES::RESTRICTED || !(_system->getSpin() != 2 || _system->getSpin() != -2)) {
    throw SerenityError("Spin-Flip TDDFT/TDHF only possible for unrestricted triplet reference!");
  }

  unsigned noa = 0, nva = 0;
  unsigned nob = 0, nvb = 0;

  Eigen::MatrixXd ca, cb;
  Eigen::VectorXd ea, eb;

  unsigned iSpin = 0;
  for_spin(_nOcc, _nVirt, _coefficients, _orbitalEnergies) {
    if (iSpin == 0) {
      noa = _nOcc_spin;
      nva = _nVirt_spin;
      ca = _coefficients_spin;
      ea = _orbitalEnergies_spin;
    }
    else {
      nob = _nOcc_spin;
      nvb = _nVirt_spin;
      cb = _coefficients_spin;
      eb = _orbitalEnergies_spin;
    }
    iSpin = iSpin + 1;
  };

  if (noa > nob) {
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(ca.rows(), noa + nvb);
    Eigen::VectorXd e = Eigen::VectorXd::Zero(noa + nvb);

    C.middleCols(0, noa) = ca.middleCols(0, noa);
    e.middleRows(0, noa) = ea.middleRows(0, noa);

    C.middleCols(noa, nvb) = cb.middleCols(nob, nvb);
    e.middleRows(noa, nvb) = eb.middleRows(nob, nvb);

    iSpin = 0;
    for_spin(_nOcc, _nVirt, _coefficients, _orbitalEnergies) {
      if (iSpin == 0) {
        _nOcc_spin = noa;
        _nVirt_spin = nvb;

        _coefficients_spin = C;
        _orbitalEnergies_spin = e;
      }
      else {
        _nOcc_spin = 0;
        _nVirt_spin = 0;

        _coefficients_spin.setZero();
        _orbitalEnergies_spin.setZero();
      }
      iSpin = iSpin + 1;
    };
  }
  else {
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(ca.rows(), nob + nva);
    Eigen::VectorXd e = Eigen::VectorXd::Zero(nob + nva);

    C.middleCols(0, nob) = cb.middleCols(0, nob);
    e.middleRows(0, nob) = eb.middleRows(0, nob);

    C.middleCols(nob, nva) = ca.middleCols(noa, nva);
    e.middleRows(nob, nva) = ea.middleRows(noa, nva);

    iSpin = 0;
    for_spin(_nOcc, _nVirt, _coefficients, _orbitalEnergies) {
      if (iSpin == 0) {
        _nOcc_spin = nob;
        _nVirt_spin = nva;

        _coefficients_spin = C;
        _orbitalEnergies_spin = e;
      }
      else {
        _nOcc_spin = 0;
        _nVirt_spin = 0;

        _coefficients_spin.setZero();
        _orbitalEnergies_spin.setZero();
      }
      iSpin = iSpin + 1;
    };
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFController<SCFMode>::rotateOrbitalsSCFInstability() {
  printBigCaption("SCF Instability Orbital Rotation");
  if (_settings.scfstab != Options::STABILITY_ANALYSIS::REAL) {
    WarningTracker::printWarning(
        " You are following roots that are not result of a real stability analysis, make sure this is intended.", true);
  }
  OutputControl::n.printf(" Following instability root : %5i\n", _settings.stabroot);
  OutputControl::n.printf(" Orbital mixing parameter   : %5.1f\n\n", _settings.stabscal);
  if (_settings.stabroot > _settings.nEigen) {
    throw SerenityError(
        "The stabroot keyword must be lower or equal to the number of roots determined in the first place.");
  }

  unsigned iSpinStart = 0;
  for_spin(_coefficients, _nOcc, _nVirt) {
    Eigen::Map<Eigen::MatrixXd> xai((*_excitationVectors)[0].col(_settings.stabroot - 1).data() + iSpinStart,
                                    _nVirt_spin, _nOcc_spin);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(_nOcc_spin + _nVirt_spin, _nOcc_spin + _nVirt_spin);
    A.topRightCorner(_nOcc_spin, _nVirt_spin) = xai.transpose();
    A.bottomLeftCorner(_nVirt_spin, _nOcc_spin) = -xai;
    Eigen::MatrixXd U = matrixExp(_settings.stabscal * A);

    _coefficients_spin *= U;
    iSpinStart = _nOcc_spin * _nVirt_spin;
  };
  _system->getActiveOrbitalController<SCFMode>()->updateOrbitals(_coefficients, _orbitalEnergies);
}

template class LRSCFController<Options::SCF_MODES::RESTRICTED>;
template class LRSCFController<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
