/**
 * @file   BrokenSymmetry.cpp
 * @author Anja Massolle
 *
 * @date   17. November 2020
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
#include "analysis/brokenSymmetry/BrokenSymmetry.h"
/* Include Serenity Internal Headers */
#include "analysis/brokenSymmetry/CorrespondingOrbitals.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/SpinPolarizedData.h"
#include "energies/EnergyContributions.h"
#include "geometry/Geometry.h"
#include "io/Filesystem.h"
#include "misc/SerenityError.h"
#include "misc/WarningTracker.h"
#include "parameters/Constants.h"
#include "scf/SCFAnalysis.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/EvaluateEnergyTask.h"
#include "tasks/FDETask.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/LocalizationTask.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {
BrokenSymmetry::BrokenSymmetry(std::vector<std::shared_ptr<SystemController>> hsSystemController, EmbeddingSettings embedding,
                               bool tsOrtho, bool allOrtho, Options::ORTHOGONALIZATION_ALGORITHMS orthogonalizationScheme,
                               std::vector<std::shared_ptr<SystemController>> bsSystemController)
  : _hsSystemController(hsSystemController),
    _embedding(embedding),
    _tsOrtho(tsOrtho),
    _allOrtho(allOrtho),
    _orthogonalizationScheme(orthogonalizationScheme),
    _bsSystemController(bsSystemController) {
  if (_tsOrtho and _allOrtho) {
    throw SerenityError("evalTsOrtho and evalAllOrtho are exclusive to each other!");
  }
  if (_hsSystemController.size() == 1) {
    if (_hsSystemController[0]->getSpin() == 0) {
      throw SerenityError("Use the spin state of the high spin system!");
    }
  }
  else if (_hsSystemController.size() == 2) {
    if ((_hsSystemController[0]->getSpin() + _hsSystemController[1]->getSpin()) == 0) {
      throw SerenityError("Use the spin state of the high spin system for the subsystems!");
    }
  }

  if (_hsSystemController.size() == 2 and _bsSystemController.size() == 0) {
    auto subsys1Geom = _hsSystemController[0]->getGeometry();
    auto subsys2Geom = _hsSystemController[1]->getGeometry();

    auto brokenSymmetrySystemActSettings = _hsSystemController[0]->getSettings();
    brokenSymmetrySystemActSettings.spin = -1 * _hsSystemController[0]->getSpin();
    brokenSymmetrySystemActSettings.charge = _hsSystemController[0]->getCharge();
    brokenSymmetrySystemActSettings.name = "brokenSymmetryActSystem";
    _bsSystemController.push_back(std::make_shared<SystemController>(subsys1Geom, brokenSymmetrySystemActSettings));

    auto brokenSymmetrySystemEnvSettings = _hsSystemController[1]->getSettings();
    brokenSymmetrySystemEnvSettings.spin = _hsSystemController[1]->getSpin();
    brokenSymmetrySystemEnvSettings.charge = _hsSystemController[1]->getCharge();
    brokenSymmetrySystemEnvSettings.name = "brokenSymmetryEnvSystem";
    _bsSystemController.push_back(std::make_shared<SystemController>(subsys2Geom, brokenSymmetrySystemEnvSettings));
  }
  else if (_hsSystemController.size() > 2) {
    throw SerenityError("Broken-Symmetry calculations with more than two spin sites are not supported!");
  }
}

void BrokenSymmetry::bsDFT(unsigned int nA, unsigned int nB, Options::ORBITAL_LOCALIZATION_ALGORITHMS locType,
                           double noThreshold) {
  if (_hsSystemController.size() > 1) {
    throw SerenityError("Standard BS-DFT calculations are only implemented for one supersystem, for two subsystems use "
                        "an embedding calculation!");
  }
  if (_hsSystemController[0]->getSpin() == 0) {
    throw SerenityError("Use the spin state of the high spin system for the supersystem!");
  }

  printSubSectionTitle("Calculate High Spin System, Spin=" + std::to_string(_hsSystemController[0]->getSpin()));

  auto hsES = _hsSystemController[0]->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  auto hsCoeff = hsES->getMolecularOrbitals()->getCoefficients();
  auto nOcc = hsES->getNOccupiedOrbitals();
  auto s = hsES->getOneElectronIntegralController()->getOverlapIntegrals();
  auto hsDensMat = hsES->getDensityMatrix();
  auto basController = _hsSystemController[0]->getBasisController();
  auto nBas = basController->getNBasisFunctions();

  /*
   * Transform the density matrix into the MO basis.
   */
  Eigen::MatrixXd hsDensMatMOalpha = (hsCoeff.alpha).transpose() * s * hsDensMat.alpha * s * hsCoeff.alpha;
  Eigen::MatrixXd hsDensMatMObeta = (hsCoeff.alpha).transpose() * s * hsDensMat.beta * s * hsCoeff.alpha;
  /*
   * Calculate the natural orbitals (NO)
   */
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(hsDensMatMOalpha + hsDensMatMObeta);
  /*
   * Convert the NO from MO in AO basis
   */
  Eigen::MatrixXd noCoeff = hsCoeff.alpha * es.eigenvectors();
  auto noOccupation = es.eigenvalues();

  /*
   * Split the NOs into DONOs (doubly occupied), SONOs (singly occupied) and UONOs
   * (unoccupied), see Chem. Phys. Lett. 608, 50-54
   */
  CoefficientMatrix<RESTRICTED> donoCoeff(basController);
  donoCoeff.setZero();
  SpinPolarizedData<RESTRICTED, unsigned int> nDONO = 0;

  CoefficientMatrix<RESTRICTED> sonoCoeff(basController);
  sonoCoeff.setZero();
  SpinPolarizedData<RESTRICTED, unsigned int> nSONO = 0;

  CoefficientMatrix<RESTRICTED> uonoCoeff(basController);
  uonoCoeff.setZero();
  SpinPolarizedData<RESTRICTED, unsigned int> nUONO = 0;

  for (unsigned int i = 0; i < noOccupation.size(); i++) {
    if (abs(noOccupation[i] - 2) < noThreshold) {
      nDONO += 1;
      donoCoeff.col(nDONO - 1) = noCoeff.col(i);
    }
    else if (abs(noOccupation[i] - 1) < noThreshold) {
      nSONO += 1;
      sonoCoeff.col(nSONO - 1) = noCoeff.col(i);
    }
    else if (abs(noOccupation[i]) < noThreshold) {
      nUONO += 1;
      uonoCoeff.col(nUONO - 1) = noCoeff.col(i);
    }
  }
  /*
   * Sanity Checks
   */
  if (nSONO != (nA + nB)) {
    throw SerenityError("The sum of nA and nB does not match the number of SONOs!");
  }
  if (nBas != (nDONO + nSONO + nUONO)) {
    throw SerenityError("The sum of DONO, SONO and UONO does not match the number of orbitals!");
  }

  /*
   * Localize SONO orbitals (LNO)
   */
  auto sonoEigenValPtr = std::make_unique<SpinPolarizedData<RESTRICTED, Eigen::VectorXd>>(Eigen::VectorXd::Zero(nBas));
  auto sonoCoeffPtr = std::make_unique<CoefficientMatrix<RESTRICTED>>(sonoCoeff);

  auto sonoOrbs = std::make_shared<OrbitalController<RESTRICTED>>(std::move(sonoCoeffPtr), basController, *sonoEigenValPtr,
                                                                  _hsSystemController[0]->getNCoreElectrons());
  auto sonoES = std::make_shared<ElectronicStructure<RESTRICTED>>(
      _hsSystemController[0]->getOneElectronIntegralController(), nSONO, _hsSystemController[0]->getNCoreElectrons());

  sonoES->setMolecularOrbitals(sonoOrbs);
  auto sonoSettings = _hsSystemController[0]->getSettings();
  sonoSettings.scfMode = Options::SCF_MODES::RESTRICTED;
  sonoSettings.charge = _hsSystemController[0]->getNElectrons<RESTRICTED>() + _hsSystemController[0]->getCharge() - 2 * nSONO;
  sonoSettings.spin = 0;
  sonoSettings.path = _hsSystemController[0]->getSystemPath() + "tmp/";
  sonoSettings.name = "SONODummySystem";

  auto sonoSys = std::make_shared<SystemController>(_hsSystemController[0]->getGeometry(), sonoSettings);
  sonoSys->setElectronicStructure<RESTRICTED>(sonoES);

  auto sonoLoc = LocalizationTask(sonoSys);
  sonoLoc.settings.locType = locType;
  sonoLoc.run();
  auto lnoCoeff = sonoSys->getActiveOrbitalController<RESTRICTED>()->getCoefficients();

  /*
   * Construct the initial guess for the broken symmetry system
   */
  CoefficientMatrix<UNRESTRICTED> bsCoeff(basController);
  for_spin(bsCoeff) {
    (bsCoeff_spin).setZero();
    if (nDONO > 0) {
      (bsCoeff_spin).leftCols(nDONO) = donoCoeff.leftCols(nDONO);
    }
    if (nUONO > 0) {
      (bsCoeff_spin).rightCols(nUONO) = uonoCoeff.leftCols(nUONO);
    }
  };
  for (unsigned int i = 0; i < nA; i++) {
    (bsCoeff.alpha).col(nDONO + i) = lnoCoeff.col(i);
    (bsCoeff.beta).col(nDONO + nA + i) = lnoCoeff.col(i);
  }
  for (unsigned int i = 0; i < nB; i++) {
    (bsCoeff.beta).col(nDONO + i) = lnoCoeff.col(i + nA);
    (bsCoeff.alpha).col(nDONO + nA + i) = lnoCoeff.col(i + nA);
  }
  /*
   * Build BS system
   */
  auto bsSysSettings = _hsSystemController[0]->getSettings();
  bsSysSettings.spin = nA - nB;
  bsSysSettings.name = "BrokenSymmetrySystem";
  _bsSystemController.push_back(std::make_shared<SystemController>(_hsSystemController[0]->getGeometry(), bsSysSettings));
  SpinPolarizedData<UNRESTRICTED, unsigned int> nOccBS;
  nOccBS.alpha = nOcc.alpha - nB;
  nOccBS.beta = nOcc.beta + nB;

  auto bsOrbs =
      std::make_shared<OrbitalController<UNRESTRICTED>>(basController, _hsSystemController[0]->getNCoreElectrons() / 2);
  bsOrbs->updateOrbitals(bsCoeff, hsES->getMolecularOrbitals()->getEigenvalues());
  auto bsES = std::make_shared<ElectronicStructure<UNRESTRICTED>>(
      _hsSystemController[0]->getOneElectronIntegralController(), nOccBS, _hsSystemController[0]->getNCoreElectrons() / 2);
  bsES->setMolecularOrbitals(bsOrbs);
  _bsSystemController[0]->setElectronicStructure(bsES);

  printSubSectionTitle("Calculate Broken-Symmetry System, Spin=" + std::to_string(_bsSystemController[0]->getSpin()));

  ScfTask<UNRESTRICTED> brokenSymmetry(_bsSystemController[0]);
  brokenSymmetry.run();
  _bsSystemController[0]->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->toHDF5(
      _bsSystemController[0]->getHDF5BaseName(), _bsSystemController[0]->getSettings().identifier);

  removeSystemFiles(_hsSystemController[0]->getSystemPath() + "tmp/SONODummySystem/", "SONODummySystem");
  std::remove((_hsSystemController[0]->getSystemPath() + "tmp/").c_str());
} /* bsDFT() */

void BrokenSymmetry::bsIsolated() {
  printSubSectionTitle("Calculate High Spin System, Spin=" +
                       std::to_string(_hsSystemController[0]->getSpin() + _hsSystemController[1]->getSpin()));
  for (unsigned int i = 0; i < _hsSystemController.size(); i++) {
    _hsSystemController[i]->getElectronicStructure<UNRESTRICTED>();
  }
  EvaluateEnergyTask<UNRESTRICTED> evalHS(_hsSystemController);
  evalHS.settings.XCfunctional = (_hsSystemController[0]->getSettings()).dft.functional;
  evalHS.settings.embedding = _embedding;
  evalHS.run();
  printSubSectionTitle("Calculate Broken-Symmetry System, Spin=" +
                       std::to_string(_bsSystemController[0]->getSpin() + _bsSystemController[1]->getSpin()));
  for (unsigned int i = 0; i < _bsSystemController.size(); i++) {
    _bsSystemController[i]->getElectronicStructure<UNRESTRICTED>();
  }
  EvaluateEnergyTask<UNRESTRICTED> evalBS(_bsSystemController);
  evalBS.settings.XCfunctional = (_bsSystemController[0]->getSettings()).dft.functional;
  evalBS.settings.embedding = _embedding;
  evalBS.run();
} /*bsIsolated()*/

void BrokenSymmetry::bsFDE() {
  // save for parallel FDE
  auto highSpinActCopySettings = _hsSystemController[0]->getSettings();
  highSpinActCopySettings.name = "highSpinActIsolated";
  highSpinActCopySettings.path = _hsSystemController[0]->getSystemPath() + "tmp/";
  auto highSpinActCopy = std::make_shared<SystemController>(_hsSystemController[0]->getGeometry(), highSpinActCopySettings);

  FDETask<UNRESTRICTED> highSpin({_hsSystemController[0]}, {_hsSystemController[1]});
  highSpin.settings.embedding = _embedding;
  highSpin.settings.calculateEnvironmentEnergy = false;
  printSubSectionTitle("Calculate High Spin System, Spin=" +
                       std::to_string(_hsSystemController[0]->getSpin() + _hsSystemController[1]->getSpin()));
  highSpin.run();

  FDETask<UNRESTRICTED> highSpinParallel({_hsSystemController[1]}, {highSpinActCopy});
  highSpinParallel.settings = highSpin.settings;
  highSpinParallel.run();

  // save for parallel FDE
  auto brokenSymmetrySystemActCopySettings = _bsSystemController[0]->getSettings();
  brokenSymmetrySystemActCopySettings.name = "brokenSymmetryActIsolated";
  brokenSymmetrySystemActCopySettings.path = _bsSystemController[0]->getSystemPath() + "tmp/";
  auto brokenSymmetryActSystemCopy =
      std::make_shared<SystemController>(_bsSystemController[0]->getGeometry(), brokenSymmetrySystemActCopySettings);

  FDETask<UNRESTRICTED> brokenSymmetry({_bsSystemController[0]}, {_bsSystemController[1]});
  brokenSymmetry.settings = highSpin.settings;
  printSubSectionTitle("Calculate Broken-Symmetry System, Spin=" +
                       std::to_string(_bsSystemController[0]->getSpin() + _bsSystemController[1]->getSpin()));
  brokenSymmetry.run();

  FDETask<UNRESTRICTED> brokenSymmetryParallel({_bsSystemController[1]}, {brokenSymmetryActSystemCopy});
  brokenSymmetryParallel.settings = brokenSymmetry.settings;
  brokenSymmetryParallel.run();
  removeSystemFiles(_hsSystemController[0]->getSystemPath() + "tmp/highSpinActIsolated/", "highSpinActIsolated");
  removeSystemFiles(_bsSystemController[0]->getSystemPath() + "tmp/brokenSymmetryActIsolated/",
                    "brokenSymmetryActIsolated");
  std::remove((_hsSystemController[0]->getSystemPath() + "tmp/").c_str());
  std::remove((_bsSystemController[0]->getSystemPath() + "tmp/").c_str());
} /*bsFDE()*/

void BrokenSymmetry::bsFAT(FreezeAndThawTaskSettings fatSettings) {
  FreezeAndThawTask<UNRESTRICTED> highSpin(_hsSystemController, {});
  highSpin.settings.embedding = _embedding;
  highSpin.settings.maxCycles = fatSettings.maxCycles;
  highSpin.settings.convThresh = fatSettings.convThresh;
  printSubSectionTitle("Calculate High Spin System, Spin=" +
                       std::to_string(_hsSystemController[0]->getSpin() + _hsSystemController[1]->getSpin()));
  highSpin.run();

  FreezeAndThawTask<UNRESTRICTED> brokenSymmetry(_bsSystemController, {});
  brokenSymmetry.settings = highSpin.settings;
  printSubSectionTitle("Calculate Broken-Symmetry System, Spin=" +
                       std::to_string(_bsSystemController[0]->getSpin() + _bsSystemController[1]->getSpin()));
  brokenSymmetry.run();
} /*bsFAT()*/

void BrokenSymmetry::tsOrtho() {
  printSubSectionTitle("Calculate the non-additive kinetic energy from orthogonalized supersystem orbitals");
  EvaluateEnergyTask<UNRESTRICTED> evalHS(_hsSystemController);
  evalHS.settings.orthogonalizationScheme = _orthogonalizationScheme;
  evalHS.settings.evalTsOrtho = _tsOrtho;
  evalHS.settings.embedding = _embedding;
  evalHS.run();
  EvaluateEnergyTask<UNRESTRICTED> evalBS(_bsSystemController);
  evalBS.settings.orthogonalizationScheme = _orthogonalizationScheme;
  evalBS.settings.evalTsOrtho = _tsOrtho;
  evalBS.settings.embedding = _embedding;
  evalBS.run();
} /*tsOrtho*/

void BrokenSymmetry::allOrtho() {
  printSubSectionTitle("Calculate all energy contributions from orthogonalized supersystem orbitals");
  Settings superSettings = _hsSystemController[0]->getSettings();
  superSettings.name = "superHS";
  superSettings.charge = 0;
  superSettings.spin = 0;
  auto superHS = std::make_shared<SystemController>(std::make_shared<Geometry>(), superSettings);

  EvaluateEnergyTask<UNRESTRICTED> evalHS(_hsSystemController, superHS);
  evalHS.settings.evalAllOrtho = true;
  evalHS.run();

  _hsSystemController = {superHS};

  superSettings = _bsSystemController[0]->getSettings();
  superSettings.name = "superBS";
  superSettings.charge = 0;
  superSettings.spin = 0;

  auto superBS = std::make_shared<SystemController>(std::make_shared<Geometry>(), superSettings);

  EvaluateEnergyTask<UNRESTRICTED> evalBS(_bsSystemController, superBS);
  evalBS.settings.evalAllOrtho = true;
  evalBS.run();

  _bsSystemController = {superBS};
} /*allOrtho*/

void BrokenSymmetry::evalJ() {
  /*
   * Calculate S*S
   */
  if (!_bsSystemController[0]->hasElectronicStructure<UNRESTRICTED>()) {
    throw SerenityError("No Broken-Symmetry state was calculated, please perform a BS-DFT calculation!");
  }
  if (_tsOrtho) {
    this->tsOrtho();
  }

  else if (_allOrtho) {
    this->allOrtho();
  }

  auto hsAnalysis = std::make_shared<SCFAnalysis<UNRESTRICTED>>(_hsSystemController);
  auto bsAnalysis = std::make_shared<SCFAnalysis<UNRESTRICTED>>(_bsSystemController);

  _s2HS = hsAnalysis->getS2();
  _s2BS = bsAnalysis->getS2();

  if (_hsSystemController.size() == 1 and _bsSystemController.size() == 1) {
    _s2UHFhs = hsAnalysis->getS2(true);
    _s2UHFbs = bsAnalysis->getS2(true);
  }

  /*
   * Calculate J coupling
   */
  auto eHs = _hsSystemController[0]->getElectronicStructure<UNRESTRICTED>()->getEnergy();
  auto eBs = _bsSystemController[0]->getElectronicStructure<UNRESTRICTED>()->getEnergy();
  double sMax = 0;

  if (_hsSystemController.size() == 1 and _bsSystemController.size() == 1) {
    sMax = 0.5 * _hsSystemController[0]->getSpin();
  }
  else if (_hsSystemController.size() == 2 and _bsSystemController.size() == 2) {
    sMax = 0.5 * (_hsSystemController[0]->getSpin() + _hsSystemController[1]->getSpin());
  }

  _j1 = (eBs - eHs) * HARTREE_TO_OOCM / (sMax * sMax);
  _j2 = (eBs - eHs) * HARTREE_TO_OOCM / (sMax * (sMax + 1));
  _j3 = (eBs - eHs) * HARTREE_TO_OOCM / (_s2HS - _s2BS);

  if (_hsSystemController.size() == 1 and _bsSystemController.size() == 1) {
    _j3UHF = (eBs - eHs) * HARTREE_TO_OOCM / (_s2UHFhs - _s2UHFbs);

    auto coeff = _bsSystemController[0]->getActiveOrbitalController<UNRESTRICTED>()->getCoefficients();
    auto nOcc = _bsSystemController[0]->getNOccupiedOrbitals<UNRESTRICTED>();
    CoefficientMatrix<RESTRICTED> bsAlpha(coeff.getBasisController());
    CoefficientMatrix<RESTRICTED> bsBeta(coeff.getBasisController());
    bsAlpha = coeff.alpha;
    bsBeta = coeff.beta;
    CorrespondingOrbitals<RESTRICTED> correspondingOrbs(bsAlpha, bsBeta, nOcc.alpha, nOcc.beta);
    _sAB = correspondingOrbs.getOverlap(nOcc.alpha);
    _j4 = (eBs - eHs) * HARTREE_TO_OOCM / (1 + pow(_sAB, 2));
  }
} /*evalJ*/

double BrokenSymmetry::getS2HS() {
  if (_s2HS == std::numeric_limits<double>::infinity()) {
    this->evalJ();
  }
  return _s2HS;
};

double BrokenSymmetry::getS2BS() {
  if (_s2BS == std::numeric_limits<double>::infinity()) {
    this->evalJ();
  }
  return _s2BS;
};

double BrokenSymmetry::getSab() {
  if (_sAB == std::numeric_limits<double>::infinity()) {
    this->evalJ();
  }
  return _sAB;
};

double BrokenSymmetry::getS2uhfHS() {
  if (_hsSystemController.size() == 1) {
    if (_s2UHFhs == std::numeric_limits<double>::infinity()) {
      this->evalJ();
    }
    return _s2UHFhs;
  }
  else {
    throw SerenityError("Orbital based S*S evaluation is only implemented for orthogonal supersystem orbitals.");
  }
};

double BrokenSymmetry::getS2uhfBS() {
  if (_hsSystemController.size() == 1) {
    if (_s2UHFbs == std::numeric_limits<double>::infinity()) {
      this->evalJ();
    }
    return _s2UHFbs;
  }
  else {
    throw SerenityError("Orbital based S*S evaluation is only implemented for orthogonal supersystem orbitals.");
  }
};

double BrokenSymmetry::getJ1() {
  if (_j1 == std::numeric_limits<double>::infinity()) {
    this->evalJ();
  }
  return _j1;
};

double BrokenSymmetry::getJ2() {
  if (_j2 == std::numeric_limits<double>::infinity()) {
    this->evalJ();
  }
  return _j2;
};

double BrokenSymmetry::getJ3() {
  if (_j3 == std::numeric_limits<double>::infinity()) {
    this->evalJ();
  }
  return _j3;
};

double BrokenSymmetry::getJ3UHF() {
  if (_hsSystemController.size() == 1) {
    if (_j3UHF == std::numeric_limits<double>::infinity()) {
      this->evalJ();
    }
    return _j3UHF;
  }
  else {
    throw SerenityError(
        "J3 with orbital based S*S evaluation is only implemented for orthogonal supersystem orbitals.");
  }
};

double BrokenSymmetry::getJ4() {
  if (_hsSystemController.size() == 1) {
    if (_j4 == std::numeric_limits<double>::infinity()) {
      this->evalJ();
    }
    return _j4;
  }
  else {
    throw SerenityError("J4 is only implemented for orthogonal supersystem orbitals.");
  }
};

} /*namespace Serenity*/