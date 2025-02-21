/**
 * @file OrbitalsIOTask.cpp
 *
 * @date Dec 9, 2020
 * @author Moritz Bensberg
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
#include "tasks/OrbitalsIOTask.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"     //Basis index map for molcas.
#include "basis/Basis.h"                           //Shell-wise sorting.
#include "basis/BasisController.h"                 //Basis.
#include "basis/CartesianToSphericalTransformer.h" //Transform from spherical to cartesian basis.
#include "data/ElectronicStructure.h"              //Constructor.
#include "data/OrbitalController.h"                //Constructor.
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h" //Overlap integrals.
#include "io/Eigen3HDF5.h"
#include "io/HDF5.h"                            //Load from HDF5.
#include "math/linearAlgebra/MatrixFunctions.h" //Cholesky orthogonalization.
#include "misc/SerenityError.h"                 //Errors.
#include "misc/WarningTracker.h"                //Warnings.
#include "system/SystemController.h"            //Access to system properties.
/* Include Std and External Headers */
#include <algorithm> //Replace within a string.
#include <fstream>   //input/output streams.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
OrbitalsIOTask<SCFMode>::OrbitalsIOTask(std::shared_ptr<SystemController> system) : _system(system) {
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::run() {
  bool updateElecStructure = false;
  auto basisController = _system->getBasisController();
  auto finalCoeffPtr = std::make_unique<CoefficientMatrix<SCFMode>>(basisController);
  auto eigenvaluesPtr = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(
      Eigen::VectorXd::Zero(basisController->getNBasisFunctions()));
  if (settings.replaceInFile && settings.fileFormat != Options::ORBITAL_FILE_TYPES::MOLCAS)
    throw SerenityError("Replacing orbitals in external files is only supported for Molcas HDF5 files."
                        "Please change the settings of the task accordingly.");
  if (settings.write) {
    if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::TURBOMOLE) {
      initTurboSortingMatrices();
      checkTurboInput();
      auto coefficients = _system->getActiveOrbitalController<SCFMode>()->getCoefficients();
      auto eigenvalues = _system->getActiveOrbitalController<SCFMode>()->getEigenvalues();
      sortMatrixAndVector(coefficients, eigenvalues);
      writeTurbomoleOrbitals(coefficients, eigenvalues);
    }
    else if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::MOLDEN) {
      if (basisController->getMaxAngularMomentum() > 4)
        throw SerenityError((std::string) "ERROR: The maximum angular momentum which Molden supports was\n" +
                            "       exceeded. Supported l= 4");
      writeMoldenOrbitals();
    }
    else {
      throw SerenityError("Writing of files currently only supported in the Turbomole and the Molden format");
    }
  }
  else {
    if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::TURBOMOLE) {
      updateElecStructure = true;
      initTurboSortingMatrices();
      checkTurboInput();
      readTurbomoleOrbitals(*finalCoeffPtr, *eigenvaluesPtr);
    }
    else if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::MOLPRO) {
      initMolproSortingMatrices();
      checkMolproInput();
      updateElecStructure = true;
      this->readMolproXMLOrbitals(*finalCoeffPtr, *eigenvaluesPtr);
    }
    else if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::MOLCAS) {
      if (SCFMode != RESTRICTED)
        throw SerenityError("Reading Molcas orbitals is currently only supported for spin-restricted calculations.");
      if (settings.replaceInFile) {
        this->replaceMolcasHDF5Orbitals(
            _system->template getElectronicStructure<SCFMode>()->getMolecularOrbitals()->getCoefficients());
      }
      else {
        this->readMolcasHDF5Orbitals(*finalCoeffPtr, *eigenvaluesPtr);
        updateElecStructure = true;
      }
    }
    else if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::MOLDEN) {
      throw SerenityError("Reading of Molden files currently not supported.");
    }
    if (updateElecStructure) {
      auto newOrbitalController = std::make_shared<OrbitalController<SCFMode>>(
          std::move(finalCoeffPtr), basisController, *eigenvaluesPtr, _system->getNCoreElectrons());
      auto targetElectronicStructure =
          std::make_shared<ElectronicStructure<SCFMode>>(newOrbitalController, _system->getOneElectronIntegralController(),
                                                         _system->template getNOccupiedOrbitals<SCFMode>());
      _system->template setElectronicStructure<SCFMode>(targetElectronicStructure);
      targetElectronicStructure->toHDF5(_system->getHDF5BaseName(), _system->getSystemIdentifier());
    }
    if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::SERENITY) {
      std::string systemId = getSerenityIDFromFile();
      auto targetElectronicStructure = std::make_shared<ElectronicStructure<SCFMode>>(
          settings.path + "/" + _system->getSystemName(), _system->getBasisController(), _system->getGeometry(),
          _system->getExternalChargeController(), systemId);
      _system->template setElectronicStructure<SCFMode>(targetElectronicStructure);
      if (settings.resetCoreOrbitals)
        _system->template getActiveOrbitalController<SCFMode>()->setCoreOrbitalsByNumber(_system->getNCoreElectrons() / 2);
      targetElectronicStructure->toHDF5(_system->getHDF5BaseName(), _system->getSystemIdentifier());
    }
  }
}

template<Options::SCF_MODES SCFMode>
std::string OrbitalsIOTask<SCFMode>::getSerenityIDFromFile() {
  std::string filePath = settings.path + "/" + _system->getSystemName();
  if (SCFMode == RESTRICTED) {
    filePath += ".orbs.res.h5";
  }
  else {
    filePath += ".orbs.unres.h5";
  }
  HDF5::Filepath name(filePath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::attribute_exists(file, "ID");
  std::string id = HDF5::load_attribute<std::string>(file, "ID");
  file.close();
  return id;
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::checkTurboInput() {
  auto basisController = _system->getBasisController();
  if (basisController->getMaxAngularMomentum() > _maxL)
    throw SerenityError((std::string) "ERROR: The maximum angular momentum for which Turbomole sorting matrices are\n" +
                        "       available was exceeded. Currently supported l= " + std::to_string(_maxL));
  if (not basisController->isPureSpherical()) {
    throw SerenityError((std::string) "ERROR: Only spherical harmonics are supported when reading Turbomole orbitals.");
  }
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::checkMolproInput() {
  auto basisController = _system->getBasisController();
  if (SCFMode != RESTRICTED)
    throw SerenityError((std::string) "ERROR: Reading Molpro-xml files is only supported for RESTRICTED orbitals.");
  if (basisController->getMaxAngularMomentum() > _maxL)
    throw SerenityError((std::string) "ERROR: The maximum angular momentum for which Molpro sorting matrices are\n" +
                        "       available was exceeded. Currently supported l= " + std::to_string(_maxL));
  if (not basisController->isPureSpherical()) {
    throw SerenityError((std::string) "ERROR: Only spherical harmonics are supported when reading Molpro orbitals.");
  }
}

template<>
void OrbitalsIOTask<RESTRICTED>::readTurbomoleOrbitals(CoefficientMatrix<RESTRICTED>& coeffs,
                                                       SpinPolarizedData<RESTRICTED, Eigen::VectorXd>& eigenvalues) {
  unsigned int nBFs = coeffs.getBasisController()->getNBasisFunctions();
  auto coeff_and_eps = readTurbomoleOrbitalFile(settings.path + "/mos", nBFs);
  coeffs.block(0, 0, nBFs, nBFs) = coeff_and_eps.first;
  eigenvalues = coeff_and_eps.second;
}

template<>
void OrbitalsIOTask<UNRESTRICTED>::readTurbomoleOrbitals(CoefficientMatrix<UNRESTRICTED>& coeffs,
                                                         SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>& eigenvalues) {
  unsigned int nBFs = coeffs.getBasisController()->getNBasisFunctions();
  auto coeff_and_eps_alpha = readTurbomoleOrbitalFile(settings.path + "/alpha", nBFs);
  coeffs.alpha.block(0, 0, nBFs, nBFs) = coeff_and_eps_alpha.first;
  eigenvalues.alpha = coeff_and_eps_alpha.second;
  auto coeff_and_eps_beta = readTurbomoleOrbitalFile(settings.path + "/beta", nBFs);
  coeffs.beta.block(0, 0, nBFs, nBFs) = coeff_and_eps_beta.first;
  eigenvalues.beta = coeff_and_eps_beta.second;
}

template<>
void OrbitalsIOTask<RESTRICTED>::writeTurbomoleOrbitals(CoefficientMatrix<RESTRICTED>& coefficients,
                                                        SpinPolarizedData<RESTRICTED, Eigen::VectorXd>& eigenvalues) {
  unsigned int nBFs = coefficients.getBasisController()->getNBasisFunctions();
  std::string filePath = settings.path + "/" + _system->getSystemName();
  std::ofstream file(filePath + "/mos");
  file << "$scfmo scfconv=7 format(4d20.14)" << std::endl
       << std::endl
       << std::endl; // scfconv=7 is the Turbomole standard for restricted
  file.close();
  writeTurbomoleOrbitalFile(filePath + "/mos", nBFs, eigenvalues, resortCoefficients(coefficients));
}

template<>
void OrbitalsIOTask<UNRESTRICTED>::writeTurbomoleOrbitals(CoefficientMatrix<UNRESTRICTED>& coefficients,
                                                          SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>& eigenvalues) {
  unsigned int nBFs = coefficients.getBasisController()->getNBasisFunctions();
  std::string filePath = settings.path + "/" + _system->getSystemName();
  std::ofstream filea(filePath + "/alpha");
  filea << "$uhfmo_alpha scfconv=8 format(4d20.14)" << std::endl
        << std::endl
        << std::endl; // scfconv=8 is the Turbomole standard for unrestricted
  filea.close();
  writeTurbomoleOrbitalFile(filePath + "/alpha", nBFs, eigenvalues.alpha, resortCoefficients(coefficients.alpha));
  std::ofstream fileb(filePath + "/beta");
  fileb << "$uhfmo_beta scfconv=8 format(4d20.14)" << std::endl
        << std::endl
        << std::endl; // scfconv=8 is the Turbomole standard for unrestricted
  fileb.close();
  writeTurbomoleOrbitalFile(filePath + "/beta", nBFs, eigenvalues.beta, resortCoefficients(coefficients.beta));
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::writeMoldenOrbitals() {
  std::ofstream file(settings.path + "/" + _system->getSystemName() + ".molden.input");
  file << "[Molden Format]" << std::endl << "[Title]" << std::endl << std::endl;
  file << "[Atoms] AU" << std::endl;
  const auto& geometry = *_system->getGeometry();
  for (unsigned int i = 0; i < geometry.getNAtoms(); i++) {
    const auto& atom = geometry.getAtoms()[i];
    const auto coords = atom->coords();
    file << geometry.getAtomSymbols()[i] << "   " << i + 1 << "   " << atom->getNuclearCharge() << " "
         << reformatNumber(coords[0], "E") << " " << reformatNumber(coords[1], "E") << " "
         << reformatNumber(coords[2], "E") << std::endl;
  }
  file << "[GTO]" << std::endl;
  bool spherical = false;
  for (unsigned int i = 0; i < geometry.getNAtoms(); i++) {
    const auto& atom = geometry.getAtoms()[i];
    const auto shells = atom->getBasisFunctions();
    file << "   " << i + 1 << " 0" << std::endl;
    for (auto shell : shells) {
      spherical = shell->isSpherical();
      unsigned int nPrim = shell->getNPrimitives();
      unsigned int l = shell->getAngularMomentum();
      switch (l) {
        case 0:
          file << "s   " << nPrim << " 1.0" << std::endl;
          break;
        case 1:
          file << "p   " << nPrim << " 1.0" << std::endl;
          break;
        case 2:
          file << "d   " << nPrim << " 1.0" << std::endl;
          break;
        case 3:
          file << "f   " << nPrim << " 1.0" << std::endl;
          break;
        case 4:
          file << "g   " << nPrim << " 1.0" << std::endl;
          break;
      }
      const auto& exp = shell->getExponents();
      const auto& contractions = shell->getNormContractions();
      Eigen::VectorXd contr = Eigen::VectorXd::Zero(nPrim);
      for (unsigned int p = 0; p < nPrim; p++) {
        contr(p) = contractions[p];
        /*
         * For a cartesian basis the normalization must be reverted. This is also done in other codes such as STDA
         * (readbasmold.f; ll. 388 - 403). In this case a cartesian Molden file is written. For a spherical
         * basis a spherical Molden file is written.
         */
        if (!spherical)
          contr(p) /= (pow(exp[p] / (2 * pow(PI, 3)), 0.25) * pow(2 * sqrt(exp[p]), (l + 1)));
        file << reformatNumber(exp[p], "E") << " ";
        if (shell->getAngularMomentum() == 2) {
          file << reformatNumber(contr(p) * sqrt(3), "E") << std::endl;
        }
        else if (shell->getAngularMomentum() == 3) {
          file << reformatNumber(contr(p) * sqrt(15), "E") << std::endl;
        }
        else if (shell->getAngularMomentum() == 4) {
          if (spherical) {
            file << reformatNumber(contr(p) * sqrt(35), "E") << std::endl;
          }
          else {
            file << reformatNumber(contr(p) * sqrt(105), "E") << std::endl;
          }
        }
        else {
          file << reformatNumber(contr(p), "E") << std::endl;
        }
      }
    }
    file << std::endl;
  }
  if (spherical) {
    file << "[5D]" << std::endl << "[7F]" << std::endl << "[9G]" << std::endl;
    initSphericalMoldenSortingMatrices();
  }
  else {
    initCartesianMoldenSortingMatrices();
  }
  file << "[MO]" << std::endl;
  auto coefficients = _system->getActiveOrbitalController<SCFMode>()->getCoefficients();
  auto eigenvalues = _system->getActiveOrbitalController<SCFMode>()->getEigenvalues();
  sortMatrixAndVector(coefficients, eigenvalues);
  const auto nocc = _system->getElectronicStructure<SCFMode>()->getNOccupiedOrbitals();
  bool isAlpha = true;
  for_spin(coefficients, eigenvalues, nocc) {
    auto C = resortCoefficients(coefficients_spin);
    for (unsigned int i = 0; i < C.cols(); i++) {
      if (eigenvalues_spin(i) == std::numeric_limits<double>::infinity()) {
        continue;
      }
      file << "Sym=  " << i + 1 << "a" << std::endl;
      file << "Ene= " << reformatNumber(eigenvalues_spin(i), "E") << std::endl;
      if (isAlpha) {
        file << "Spin= Alpha" << std::endl;
      }
      else {
        file << "Spin= Beta" << std::endl;
      }
      if (i < nocc_spin && SCFMode == Options::SCF_MODES::RESTRICTED) {
        file << "Occup= 2.000000" << std::endl;
      }
      else if (i < nocc_spin && SCFMode == Options::SCF_MODES::UNRESTRICTED) {
        file << "Occup= 1.000000" << std::endl;
      }
      else {
        file << "Occup= 0.000000" << std::endl;
      }
      for (unsigned int j = 0; j < C.cols(); j++) {
        file << j + 1 << "  " << reformatNumber(C(j, i), "E") << std::endl;
      }
    }
    isAlpha = false;
  };
  file.close();
}

template<Options::SCF_MODES SCFMode>
double OrbitalsIOTask<SCFMode>::getTurboEigenvalue(std::istringstream& iss) {
  std::string word;
  iss >> word;
  std::replace(word.begin(), word.end(), 'D', 'e');
  return std::stod(word.substr(11));
}
template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::checkTurboNAOs(std::istringstream& iss, unsigned int nBFs) {
  std::string word;
  iss >> word;
  unsigned int nAOs = std::stoi(word.substr(6));
  if (nAOs != nBFs) {
    throw SerenityError((std::string) "ERROR: Inconsistent number of basis functions in file!");
  }
}
template<Options::SCF_MODES SCFMode>
Eigen::VectorXd OrbitalsIOTask<SCFMode>::getTurboOrbital(std::ifstream& input, unsigned int nBFs, unsigned int nPerLine,
                                                         unsigned int length) {
  unsigned int nLines = nBFs / nPerLine;
  unsigned int nRemain = nBFs % nPerLine;
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(nBFs);
  unsigned int counter = 0;
  std::string line;
  for (unsigned int iLine = 0; iLine < nLines; ++iLine) {
    if (getline(input, line)) {
      splitTurboOrbLine(line, nPerLine, length, coeffs, counter);
    }
    else {
      throw SerenityError((std::string) "ERROR: Unexpected file format in file!");
    }
  }
  if (nRemain > 0) {
    if (getline(input, line)) {
      splitTurboOrbLine(line, nRemain, length, coeffs, counter);
    }
    else {
      throw SerenityError((std::string) "ERROR: Unexpected file format in file!");
    }
  }
  return coeffs;
}
template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::splitTurboOrbLine(std::string& line, unsigned int nPerLine, unsigned int length,
                                                Eigen::VectorXd& coeffs, unsigned int& counter) {
  unsigned int start = 0;
  for (unsigned int iCoeff = 0; iCoeff < nPerLine; ++iCoeff) {
    std::string coeff_str = line.substr(start, length);
    std::replace(coeff_str.begin(), coeff_str.end(), 'D', 'e');
    coeffs(counter) = std::stod(coeff_str);
    ++counter;
    start += length;
  } // for iCoeff
}

template<Options::SCF_MODES SCFMode>
std::string OrbitalsIOTask<SCFMode>::reformatNumber(double number, std::string exponentSymbol) {
  number = round(number * 1e13) / 1e13;
  int exponent = floor(log10(std::abs(number)));
  double base = number / pow(10, exponent);
  if (abs(number) < 10e-14) {
    exponent = 0;
    base = 0.0;
  }
  base /= 10.0;
  exponent += 1;
  std::stringstream exp;
  if (exponent < -9)
    exp << exponentSymbol << exponent;
  else if (exponent > -10 && exponent < 0)
    exp << exponentSymbol << "-0" << abs(exponent);
  else if (exponent > -1 && exponent < 10)
    exp << exponentSymbol << "+0" << exponent;
  else if (exponent > 9)
    exp << exponentSymbol << "+" << exponent;
  std::stringstream outputss;
  outputss << std::fixed << std::setprecision(14) << base << exp.str();
  std::string output = outputss.str();
  if (base < 0) {
    output.erase(0, 1);
    output[0] = '-';
  }
  return output;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd OrbitalsIOTask<SCFMode>::resortCoefficients(const Eigen::MatrixXd& coefficients) {
  auto basis = _system->getBasisController()->getBasis();
  unsigned int startRowOrig = 0;
  unsigned int startRowTarget = 0;
  unsigned int nOrbs = _system->getBasisController()->getNBasisFunctions();
  Eigen::MatrixXd finalCoefficients = Eigen::MatrixXd::Zero(nOrbs, nOrbs);
  for (const auto& shell : basis) {
    unsigned int angularMomentum = shell->getAngularMomentum();
    unsigned int nContractedOrig = shell->getNContracted();
    unsigned int nContractedTarget = shell->getNContracted();
    if (settings.write) {
      auto sortMatrix = *_sortMatrices[angularMomentum];
      finalCoefficients.block(startRowTarget, 0, nContractedTarget, nOrbs) =
          (sortMatrix.transpose() * coefficients.block(startRowOrig, 0, nContractedOrig, nOrbs)).eval();
    }
    else {
      finalCoefficients.block(startRowTarget, 0, nContractedTarget, nOrbs) =
          (*_sortMatrices[angularMomentum] * coefficients.block(startRowOrig, 0, nContractedOrig, nOrbs)).eval();
    }
    startRowOrig += nContractedOrig;
    startRowTarget += nContractedTarget;
  } // for shell
  return finalCoefficients;
}

template<Options::SCF_MODES SCFMode>
std::pair<Eigen::MatrixXd, Eigen::VectorXd> OrbitalsIOTask<SCFMode>::readTurbomoleOrbitalFile(std::string filePath,
                                                                                              unsigned int nBFs) {
  std::string line;
  std::ifstream input(filePath);
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::Zero(nBFs);
  Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(nBFs, nBFs);
  // Read format
  if (getline(input, line)) {
    std::string word, dummy;
    std::istringstream iss(line);
    iss >> dummy;
    iss >> dummy;
    iss >> word;
    if (word != "format(4d20.14)") {
      throw SerenityError((std::string) "ERROR: Unexpected file format in file " + filePath);
    }
  }
  else {
    throw SerenityError((std::string) "ERROR: Unexpected file format in file " + filePath);
  }
  unsigned int length = 20; // 2+14+4
  unsigned int nPerLine = 4;
  // skip two lines
  if (not getline(input, line))
    throw SerenityError((std::string) "ERROR: Unexpected file format in file " + filePath);
  if (not getline(input, line))
    throw SerenityError((std::string) "ERROR: Unexpected file format in file " + filePath);
  unsigned int iOrb = 0;
  while (getline(input, line)) {
    std::string word, dummy;
    std::istringstream iss(line);
    // Skip the two comment lines
    iss >> dummy;
    if (dummy == "$end")
      break;
    iss >> dummy;
    // get eigenvalue and check the number of AOs
    eigenvalues(iOrb) = getTurboEigenvalue(iss);
    checkTurboNAOs(iss, nBFs);
    // get the coefficients
    coefficients.col(iOrb) = getTurboOrbital(input, nBFs, nPerLine, length);
    ++iOrb;
  }
  if (iOrb < nBFs) {
    WarningTracker::printWarning((std::string) "Warning: The number of orbitals read from the mos-file does\n" +
                                     "         not correspond to the total number of orbitals in\n" + "         the system.",
                                 true);
  }
  coefficients = resortCoefficients(coefficients);
  checkNorm(coefficients);
  coefficients = orthogonalize_chol(coefficients, _system->getOneElectronIntegralController()->getOverlapIntegrals());
  return std::make_pair(coefficients, eigenvalues);
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::writeTurbomoleOrbitalFile(std::string filePath, unsigned int nBFs,
                                                        const Eigen::VectorXd& eigenvalues,
                                                        const Eigen::MatrixXd& coefficients) {
  std::ofstream file;
  file.open(filePath, std::ios_base::app);
  unsigned int nPerLine = 4;
  unsigned int nLines = nBFs / nPerLine;
  unsigned int nRemain = nBFs % nPerLine;
  for (unsigned int i = 0; i < nBFs; i++) {
    file << " " << i + 1 << " a eigenvalue=" << reformatNumber(eigenvalues(i)) << " nsaos=" << nBFs;
    Eigen::VectorXd coeffCol = coefficients.col(i);
    for (unsigned int j = 0; j < nLines; j++) {
      file << std::endl << reformatNumber(coeffCol(j * nPerLine + 0)) << reformatNumber(coeffCol(j * nPerLine + 1));
      file << reformatNumber(coeffCol(j * nPerLine + 2)) << reformatNumber(coeffCol(j * nPerLine + 3));
    }
    if (nRemain > 0) {
      file << std::endl;
      for (unsigned int k = 0; k < nRemain; k++) {
        file << reformatNumber(coeffCol(nLines * nPerLine + k));
      }
    }
    file << std::endl;
  }
  file << "$end" << std::endl;
  file.close();
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::checkNorm(const Eigen::MatrixXd& coefficients) {
  unsigned int nOrbs = coefficients.cols();
  const MatrixInBasis<RESTRICTED> s = _system->getOneElectronIntegralController()->getOverlapIntegrals();
  const Eigen::MatrixXd s_mo = (coefficients.transpose() * s * coefficients).eval();
  for (unsigned int iOrb = 0; iOrb < nOrbs; ++iOrb) {
    double check = s_mo.col(iOrb).array().abs().sum() - s_mo(iOrb, iOrb);
    if (check > 1e-6) {
      WarningTracker::printWarning((std::string) "Warning: The orbital with index " + std::to_string(iOrb + 1) + "\n" +
                                       "         is not orthogonal within the given basis set!\n" +
                                       "         Absolute sum of off-diagonal overlap elements:\n" + "         " +
                                       std::to_string(check),
                                   true);
    }
    if (std::fabs(s_mo(iOrb, iOrb) - 1) > 1e-6) {
      WarningTracker::printWarning((std::string) "Warning: The orbital with index " + std::to_string(iOrb + 1) + "\n" +
                                       "         is not normalized within the given basis set!\n" +
                                       "         Diagonal overlap element: " + std::to_string(s_mo(iOrb, iOrb)),
                                   true);
    }
  } // for iOrb
}

template<Options::SCF_MODES SCFMode>
std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> OrbitalsIOTask<SCFMode>::_sortMatrices = {
    {0, std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Identity(1, 1))},
    {1, nullptr},
    {2, nullptr},
    {3, nullptr},
    {4, nullptr},
    {5, nullptr},
    {6, nullptr},
    {7, nullptr},
    {8, nullptr},
    {9, nullptr},
};
template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::initTurboSortingMatrices() {
  /*
   * These matrices were generated by comparing Hartree-Fock orbital
   * coefficients between Serenity and Turbomole.
   */
  // Adjust according to the supported value of l.
  _maxL = 4;
  auto p = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(3, 3));
  // clang-format off
  // Turbomole ordering pz px py
  *p << 0, 1, 0,
        0, 0, 1,
        1, 0, 0;
  _sortMatrices[1] = p;
  auto d = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(5, 5));
  *d << 0, 0, 0, 1, 0,
        0, 0, 1, 0, 0,
        1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 0, 0, 1;
  _sortMatrices[2] = d;
  auto f = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(7, 7));
  *f << 0, 0, 0, 0, 0, 0,-1,
        0, 0, 0, 1, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 1, 0;
  _sortMatrices[3] = f;
  auto g = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(9, 9));
  *g << 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0,-1, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0,-1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1;
  _sortMatrices[4] = g;
  // clang-format on
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::initMolproSortingMatrices() {
  /*
   * These matrices were generated by comparing Hartree-Fock orbital
   * coefficients between Serenity and Molpro.
   */
  // Adjust according to the supported value of l.
  _maxL = 4;
  auto p = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(3, 3));
  // clang-format off
  *p << 0, 1, 0,
        0, 0, 1,
        1, 0, 0;
  _sortMatrices[1] = p;
  auto d = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(5, 5));
  *d << 0, 1, 0, 0, 0,
        0, 0, 0, 0, 1,
        1, 0, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0;
  _sortMatrices[2] = d;
  auto f = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(7, 7));
  *f << 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0, 0;
  _sortMatrices[3] = f;
  auto g = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(9, 9));
  *g << 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0;
  _sortMatrices[4] = g;
  // clang-format on
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::initSphericalMoldenSortingMatrices() {
  /*
   * These matrices were generated by comparing Hartree-Fock orbital
   * coefficients between Serenity and Molden.
   */
  // Adjust according to the supported value of l.
  _maxL = 4;
  auto p = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(3, 3));
  // clang-format off
  *p << 0, 1, 0,
        0, 0, 1,
        1, 0, 0;
  _sortMatrices[1] = p;
  auto d = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(5, 5));
  *d << 0, 0, 0, 0, 1,
        0, 0, 1, 0, 0,
        1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 0, 1, 0;
  _sortMatrices[2] = d;
  auto f = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(7, 7));
  *f << 0, 0, 0, 0, 0, 0,-1,
        0, 0, 0, 0, 1, 0, 0,
        0, 0, 1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0,-1, 0;
  _sortMatrices[3] = f;
  auto g = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(9, 9));
  *g << 0, 0, 0, 0, 0, 0, 0, 0,-1,
        0, 0, 0, 0, 0, 0,-1, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,-1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,-1, 0;
  _sortMatrices[4] = g;
  // clang-format on
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::initCartesianMoldenSortingMatrices() {
  /*
   * These matrices were generated by comparing Hartree-Fock orbital
   * coefficients between Serenity and Molden.
   */
  // Adjust according to the supported value of l.
  _maxL = 4;
  auto p = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(3, 3));
  // clang-format off
  *p << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;
  _sortMatrices[1] = p;
  auto d = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(6, 6));
  *d << 1, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 1, 0, 0, 0;
  _sortMatrices[2] = d;
  auto f = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(10, 10));
  *f << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
  _sortMatrices[3] = f;
  auto g = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(15, 15));
  *g << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  _sortMatrices[4] = g;
  // clang-format on
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::readMolproXMLOrbitals(CoefficientMatrix<SCFMode>& coeffs,
                                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues) {
  std::ifstream input(settings.path);
  if (SCFMode == UNRESTRICTED)
    throw SerenityError("ERROR: The reading of Molpro-xml files is only supported for RESTRICTED calculations");
  const unsigned int nBFs = coeffs.getBasisController()->getNBasisFunctions();
  this->skipToOrbitalDefinitionMolproXML(input, settings.path);
  bool orbitalDefinitionEnd = false;
  unsigned int iOrb = 0;
  while (true) {
    if (input.peek() == EOF)
      break;
    const double eigenvalue = this->getMolproXMLOrbitalEigenvalue(input, orbitalDefinitionEnd);
    if (orbitalDefinitionEnd)
      break;
    const Eigen::VectorXd coefficients = this->getMolproXMLOrbitalCoefficients(input, nBFs);
    this->assignCoefficentsAndEigenvalue(coefficients, coeffs, eigenvalues, eigenvalue, true, iOrb);
    ++iOrb;
  }
  if (!orbitalDefinitionEnd)
    throw SerenityError("ERROR: Orbital definition in Molpro-xml file ended unexpectedly.");
  for_spin(coeffs) {
    coeffs_spin = this->resortCoefficients(coeffs_spin);
    checkNorm(coeffs_spin);
  };
}

template<>
void OrbitalsIOTask<RESTRICTED>::assignCoefficentsAndEigenvalue(const Eigen::VectorXd& orbCoeff,
                                                                CoefficientMatrix<RESTRICTED>& coeffMatrix,
                                                                SpinPolarizedData<RESTRICTED, Eigen::VectorXd>& eigenvalues,
                                                                const double eigenvalue, const bool isAlpha,
                                                                const unsigned int iOrb) {
  (void)isAlpha;
  eigenvalues(iOrb) = eigenvalue;
  coeffMatrix.col(iOrb) = orbCoeff;
}

template<>
void OrbitalsIOTask<UNRESTRICTED>::assignCoefficentsAndEigenvalue(const Eigen::VectorXd& orbCoeff,
                                                                  CoefficientMatrix<UNRESTRICTED>& coeffMatrix,
                                                                  SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>& eigenvalues,
                                                                  const double eigenvalue, const bool isAlpha,
                                                                  const unsigned int iOrb) {
  if (isAlpha) {
    eigenvalues.alpha(iOrb) = eigenvalue;
    coeffMatrix.alpha.col(iOrb) = orbCoeff;
  }
  else {
    eigenvalues.beta(iOrb) = eigenvalue;
    coeffMatrix.beta.col(iOrb) = orbCoeff;
  }
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::skipToOrbitalDefinitionMolproXML(std::ifstream& input, std::string& filePath) {
  std::string line, word;
  bool moDefinitionFound = false;
  try {
    while (std::getline(input, line)) {
      std::istringstream iss(line);
      iss >> word;
      if (word == "<orbitals") {
        moDefinitionFound = true;
        break;
      }
    }
    if (moDefinitionFound) {
      std::string angularMomentumMode, spinMode, dummy, symmetry;
      std::istringstream iss(line);
      iss >> dummy;
      iss >> dummy;
      iss >> angularMomentumMode;
      iss >> spinMode;
      if (angularMomentumMode != "angular=\"spherical\"" && this->_system->getBasisController()->isPureSpherical())
        std::cout << "ERROR: The system definition expects a spherical basis set which appears to not be encoded in "
                     "the XML file."
                  << std::endl;

      if (angularMomentumMode != "angular=\"cartesian\"" && this->_system->getBasisController()->isPureCartesian())
        std::cout << "ERROR: The system definition expects a Cartesian basis set which appears to not be encoded in "
                     "the XML file."
                  << std::endl;

      if (spinMode != "spin=\"closed\"" && SCFMode == RESTRICTED)
        std::cout << "ERROR: The system definition expects a closed-shell orbital set which appears to not be encoded "
                     "in the XML file."
                  << std::endl;

      if (spinMode != "spin=\"open\"" && SCFMode == UNRESTRICTED)
        std::cout << "ERROR: The system definition expects an open-shell orbital set which appears to not be encoded "
                     "in the XML file."
                  << std::endl;

      std::getline(input, line);
      std::istringstream iss2(line);
      iss2 >> dummy;
      iss2 >> symmetry;
      if (symmetry != "state_symmetry=\"1\"")
        std::cout << "ERROR: Serenity does not support any symmetry. However, the orbitals encoded in the file appear "
                     "to use symmetry."
                  << std::endl;
      // Expects <orbitals ... > in either two or three lines
      if (line.substr(line.length() - 1) != ">")
        std::getline(input, line);
    }
  }
  catch (...) {
    throw SerenityError("ERROR: Unknown file format in Molpro-xml file: " + filePath);
  }

  if (!moDefinitionFound) {
    throw SerenityError("ERROR: No molecular orbitals defined in Molpro-xml file: " + filePath);
  }
}

template<Options::SCF_MODES SCFMode>
double OrbitalsIOTask<SCFMode>::getMolproXMLOrbitalEigenvalue(std::ifstream& input, bool& orbitalDefinitionEnd) {
  /*
   * Expected format:
   * <orbital occupation="2.0" energy="-20.546221" symmetryID="1">
   */
  std::string line, orbital, occupation, eigenvalueDefinition, eigenvalueString;
  std::getline(input, line);
  std::istringstream iss(line);
  iss >> orbital;
  if (orbital == "</orbitals>") {
    orbitalDefinitionEnd = true;
    return 0.0;
  }
  orbitalDefinitionEnd = false;
  if (orbital != "<orbital")
    throw SerenityError("ERROR: Unexpected file format in Molpro-xml file for orbital definition. Expected: '<orbital' "
                        "at the begging of the line, but got: " +
                        orbital);
  iss >> occupation;
  iss >> eigenvalueDefinition;
  iss = std::istringstream(eigenvalueDefinition);
  std::getline(iss, eigenvalueString, '=');
  std::getline(iss, eigenvalueString, '=');
  eigenvalueString.erase(0, 1);
  eigenvalueString.erase(eigenvalueString.size() - 1);
  double eigenvalue = std::stod(eigenvalueString);
  return eigenvalue;
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd OrbitalsIOTask<SCFMode>::getMolproXMLOrbitalCoefficients(std::ifstream& input, unsigned int nBFs) {
  std::string line, coefficientString, test;
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(nBFs);
  unsigned int iBF = 0;
  try {
    while (true) {
      std::getline(input, line);
      std::istringstream iss(line);
      std::istringstream iss2(line);
      iss2 >> test;
      if (test == "</orbital>")
        break;
      while (iss >> coefficientString) {
        coeffs(iBF) = std::stod(coefficientString);
        iBF++;
      }
    } // while iBF < nBFs
  }
  catch (...) {
    throw SerenityError("ERROR: Unexpected file format in Molpro-xml file for orbital definition.");
  }
  if (iBF < nBFs)
    throw SerenityError("ERROR: The basis set in the Molpro-xml file has a different dimension than expected! The "
                        "basis set may be incorrect!");

  return coeffs;
}

/*****************************************************/
/**************** Molcas HDF5 reading ****************/
/*****************************************************/

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::readMolcasHDF5Orbitals(CoefficientMatrix<SCFMode>& coeffs,
                                                     SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues) {
  const std::string filePath = this->settings.path + "/" + _system->getSystemName() + ".scf.h5";
  this->checkMolcasBasisDefinition(filePath);
  this->readMolcasEigenvalues(filePath, eigenvalues);
  this->readMolcasCoefficients(filePath, coeffs);
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::checkMolcasBasisDefinition(const std::string& filePath) {
  const std::string molcasBasFuncLabel = "NBAS";
  auto basisController = _system->getBasisController();
  const unsigned int nBasFuncSerenity = basisController->getNBasisFunctions();
  HDF5::Filepath name(filePath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::attribute_exists(file, molcasBasFuncLabel);
  Eigen::MatrixXi tmp = Eigen::MatrixXi::Zero(1, 1);
  EigenHDF5::load_attribute(file, molcasBasFuncLabel, tmp);
  const unsigned int nBasFuncMolcas = tmp(0, 0);
  file.close();
  if (nBasFuncMolcas != nBasFuncSerenity) {
    throw SerenityError("Different number of basis functions in Serenity's basis definition and molcas basis\n"
                        "definition. Make sure you are using the same basis set and spherical harmonics\n"
                        "(spherical catesian) settings.");
  }
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::readMolcasEigenvalues(const std::string& filePath,
                                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues) {
  const std::string molcasMOEigenvalueLabel = "MO_ENERGIES";
  HDF5::Filepath name(filePath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, molcasMOEigenvalueLabel);
  for_spin(eigenvalues) {
    HDF5::load(file, molcasMOEigenvalueLabel, eigenvalues_spin);
  };
  file.close();
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::readMolcasCoefficients(const std::string& filePath, CoefficientMatrix<SCFMode>& coeffs) {
  const std::string molcasMOCoefficientLabel = "MO_VECTORS";
  HDF5::Filepath name(filePath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, molcasMOCoefficientLabel);
  unsigned int nBasFunc = _system->getBasisController()->getNBasisFunctions();
  const auto sortingMatrix = this->getMolcasSortingMatrix(filePath);
  for_spin(coeffs) {
    Eigen::VectorXd buffer = Eigen::VectorXd::Zero(nBasFunc * nBasFunc);
    HDF5::load(file, molcasMOCoefficientLabel, buffer);
    const Eigen::MatrixXd unsortedCoefficients = Eigen::Map<Eigen::MatrixXd>(buffer.data(), nBasFunc, nBasFunc);
    coeffs_spin = sortingMatrix.template cast<double>() * unsortedCoefficients;
  };
  file.close();
}
template<Options::SCF_MODES SCFMode>
Eigen::MatrixXi OrbitalsIOTask<SCFMode>::getMolcasSortingMatrix(const std::string& filePath) {
  const std::string basisFunctionSortingLabel = "BASIS_FUNCTION_IDS";
  const unsigned int nBasFunc = _system->getBasisController()->getNBasisFunctions();
  HDF5::Filepath name(filePath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, basisFunctionSortingLabel);
  // The raw index matrix from molcas. Ordering: Atom-center, shell-type index on center, angular momentum, ang.
  // orientation. E.g.
  //  1, 1, 0, 0,   --> Atom1, 1s
  //  1, 2, 0, 0,   --> Atom1, 2s
  //  1, 1, 1, -1,  --> Atom1, first p, m=-1
  Eigen::MatrixXi rawMatrix = Eigen::MatrixXi::Zero(nBasFunc, 4);
  EigenHDF5::load(file, basisFunctionSortingLabel, rawMatrix);
  if (rawMatrix.rows() != nBasFunc)
    throw SerenityError("Inconsistent orbital indexing in the molcas orbital file.");
  Eigen::MatrixXi sortingMatrix = Eigen::MatrixXi::Zero(nBasFunc, nBasFunc);
  for (unsigned int iRow = 0; iRow < rawMatrix.rows(); ++iRow) {
    const int atomIndex = rawMatrix(iRow, 0) - 1;
    if (atomIndex < 0)
      throw SerenityError("Expected indexing starting with 1, but found an atom with index 0.");
    const unsigned int basisFunctionIndex =
        this->getBasisFunctionIndex(atomIndex, rawMatrix(iRow, 2), rawMatrix(iRow, 1), rawMatrix(iRow, 3));
    sortingMatrix(basisFunctionIndex, iRow) = 1;
  }
  return sortingMatrix;
}
template<Options::SCF_MODES SCFMode>
unsigned int OrbitalsIOTask<SCFMode>::getBasisFunctionIndex(const unsigned int atomIndex, const unsigned int angularMomentum,
                                                            const unsigned int iShellOfMomentum, const int orientation) {
  const auto atomCenteredBasisController = _system->getAtomCenteredBasisController();
  const auto shellIndices = atomCenteredBasisController->getBasisIndicesRed()[atomIndex];
  const auto shells = atomCenteredBasisController->getBasis();
  unsigned int counter = 0;
  unsigned int shellIndex = shells.size();
  for (unsigned int iShell = shellIndices.first; iShell < shellIndices.second; ++iShell) {
    const auto& shell = shells[iShell];
    if (shell->getAngularMomentum() == angularMomentum)
      counter++;
    if (counter == iShellOfMomentum) {
      shellIndex = iShell;
      break;
    }
  }
  if (shellIndex >= shells.size())
    throw SerenityError((std::string) "Inconsistent indexing for shells on atoms in the molcas orbital file. Atom " +
                        std::to_string(atomIndex) + " has no shell with\nangular momentum " +
                        std::to_string(angularMomentum) + " and index " + std::to_string(iShellOfMomentum) + ".");
  // the index will be firstFuncIndex + l + m
  // E.g. l=2, m=1 --> n + 2 + 1 = n + 3
  // for l=2 the complete shell is ordered as -2 -1 0 1 2.
  // In this array, m=1 has the index 3.
  return atomCenteredBasisController->extendedIndex(shellIndex) + angularMomentum + orientation;
}

template<Options::SCF_MODES SCFMode>
void OrbitalsIOTask<SCFMode>::replaceMolcasHDF5Orbitals(const CoefficientMatrix<SCFMode>& coeffs) {
  const std::string backupFilePath = this->settings.path + "/" + _system->getSystemName() + "_backup.scf.h5";
  const std::string filePath = this->settings.path + "/" + _system->getSystemName() + ".scf.h5";
  // Copy the original file
  {
    std::ifstream src(filePath, std::ios::binary);
    std::ofstream dst(backupFilePath, std::ios::binary);
    dst << src.rdbuf();
    dst.flush();
  }
  this->checkMolcasBasisDefinition(filePath);
  const auto sortingMatrix = this->getMolcasSortingMatrix(filePath);
  const std::string molcasMOCoefficientLabel = "MO_VECTORS";
  // Replace in the copied file.
  HDF5::Filepath name(filePath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDWR);
  HDF5::dataset_exists(file, molcasMOCoefficientLabel);
  for_spin(coeffs) {
    Eigen::MatrixXd tmp = sortingMatrix.transpose().template cast<double>() * coeffs_spin;
    Eigen::VectorXd buffer = Eigen::Map<Eigen::VectorXd>(tmp.data(), tmp.rows() * tmp.cols());
    HDF5::replace(file, molcasMOCoefficientLabel, buffer);
  };
  file.close();
}

template<Options::SCF_MODES SCFMode>
OrbitalsIOTask<SCFMode>::~OrbitalsIOTask() = default;

template class OrbitalsIOTask<Options::SCF_MODES::RESTRICTED>;
template class OrbitalsIOTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
