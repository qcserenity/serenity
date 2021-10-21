/**
 * @file ReadOrbitalsTask.cpp
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
#include "tasks/ReadOrbitalsTask.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"                             //Sehll-wise sorting.
#include "basis/BasisController.h"                   //Basis.
#include "data/ElectronicStructure.h"                //Constructor.
#include "data/OrbitalController.h"                  //Constructor.
#include "integrals/OneElectronIntegralController.h" //Overlap integrals.
#include "io/HDF5.h"                                 //Load from HDF5.
#include "math/linearAlgebra/MatrixFunctions.h"      //Cholesky orthogonalization.
#include "misc/SerenityError.h"                      //Errors.
#include "misc/WarningTracker.h"                     //Warnings.
#include "system/SystemController.h"                 //Access to system properties.
/* Include Std and External Headers */
#include <algorithm> //Replace within a string.
#include <fstream>   //input/output streams.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ReadOrbitalsTask<SCFMode>::ReadOrbitalsTask(std::shared_ptr<SystemController> system) : _system(system) {
}

template<Options::SCF_MODES SCFMode>
void ReadOrbitalsTask<SCFMode>::run() {
  bool updateElecStructure = false;
  auto basisController = _system->getBasisController();
  auto finalCoeffPtr = std::make_unique<CoefficientMatrix<SCFMode>>(basisController);
  auto eigenvaluesPtr = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(
      Eigen::VectorXd::Zero(basisController->getNBasisFunctions()));
  if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::TURBOMOLE) {
    updateElecStructure = true;
    initTurboSortingMatrices();
    checkTurboInput();
    readTurbomoleOrbitals(*finalCoeffPtr, *eigenvaluesPtr);
  }
  if (updateElecStructure) {
    auto newOrbitalController = std::make_shared<OrbitalController<SCFMode>>(
        std::move(finalCoeffPtr), basisController, std::move(eigenvaluesPtr), _system->getNCoreElectrons());
    auto targetElectronicStructure =
        std::make_shared<ElectronicStructure<SCFMode>>(newOrbitalController, _system->getOneElectronIntegralController(),
                                                       _system->template getNOccupiedOrbitals<SCFMode>());
    _system->template setElectronicStructure<SCFMode>(targetElectronicStructure);
    targetElectronicStructure->toHDF5(_system->getHDF5BaseName(), _system->getSystemIdentifier());
  }
  if (settings.fileFormat == Options::ORBITAL_FILE_TYPES::SERENITY) {
    std::string systemId = getSerenityIDFromFile();
    auto targetElectronicStructure = std::make_shared<ElectronicStructure<SCFMode>>(
        settings.path + "/" + _system->getSystemName(), _system->getBasisController(), _system->getGeometry(), systemId);
    _system->template setElectronicStructure<SCFMode>(targetElectronicStructure);
    if (settings.resetCoreOrbitals)
      _system->template getActiveOrbitalController<SCFMode>()->setCoreOrbitalsByNumber(_system->getNCoreElectrons() / 2);
    targetElectronicStructure->toHDF5(_system->getHDF5BaseName(), _system->getSystemIdentifier());
  }
}

template<Options::SCF_MODES SCFMode>
std::string ReadOrbitalsTask<SCFMode>::getSerenityIDFromFile() {
  std::string filePath = settings.path + "/" + _system->getSystemName();
  if (SCFMode == RESTRICTED) {
    filePath += ".orbs.res.h5";
  }
  else {
    filePath += ".orbs.unres.h5";
  }
  HDF5::Filepath name(filePath);
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5::attribute_exists(file, "ID");
  std::string id = HDF5::load_attribute<std::string>(file, "ID");
  file.close();
  return id;
}

template<Options::SCF_MODES SCFMode>
void ReadOrbitalsTask<SCFMode>::checkTurboInput() {
  auto basisController = _system->getBasisController();
  if (basisController->getMaxAngularMomentum() > _maxL)
    throw SerenityError((std::string) "ERROR: The maximum angular momentum for which Turbomole sorting matrices are\n" +
                        "       available was exceeded. Currently supported l= " + std::to_string(_maxL));
  if (not basisController->isPureSpherical()) {
    throw SerenityError((std::string) "ERROR: Only spherical harmonics are supported when reading Turbomole orbitals.");
  }
}

template<>
void ReadOrbitalsTask<RESTRICTED>::readTurbomoleOrbitals(CoefficientMatrix<RESTRICTED>& coeffs,
                                                         SpinPolarizedData<RESTRICTED, Eigen::VectorXd>& eigenvalues) {
  unsigned int nBFs = coeffs.getBasisController()->getNBasisFunctions();
  auto coeff_and_eps = readTurbomoleOrbitalFile(settings.path + "/mos", nBFs);
  coeffs.block(0, 0, nBFs, nBFs) = coeff_and_eps.first;
  eigenvalues = coeff_and_eps.second;
}

template<>
void ReadOrbitalsTask<UNRESTRICTED>::readTurbomoleOrbitals(CoefficientMatrix<UNRESTRICTED>& coeffs,
                                                           SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>& eigenvalues) {
  unsigned int nBFs = coeffs.getBasisController()->getNBasisFunctions();
  auto coeff_and_eps_alpha = readTurbomoleOrbitalFile(settings.path + "/alpha", nBFs);
  coeffs.alpha.block(0, 0, nBFs, nBFs) = coeff_and_eps_alpha.first;
  eigenvalues.alpha = coeff_and_eps_alpha.second;
  auto coeff_and_eps_beta = readTurbomoleOrbitalFile(settings.path + "/beta", nBFs);
  coeffs.beta.block(0, 0, nBFs, nBFs) = coeff_and_eps_beta.first;
  eigenvalues.beta = coeff_and_eps_beta.second;
}

template<Options::SCF_MODES SCFMode>
double ReadOrbitalsTask<SCFMode>::getTurboEigenvalue(std::istringstream& iss) {
  std::string word;
  iss >> word;
  std::replace(word.begin(), word.end(), 'D', 'e');
  return std::stod(word.substr(11));
}
template<Options::SCF_MODES SCFMode>
void ReadOrbitalsTask<SCFMode>::checkTurboNAOs(std::istringstream& iss, unsigned int nBFs) {
  std::string word;
  iss >> word;
  unsigned int nAOs = std::stoi(word.substr(6));
  if (nAOs != nBFs) {
    throw SerenityError((std::string) "ERROR: Inconsistent number of basis functions in file!");
  }
}
template<Options::SCF_MODES SCFMode>
Eigen::VectorXd ReadOrbitalsTask<SCFMode>::getTurboOrbital(std::ifstream& input, unsigned int nBFs,
                                                           unsigned int nPerLine, unsigned int length) {
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
void ReadOrbitalsTask<SCFMode>::splitTurboOrbLine(std::string& line, unsigned int nPerLine, unsigned int length,
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
void ReadOrbitalsTask<SCFMode>::resortTurboCoefficients(Eigen::MatrixXd& coefficients) {
  auto basis = _system->getBasisController()->getBasis();
  unsigned int startRow = 0;
  unsigned int nOrbs = coefficients.cols();
  for (const auto& shell : basis) {
    unsigned int angularMomentum = shell->getAngularMomentum();
    unsigned int nContracted = shell->getNContracted();
    const Eigen::MatrixXd& rotationMatrix = *_turboSortMatrices[angularMomentum];
    coefficients.block(startRow, 0, nContracted, nOrbs) =
        (rotationMatrix * coefficients.block(startRow, 0, nContracted, nOrbs)).eval();
    startRow += nContracted;
  } // for shell
}

template<Options::SCF_MODES SCFMode>
std::pair<Eigen::MatrixXd, Eigen::VectorXd> ReadOrbitalsTask<SCFMode>::readTurbomoleOrbitalFile(std::string filePath,
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
  resortTurboCoefficients(coefficients);
  checkNorm(coefficients);
  coefficients = orthogonalize_chol(coefficients, _system->getOneElectronIntegralController()->getOverlapIntegrals());
  return std::make_pair(coefficients, eigenvalues);
}

template<Options::SCF_MODES SCFMode>
void ReadOrbitalsTask<SCFMode>::checkNorm(const Eigen::MatrixXd& coefficients) {
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
std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> ReadOrbitalsTask<SCFMode>::_turboSortMatrices = {
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
void ReadOrbitalsTask<SCFMode>::initTurboSortingMatrices() {
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
  _turboSortMatrices[1] = p;
  auto d = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(5, 5));
  *d << 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0,
      1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 0, 0, 1;
  _turboSortMatrices[2] = d;
  auto f = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(7, 7));
  *f << 0, 0, 0, 0, 0, 0,-1,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0;
  _turboSortMatrices[3] = f;
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
  _turboSortMatrices[4] = g;
  // clang-format on
}

template<Options::SCF_MODES SCFMode>
ReadOrbitalsTask<SCFMode>::~ReadOrbitalsTask() = default;

template class ReadOrbitalsTask<Options::SCF_MODES::RESTRICTED>;
template class ReadOrbitalsTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
