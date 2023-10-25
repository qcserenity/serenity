/**
 * @file OrbitalsIOTask.h
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

#ifndef TASKS_ORBITALSIOTASK_H_
#define TASKS_ORBITALSIOTASK_H_
/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <string>

namespace Serenity {
/* Forward declarations */
class SystemController;
class Geometry;
using namespace Serenity::Reflection;
struct OrbitalsIOTaskSettings {
  OrbitalsIOTaskSettings()
    : fileFormat(Options::ORBITAL_FILE_TYPES::TURBOMOLE), path("."), resetCoreOrbitals(false), replaceInFile(false), write(false) {
  }
  REFLECTABLE((Options::ORBITAL_FILE_TYPES)fileFormat, (std::string)path, (bool)resetCoreOrbitals, (bool)replaceInFile,
              (bool)write)
};

/**
 * @class
 * @brief Reads orbitals from external files and assigns them to the given system. Can also write files.
 *
 * Currently supported:\n
 *    - turbomole ASCII-MO files, spherical harmonics.\n
 *    - serenity HDF5 files.
 *    - molpro xml files.
 *    - molcas HDF5 files.
 *    - Molden files for both spherical and cartesian harmonics.
 */
template<Options::SCF_MODES SCFMode>
class OrbitalsIOTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param system The system controller.
   */
  OrbitalsIOTask(std::shared_ptr<SystemController> system);
  /**
   * @brief Default destructor.
   */
  virtual ~OrbitalsIOTask();
  /**
   * @brief Execute the task.
   */
  void run();
  /**
   * @brief The settings.
   *   -fileformat:        The file format of the MO files.
   *   -path:              Path to the directory containing the file to be read.
   *   -resetCoreOrbitals: Reset the core orbital assignment according to the geometry and eigenvalues.
   *   -replaceInFile:     If true, copies the original orbital file and replaces the coefficients in this file
   *                       by Serenity's coefficients. (Only implemented for Molcas HDF5 files).
   *   -write              Write files in either Turbomole or Molden format.
   */
  OrbitalsIOTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _system;
  /*
   * Read the coefficients and eigenvalues from a turbomole file to the given parameters.
   * Takes care of RESTRICTED vs. UNRESTRICTED
   */
  void readTurbomoleOrbitals(CoefficientMatrix<SCFMode>& coeffs, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);
  /*
   * Write the coefficients and eigenvalues from Serenity to a Turbomole File.
   * Takes care of RESTRICTED vs. UNRESTRICTED
   */
  void writeTurbomoleOrbitals(CoefficientMatrix<SCFMode>& coefficients, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);
  /*
   * Write a Molden File.
   */
  void writeMoldenOrbitals();
  /*
   * Read a Turbomole orbital file.
   */
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> readTurbomoleOrbitalFile(std::string filePath, unsigned int nBFs);
  /*
   * Write a Turbomole orbital file.
   */
  void writeTurbomoleOrbitalFile(std::string filePath, unsigned int nBFs, const Eigen::VectorXd& eigenvalues,
                                 const Eigen::MatrixXd& coefficients);
  /**
   * Get the next eigenvalue.
   */
  double getTurboEigenvalue(std::istringstream& iss);
  /*
   * Check the number of AOs correspond to the number of coefficients for each orbital.
   */
  void checkTurboNAOs(std::istringstream& iss, unsigned int nBFs);
  /*
   * Get the next orbital from the file.
   */
  Eigen::VectorXd getTurboOrbital(std::ifstream& input, unsigned int nBFs, unsigned int nPerLine, unsigned int length);
  /*
   * Get the next set of coefficients from the line.
   */
  void splitTurboOrbLine(std::string& line, unsigned int nPerLine, unsigned int length, Eigen::VectorXd& coeffs,
                         unsigned int& counter);
  /*
   * Sort the turbomole coefficients to fit the basis function ordering in Serenity.
   */
  Eigen::MatrixXd resortCoefficients(const Eigen::MatrixXd& coefficients);
  /*
   * Reformat Number to ASCII format.
   */
  std::string reformatNumber(double number, std::string exponentSymbol = "D");
  /*
   * Initialize sorting matrices for Turbomole.
   */
  void initTurboSortingMatrices();
  /*
   * Initialize spherical sorting matrices for Molden.
   */
  void initSphericalMoldenSortingMatrices();
  /*
   * Initialize cartesian sorting matrices for Molden.
   */
  void initCartesianMoldenSortingMatrices();
  /*
   * Check if the task input allows reading of turbomole orbitals.
   */
  void checkTurboInput();
  /*
   * Check if the orbitals are normalized and orthogonal with respect to the overlap matrix in
   * Serenity.
   */
  void checkNorm(const Eigen::MatrixXd& coefficients);
  /*
   * The map that contains the sorting matrices for turbomole orbitals..
   */
  static std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> _sortMatrices;
  /*
   * The maximum angular momentum tolerated for turbomole orbitals. This is overridden in
   * initTurboSortingMatrices().
   */
  unsigned int _maxL = 0;

  std::string getSerenityIDFromFile();

  /**
   * **** Molpro xml files ****
   */
  void readMolproXMLOrbitals(CoefficientMatrix<SCFMode>& coeffs, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);
  /*
   * Get the next set of orbital coefficients form the input stream.
   */
  Eigen::VectorXd getMolproXMLOrbitalCoefficients(std::ifstream& input, unsigned int nBFs);
  /*
   * Get the next orbital eigenvalue from file.
   */
  double getMolproXMLOrbitalEigenvalue(std::ifstream& input, bool& orbitalDefinitionEnd);
  /*
   * Skip the stream to the orbital definitions.
   */
  void skipToOrbitalDefinitionMolproXML(std::ifstream& input, std::string& filePath);
  /*
   * Assign coefficents and eigenvalues to the coefficient matrix and eigenvalue vector (according to spin settings)
   */
  void assignCoefficentsAndEigenvalue(const Eigen::VectorXd& orbCoeff, CoefficientMatrix<SCFMode>& coeffMatrix,
                                      SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues, const double eigenvalue,
                                      const bool isAlpha, const unsigned int iOrb);
  /*
   * Initialize sorting matrices for Molpro.
   */
  void initMolproSortingMatrices();
  /*
   * Check the input for reading Molpro orbitals from xml files.
   */
  void checkMolproInput();

  /*
   * Read orbitals from molcas HDF5-files.
   */
  void readMolcasHDF5Orbitals(CoefficientMatrix<SCFMode>& coeffs, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);
  /*
   * Replace the orbitals defined in a Molcas HDF5 file by Serenity's cofficients.
   */
  void replaceMolcasHDF5Orbitals(const CoefficientMatrix<SCFMode>& coeffs);
  /*
   * Ensure that the data sets have the same dimensions in Serenity and Molcas.
   */
  void checkMolcasBasisDefinition(const std::string& filePath);
  /*
   * Read the orbital eigenvalues.
   */
  void readMolcasEigenvalues(const std::string& filePath, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);
  /*
   * Read and sort the orbital coefficients.
   */
  void readMolcasCoefficients(const std::string& filePath, CoefficientMatrix<SCFMode>& coeffs);
  /*
   * Get the permutation matrix for the AO ordering.
   */
  Eigen::MatrixXi getMolcasSortingMatrix(const std::string& filePath);
  /*
   * Get the basis function index characterized by its atom-index, angular momentum, ang. orientation, and index of the
   * shell with the given angular momentum on the atom.
   */
  unsigned int getBasisFunctionIndex(const unsigned int atomIndex, const unsigned int angularMomentum,
                                     const unsigned int iShellOfMomentum, const int orientation);
};

} /* namespace Serenity */

#endif /* TASKS_ORBITALSIOTASK_H_ */
