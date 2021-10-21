/**
 * @file ReadOrbitalsTask.h
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

#ifndef TASKS_READORBITALSTASK_H_
#define TASKS_READORBITALSTASK_H_
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
using namespace Serenity::Reflection;
struct ReadOrbitalsTaskSettings {
  ReadOrbitalsTaskSettings() : fileFormat(Options::ORBITAL_FILE_TYPES::TURBOMOLE), path("."), resetCoreOrbitals(false) {
  }
  REFLECTABLE((Options::ORBITAL_FILE_TYPES)fileFormat, (std::string)path, (bool)resetCoreOrbitals)
};

/**
 * @class
 * @brief Reads orbitals from external files and assigns them to the given system.
 *
 * Currently supported:\n
 *    - turbomole ASCII-MO files, spherical harmonics.\n
 *    - serenity HDF5 files.
 */
template<Options::SCF_MODES SCFMode>
class ReadOrbitalsTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param system The system controller.
   */
  ReadOrbitalsTask(std::shared_ptr<SystemController> system);
  /**
   * @brief Default destructor.
   */
  virtual ~ReadOrbitalsTask();
  /**
   * @brief Execute the task.
   */
  void run();
  /**
   * @brief The settings.
   *   -turbomole:         Read ASCII turbomole MO files.
   *   -path:              Path to the directory containing the file to be read.
   *   -resetCoreOrbitals: Reset the core orbital assignement according to the geometry and eigenvalues.
   */
  ReadOrbitalsTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _system;
  /*
   * Read the coefficients and eigenvalues from a turbomole file to the given parameters.
   * Takes care of RESTRICTED vs. UNRESTRICTED
   */
  void readTurbomoleOrbitals(CoefficientMatrix<SCFMode>& coeffs, SpinPolarizedData<SCFMode, Eigen::VectorXd>& eigenvalues);
  /*
   * Read a turbomole orbital file.
   */
  std::pair<Eigen::MatrixXd, Eigen::VectorXd> readTurbomoleOrbitalFile(std::string filePath, unsigned int nBFs);
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
  void resortTurboCoefficients(Eigen::MatrixXd& coefficients);
  /*
   * Create all sorting matrices.
   */
  void initTurboSortingMatrices();
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
  static std::map<unsigned int, std::shared_ptr<Eigen::MatrixXd>> _turboSortMatrices;
  /*
   * The maximum angular momentum tolerated for turbomole orbitals. This is overridden in
   * initTurboSortingMatrices().
   */
  unsigned int _maxL = 0;

  std::string getSerenityIDFromFile();
};

} /* namespace Serenity */

#endif /* TASKS_READORBITALSTASK_H_ */
