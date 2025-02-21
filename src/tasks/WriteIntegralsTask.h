/**
 * @file   WriteIntegralsTask.h
 *
 * @date   Nov. 10, 2023
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

#ifndef SERENITY_WRITEINTEGRALSTASK_H
#define SERENITY_WRITEINTEGRALSTASK_H

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
template<Options::SCF_MODES>
class SPMatrix;

struct WriteIntegralsTaskSettings {
  WriteIntegralsTaskSettings() : hCoreIntegrals(false), fileFormat(Options::INTEGRAL_FILE_TYPES::HDF5) {
  }
  REFLECTABLE((bool)hCoreIntegrals, (Options::INTEGRAL_FILE_TYPES)fileFormat)
};

/**
 * @class WriteIntegralsTask WriteIntegralsTask.h
 * @brief A task which write integrals to file.
 */
class WriteIntegralsTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param system The system controller.
   */
  WriteIntegralsTask(std::shared_ptr<SystemController> system);
  /**
   * @brief Default destructor.
   */
  virtual ~WriteIntegralsTask() = default;
  /**
   * @brief Execute the task.
   */
  void run();
  /**
   * @brief The settings.
   *   -hCoreIntegrals     If True, the core hamiltonian integrals are written. By default, false.
   *   -fileFormat:        The file format for the integral files.
   */
  WriteIntegralsTaskSettings settings;

 private:
  std::shared_ptr<SystemController> _system;
  ///@brief Write HCore integrals to file.
  void writeHCoreIntegrals();
  ///@brief Write matrix to file.
  template<Options::SCF_MODES SCFMode>
  void writeFile(const std::string& fileBaseName, const SPMatrix<SCFMode>& integrals) const;
  ///@brief Write matrix to ASCII file.
  template<Options::SCF_MODES SCFMode>
  void writeAsciiFile(const std::string& fileBaseName, const SPMatrix<SCFMode>& integrals) const;
  ///@brief Write matrix to HDF5 file.
  template<Options::SCF_MODES SCFMode>
  void writeHDF5File(const std::string& fileBaseName, const SPMatrix<SCFMode>& integrals) const;
  ///@brief Get name tags for alpha and beta polarization.
  template<Options::SCF_MODES SCFMode>
  SpinPolarizedData<SCFMode, std::string> getFileSuffix() const;
};

} /* namespace Serenity */

#endif // SERENITY_WRITEINTEGRALSTASK_H
