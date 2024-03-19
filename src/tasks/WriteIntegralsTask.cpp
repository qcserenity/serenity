/**
 * @file   WriteIntegralsTask.cpp
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

/* Include Class Header*/
#include "tasks/WriteIntegralsTask.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutputStream.h"  // Filtered output streams.
#include "io/HDF5.h"                   // Write HDF5 files.
#include "potentials/HCorePotential.h" // Calculate HCore integrals.
#include "system/SystemController.h"   // Access system name.
/* Include Std and External Headers */
#include <iomanip> // Format ASCII output files.

namespace Serenity {

WriteIntegralsTask::WriteIntegralsTask(std::shared_ptr<SystemController> system) : _system(system) {
}

void WriteIntegralsTask::run() {
  if (this->settings.hCoreIntegrals) {
    this->writeHCoreIntegrals();
  }
}

void WriteIntegralsTask::writeHCoreIntegrals() {
  HCorePotential<RESTRICTED> hCorePotential(this->_system);
  const std::string baseName = "hcore";
  this->writeFile<RESTRICTED>(baseName, hCorePotential.getMatrix());
}

template<Options::SCF_MODES SCFMode>
void WriteIntegralsTask::writeFile(const std::string& fileBaseName, const SPMatrix<SCFMode>& integrals) const {
  switch (this->settings.fileFormat) {
    case Options::INTEGRAL_FILE_TYPES::HDF5: {
      this->writeHDF5File(fileBaseName, integrals);
      break;
    }
    case Options::INTEGRAL_FILE_TYPES::ASCII: {
      this->writeAsciiFile(fileBaseName, integrals);
      break;
    }
  }
}

template<Options::SCF_MODES SCFMode>
void WriteIntegralsTask::writeAsciiFile(const std::string& fileBaseName, const SPMatrix<SCFMode>& integrals) const {
  const std::string fileName = this->_system->getSystemName() + "." + fileBaseName + ".dat";
  OutputControl::nOut << "Writing integrals to file " << fileName << "\n" << std::endl;
  std::ofstream ofstream;
  ofstream.open(fileName);
  for_spin(integrals) {
    ofstream << std::scientific << std::setprecision(12) << integrals_spin << std::endl;
  };
  ofstream.close();
}

template<Options::SCF_MODES SCFMode>
void WriteIntegralsTask::writeHDF5File(const std::string& fileBaseName, const SPMatrix<SCFMode>& integrals) const {
  const auto suffix = this->getFileSuffix<SCFMode>();
  const std::string fileName = this->_system->getSystemName() + "." + fileBaseName + ".h5";
  OutputControl::nOut << "Writing integrals to file " << fileName << "\n" << std::endl;
  HDF5::H5File file(fileName.c_str(), H5F_ACC_TRUNC);
  for_spin(integrals, suffix) {
    HDF5::save(file, fileBaseName + suffix_spin, integrals_spin);
  };
  file.close();
}

template<>
SpinPolarizedData<RESTRICTED, std::string> WriteIntegralsTask::getFileSuffix() const {
  return "";
}

template<>
SpinPolarizedData<UNRESTRICTED, std::string> WriteIntegralsTask::getFileSuffix() const {
  SpinPolarizedData<UNRESTRICTED, std::string> suffix;
  suffix.alpha = ".alpha";
  suffix.beta = ".beta";
  return suffix;
}

} /* namespace Serenity */
