/**
 * @file HDF5.h
 *
 * @date Jun 8, 2016
 * @author Jan Unsleber
 *
 *
 * A header to bundle all the HDF5 support and generate a single Namespace for it.
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

#ifndef BASICS_IO_HDF5_H_
#define BASICS_IO_HDF5_H_

/* Include Serenity Internal Headers */
#include "io/Eigen3HDF5.h"
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <sys/stat.h>
#include <iostream>

namespace Serenity {
/**
 * @namespace HDF5 HDF5.h
 * @brief Namespace for all additional HDF5 related functionalities.
 */
namespace HDF5 {
using namespace EigenHDF5;
using namespace H5;

/**
 * @class Filepath HDF5.h
 * @brief A class that holds a HDF5 file file path.
 *
 * The reason this class exists is, that it is a central place to check for file existence.
 * Furthermore, functions can then require a file path which is way more descriptive than
 * a plain string.
 */
class Filepath : public std::string {
 public:
  /**
   * @brief Constructor.
   * @param path
   */
  Filepath(std::string path) : std::string(path) {
    // All HDF5 file from Serenity have this extension!
    assert(path.substr(path.length() - 3) == ".h5");
    // Check if the file exists
    struct stat buffer;
    if (stat(path.c_str(), &buffer) != 0) {
      throw SerenityError("File " + path + " not found");
    }
  }
  /**
   * @brief Default destructor.
   */
  virtual ~Filepath() = default;
};

/**
 * @brief Checks for a dataset in a HDF5 file.
 * @param file The HDF5 file.
 * @param dataSetName The needed dataset.
 */
inline void dataset_exists(H5::H5File file, const std::string& dataSetName) {
  Exception::dontPrint();
  try {
    DataSet test = file.openDataSet(dataSetName.c_str());
  }
  catch (GroupIException& not_found_error) {
    throw SerenityError("Dataset " + dataSetName + " in " + file.getFileName() + " not found");
  }
  catch (FileIException& not_found_error) {
    throw SerenityError("File " + file.getFileName() + " not found");
  }
}
/**
 * @brief Checks for an attribute in a HDF5 file.
 * @param file The HDF5 file.
 * @param dataSetName The needed dataset.
 */
inline void attribute_exists(H5::H5File file, const std::string& attributeName) {
  Exception::dontPrint();
  try {
    Attribute test = file.openAttribute(attributeName.c_str());
  }
  catch (AttributeIException& not_found_error) {
    throw SerenityError("Attribute " + attributeName + " in " + file.getFileName() + " not found");
  }
  catch (FileIException& not_found_error) {
    throw SerenityError("File " + file.getFileName() + " not found");
  }
}
/**
 * @brief Checks if an attribute in a HDF5 file matches the
 *        configuration of the System.
 * @param file The HDF5 file.
 * @param attributeName The name of the attribute to be checked.
 * @param reference The value within the current Configuration.
 */
template<typename T>
void check_attribute(H5::H5File file, const std::string& attributeName, const T& reference) {
  T value;
  load_scalar_attribute(file, attributeName, value);
  if (value != reference) {
    throw SerenityError(attributeName + " from file " + file.getFileName() + "does not match the current system");
  }
}
/**
 * @brief Load an attribute from HDF5 file.
 * @param file The HDF5 file.
 * @param attributeName The name of the attribute to be loaded.
 * @return The attribute value.
 */
template<typename T>
T load_attribute(H5::H5File file, const std::string& attributeName) {
  T value;
  load_scalar_attribute(file, attributeName, value);
  return value;
}

/**
 * @brief Loads a vector from file.
 * @param file The HDF5 file.
 * @param data_set_name The name of the dataset.
 * @param strings The vector to be filled with the data.
 */
inline void load_std_vector(H5::H5File file, const std::string& data_set_name, std::vector<std::string>& strings) {
  strings.resize(0);
  DataSet datset = file.openDataSet(data_set_name.c_str());
  DataSpace dataspace = getSpace(datset);

  const std::size_t ndims = dataspace.getSimpleExtentNdims();
  if (ndims != 1)
    throw SerenityError("Number of dimensions in HDF5 Dataspace does not match requirements.");
  hsize_t dimensions[1];
  dataspace.getSimpleExtentDims(dimensions);
  strings.resize(dimensions[0]);
  H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
  std::vector<char*> rdata(dimensions[0]);
  rdata.push_back(nullptr);
  datset.read(rdata.data(), datatype, dataspace);
  dataspace.close();
  datset.close();
  for (hsize_t ii = 0; ii < dimensions[0]; ii++) {
    strings[ii] = std::string(rdata[ii]);
  }
}

/**
 * @brief Saves a vector to file.
 * @param file The HDF5 file.
 * @param data_set_name The name of the dataset.
 * @param strings A vector of strings to be saved.
 */
inline void save_std_vector(H5::H5File file, const std::string& data_set_name, const std::vector<const char*>& strings) {
  //  one dimension
  hsize_t str_dimsf[1]{strings.size()};
  H5::DataSpace dataspace(1, str_dimsf);

  // Variable length string
  H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
  H5::DataSet str_dataset = file.createDataSet(data_set_name.c_str(), datatype, dataspace);

  str_dataset.write(strings.data(), datatype);
}

/**
 * @brief Overload for strings. See above.
 */
inline void save_std_vector(H5::H5File file, const std::string& data_set_name, const std::vector<std::string>& strings) {
  // HDF5 only understands vector of char*
  std::vector<const char*> arr_c_str;
  for (unsigned ii = 0; ii < strings.size(); ++ii) {
    arr_c_str.push_back(strings[ii].c_str());
  }

  save_std_vector(file, data_set_name, arr_c_str);
}

/**
 * @brief Saves an attribute to a HDF5 file. The attribute is
 *        converted to a string before writing it into file.
 * @param file The HDF5 file.
 * @param name The attribute name under which the attribute will
 *        be accessible within the file.
 * @param value The value the attribute should have.
 */
template<typename T>
inline void save_attribute(H5::H5File file, const std::string& name, const T& value) {
  save_scalar_attribute(file, name, std::to_string(value));
}

/**
 * @brief Overload for strings. See above.
 */
template<>
inline void save_attribute(H5::H5File file, const std::string& name, const std::string& value) {
  save_scalar_attribute(file, name, value);
}

} // namespace HDF5
} // namespace Serenity
#endif /* BASICS_IO_HDF5_H_ */
