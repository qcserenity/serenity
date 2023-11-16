/**
 * @file   Filesystem.h
 *
 * @date   October 31, 2019
 * @author Jan Unsleber
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

#ifndef IO_FILESYSTEM_H_
#define IO_FILESYSTEM_H_

/* Include Std and External Headers */
#include <string>

namespace Serenity {

/**
 * @brief Removes any Serenity-system files at the given path and name.
 * @param path        The path to the system directory.
 * @param systemName  The file base name.
 */
void removeSystemFiles(std::string path, std::string systemName, std::string basisLabel = "");
/**
 * @brief Checks if a directory or path exists.
 * @param path The path to check for.
 * @return true  If the path exists.
 * @return false If the path does not exist.
 */
bool directoryExists(const std::string& path);
/**
 * @brief Generates a path (set of nested directories) if it is not available.
 * @param path The path to generate.
 * @return true  If the path was created, or exists.
 * @return false If the path could not be created.
 */
bool makePath(const std::string& path);

} /* namespace Serenity */
#endif /* IO_FILESYSTEM_H_ */
