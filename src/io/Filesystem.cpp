/**
 * @file   Filesystem.cpp
 *
 * @date   October 31, 2019
 * @author Jan Unsleber
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#include <iostream>
#include <string>
#include <sys/stat.h>
#include <errno.h>
#if defined(_WIN32)
#include <direct.h>
#endif

namespace Serenity {

bool directoryExists(const std::string& path)
{
#if defined(_WIN32)
  struct _stat info;
  if (_stat(path.c_str(), &info) != 0){
    return false;
  }
  return (info.st_mode & _S_IFDIR) != 0;
#else
  struct stat info;
  if (stat(path.c_str(), &info) != 0){
    return false;
  }
  return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool makePath(const std::string& path){
#if defined(_WIN32)
  int ret = _mkdir(path.c_str());
#else
  mode_t mode = 0755;
  int ret = mkdir(path.c_str(), mode);
#endif
  if (ret == 0)
    return true;

  switch (errno){
    case ENOENT:
      // parent didn't exist, try to create it
      {
      size_t pos = path.find_last_of('/');
      if (pos == std::string::npos)
#if defined(_WIN32)
        pos = path.find_last_of('\\');
      if (pos == std::string::npos)
#endif
        return false;
      if (!makePath( path.substr(0, pos) ))
        return false;
      }
    // now, try to create again
#if defined(_WIN32)
      return 0 == _mkdir(path.c_str());
#else
      return 0 == mkdir(path.c_str(), mode);
#endif
    case EEXIST:
      // done!
      return directoryExists(path);
    default:
      return false;
  } /* end of switch */
}

} /* namespace Serenity */
