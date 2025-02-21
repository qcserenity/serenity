/**
 * @file   Filesystem.cpp
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
/* Include Std and External Headers */
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <iostream>
#include <string>
#if defined(_WIN32)
#include <direct.h>
#endif

namespace Serenity {

void removeSystemFiles(std::string path, std::string systemName, std::string basisLabel) {
  /* TODO
   * Replace with std::filesystem when c++17 is available
   */
  std::cout << "\n" << path << std::endl;
  std::string baseName = path + systemName;
  std::remove((baseName + ".settings").c_str());
  std::remove((baseName + ".xyz").c_str());
  std::remove((baseName + ".orbs.res.h5").c_str());
  std::remove((baseName + ".orbs.unres.h5").c_str());
  std::remove((baseName + ".FockMatrix.res.h5").c_str());
  std::remove((baseName + ".FockMatrix.unres.h5").c_str());
  std::remove((baseName + ".energies.res").c_str());
  std::remove((baseName + ".energies.unres").c_str());
  std::remove((baseName + ".dmat.res.h5").c_str());
  std::remove((baseName + ".dmat.unres.h5").c_str());
  std::remove((baseName + ".basis.h5").c_str());
  std::remove((baseName + ".hess.h5").c_str());
  std::remove((baseName + ".elecPotInts.h5").c_str());
  std::remove((baseName + "_lrscf.settings").c_str());
  std::remove((baseName + "_lrscf.iso.res.h5").c_str());
  std::remove((baseName + "_lrscf.iso.unres.h5").c_str());
  std::remove((baseName + "_lrscf.fdeu.res.h5").c_str());
  std::remove((baseName + "_lrscf.fdeu.unres.h5").c_str());
  std::remove((baseName + "_lrscf.fdec.res.h5").c_str());
  std::remove((baseName + "_lrscf.fdec.unres.h5").c_str());
  std::remove((baseName + "_lrscf_resp.iso.res.h5").c_str());
  std::remove((baseName + "_lrscf_resp.iso.unres.h5").c_str());
  std::remove((baseName + "_lrscf_resp.fdeu.res.h5").c_str());
  std::remove((baseName + "_lrscf_resp.fdeu.unres.h5").c_str());
  std::remove((baseName + "_lrscf_resp.fdec.res.h5").c_str());
  std::remove((baseName + "_lrscf_resp.fdec.unres.h5").c_str());
  std::remove((baseName + ".exspectrum.txt").c_str());
  std::remove((baseName + ".transitioncharges.txt").c_str());
  std::remove((baseName + ".cd.h-ERF.h5").c_str());
  std::remove((baseName + ".cd.h-ERF_EXPANDED.h5").c_str());
  std::remove((baseName + ".cd.h.h5").c_str());
  std::remove((baseName + ".cd.h_EXPANDED.h5").c_str());
  std::remove((baseName + ".cd.c.h5").c_str());
  std::remove((baseName + ".cd.c_EXPANDED.h5").c_str());
  std::remove((baseName + ".cd.o.h5").c_str());
  std::remove((baseName + ".cd.o_EXPANDED.h5").c_str());
  std::remove((baseName + ".cd.AO.h5").c_str());
  std::remove((baseName + "_cc2_dens.res.h5").c_str());
  std::remove((baseName + "_cc2_dens.unres.h5").c_str());
  std::remove((baseName + "_cc2_dens.fdec.res.h5").c_str());
  std::remove((baseName + "_cc2_dens.fdec.unres.h5").c_str());
  std::remove((path + "ACD-" + basisLabel).c_str());
  std::remove((path + "ACD-" + basisLabel + "-ERF").c_str());
  std::remove((path + "ACCD-" + basisLabel).c_str());
  std::remove((path + "ACCD-" + basisLabel + "-ERF").c_str());
  std::remove((path).c_str());
  std::remove("WARNING");
}

bool directoryExists(const std::string& path) {
#if defined(_WIN32)
  struct _stat info;
  if (_stat(path.c_str(), &info) != 0) {
    return false;
  }
  return (info.st_mode & _S_IFDIR) != 0;
#else
  struct stat info;
  if (stat(path.c_str(), &info) != 0) {
    return false;
  }
  return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool makePath(const std::string& inpt) {
  std::string path = inpt;
#if defined(_WIN32)
  if (path[path.size() - 1] == '\\')
    path = path.substr(0, path.size() - 1);
#else
  if (path[path.size() - 1] == '/')
    path = path.substr(0, path.size() - 1);
#endif
#if defined(_WIN32)
  int ret = _mkdir(path.c_str());
#else
  mode_t mode = 0755;
  int ret = mkdir(path.c_str(), mode);
#endif
  if (ret == 0)
    return true;

  switch (errno) {
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
        if (!makePath(path.substr(0, pos)))
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
