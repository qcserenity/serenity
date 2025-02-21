/**
 * @file LRSCFRestart.cpp
 *
 * @date Jan 6, 2020
 * @author Niklas Niemeyer
 *
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
#include "postHF/LRSCF/Tools/LRSCFRestart.h"
/* Include Serenity Internal Headers */
#include "io/Filesystem.h"
#include "io/HDF5.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

/* Include Std and External Headers */

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LRSCFRestart<SCFMode>::LRSCFRestart(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                    const LRSCFTaskSettings& settings, Options::LRSCF_TYPE type)
  : _lrscf(lrscf), _settings(settings), _type(type) {
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<std::vector<Eigen::MatrixXd>> LRSCFRestart<SCFMode>::fetchEigenpairs(Eigen::VectorXd& eigenvalues) {
  printBigCaption("Restart from Disk");
  try {
    std::vector<Eigen::MatrixXd> eigenvectors(2, Eigen::MatrixXd(0, 0));
    Eigen::MatrixXd tmp;

    unsigned iStart = 0;
    for (auto& lrscf : _lrscf) {
      auto sysSettings = lrscf->getSysSettings();

      // prepare hdf5 file
      std::string fName = ((sysSettings.load == "") ? sysSettings.path : sysSettings.load) + sysSettings.name + "_lrscf.";
      if (_type == Options::LRSCF_TYPE::ISOLATED) {
        fName += "iso.";
      }
      else if (_type == Options::LRSCF_TYPE::UNCOUPLED) {
        fName += "fdeu.";
      }
      else {
        fName += "fdec.";
      }
      fName += (SCFMode == RESTRICTED) ? "res.h5" : "unres.h5";

      // get lrscf dimensions
      auto no = lrscf->getNOccupied();
      auto nv = lrscf->getNVirtual();
      unsigned nDimI = 0;
      for_spin(no, nv) {
        nDimI += no_spin * nv_spin;
      };

      // try to find a converged solution first and then non-converged
      if (std::ifstream(fName).good()) {
        printf("  Restarting from converged solution found in:\n\n   $  %-20s\n\n", fName.c_str());
      }
      else if (std::ifstream(fName + ".tmp").good()) {
        fName += ".tmp";
        printf("  Restarting from non-converged solution found in:\n\n   $  %-20s\n\n", fName.c_str());
      }
      else {
        throw SerenityError("Could not find any old files.");
      }

      // initialize file
      HDF5::H5File file(fName.c_str(), H5F_ACC_RDONLY);

      // fill up eigenvectors with loaded data
      for (unsigned iSet = 0; iSet < 2; ++iSet) {
        HDF5::dataset_exists(file, (iSet == 0) ? "RIGHT" : "LEFT");
        HDF5::load(file, (iSet == 0) ? "RIGHT" : "LEFT", tmp);
        if (tmp.rows() != nDimI) {
          throw SerenityError("The dimension of your loaded eigenpairs (" + std::to_string(tmp.rows()) +
                              ") does not match with this response problem (" + std::to_string(nDimI) + ").");
        }
        eigenvectors[iSet].conservativeResize(iStart + nDimI, tmp.cols());
        eigenvectors[iSet].middleRows(iStart, nDimI) = tmp;
      }

      HDF5::dataset_exists(file, "EIGENVALUES");
      HDF5::load(file, "EIGENVALUES", eigenvalues);
      iStart += nDimI;
      file.close();
    } /* loop over lrscf controller */

    printf("  Successfully loaded %3i eigenpairs from disk.\n\n", (int)eigenvectors[0].cols());
    return std::make_shared<std::vector<Eigen::MatrixXd>>(eigenvectors);
  }
  catch (...) {
    printf("  Will continue from scratch.\n\n");
    return nullptr;
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFRestart<SCFMode>::storeConvergedSolution(std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors,
                                                   Eigen::VectorXd eigenvalues) {
  unsigned iStart = 0;
  for (auto& lrscf : _lrscf) {
    // remove temporary eigenvector files
    std::string fName = lrscf->getSys()->getSystemPath() + lrscf->getSys()->getSystemName() + "_lrscf.";
    if (_type == Options::LRSCF_TYPE::ISOLATED) {
      fName += "iso.";
    }
    else if (_type == Options::LRSCF_TYPE::UNCOUPLED) {
      fName += "fdeu.";
    }
    else {
      fName += "fdec.";
    }
    fName += (SCFMode == RESTRICTED) ? "res.h5.tmp" : "unres.h5.tmp";
    std::remove(fName.c_str());
    // store converged solution
    if (_type != Options::LRSCF_TYPE::COUPLED) {
      lrscf->setSolution(eigenvectors, std::make_shared<Eigen::VectorXd>(eigenvalues), _type);
    }
    else {
      auto no = lrscf->getNOccupied();
      auto nv = lrscf->getNVirtual();
      auto vec = std::make_shared<std::vector<Eigen::MatrixXd>>(2);
      unsigned nDimI = 0;
      for_spin(no, nv) {
        nDimI += no_spin * nv_spin;
      };
      (*vec)[0] = (*eigenvectors)[0].middleRows(iStart, nDimI);
      (*vec)[1] = (*eigenvectors)[1].middleRows(iStart, nDimI);
      iStart += nDimI;
      lrscf->setSolution(vec, std::make_shared<Eigen::VectorXd>(eigenvalues), _type);
    }
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFRestart<SCFMode>::storeConvergedResponse(std::shared_ptr<std::vector<Eigen::MatrixXd>> solutionvectors,
                                                   Eigen::VectorXd frequencies) {
  unsigned iStart = 0;
  for (auto& lrscf : _lrscf) {
    // not coupled means there is only one LRSCFController in _lrscf
    if (_type != Options::LRSCF_TYPE::COUPLED) {
      std::string fileName = lrscf->getSys()->getSystemPath() + lrscf->getSys()->getSystemName() + "_lrscf_resp.";
      if (_type == Options::LRSCF_TYPE::ISOLATED) {
        fileName += "iso.";
      }
      else if (_type == Options::LRSCF_TYPE::UNCOUPLED) {
        fileName += "fdeu.";
      }
      fileName += (SCFMode == RESTRICTED) ? "res." : "unres.";
      fileName += "h5";

      HDF5::H5File file(fileName, H5F_ACC_TRUNC);
      HDF5::save_scalar_attribute(file, "ID", lrscf->getSys()->getSystemIdentifier());
      HDF5::save(file, "X+Y", (*solutionvectors)[0]);
      HDF5::save(file, "X-Y", (*solutionvectors)[1]);
      HDF5::save(file, "frequencies", frequencies);
      file.close();
    }
    // this is the coupled case
    else {
      auto no = lrscf->getNOccupied();
      auto nv = lrscf->getNVirtual();
      auto vec = std::make_shared<std::vector<Eigen::MatrixXd>>(2);
      unsigned nDimI = 0;
      for_spin(no, nv) {
        nDimI += no_spin * nv_spin;
      };
      (*vec)[0] = (*solutionvectors)[0].middleRows(iStart, nDimI);
      if (_settings.method == Options::LR_METHOD::TDDFT) {
        (*vec)[1] = (*solutionvectors)[1].middleRows(iStart, nDimI);
      }
      iStart += nDimI;
      std::string fileName = lrscf->getSys()->getSystemPath() + lrscf->getSys()->getSystemName() + "_lrscf_resp.fdec.";
      fileName += (SCFMode == RESTRICTED) ? "res." : "unres.";
      fileName += "h5";

      HDF5::H5File file(fileName, H5F_ACC_TRUNC);
      HDF5::save_scalar_attribute(file, "ID", lrscf->getSys()->getSystemIdentifier());
      HDF5::save(file, "X+Y", (*vec)[0]);
      HDF5::save(file, "X-Y", (*vec)[1]);
      HDF5::save(file, "frequencies", frequencies);
      file.close();
    }
  }
}

template<Options::SCF_MODES SCFMode>
std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> LRSCFRestart<SCFMode>::getWriteToDisk() {
  std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> func =
      [&](std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd& eigenvalues) {
        Timings::takeTime("LRSCF -              Disk I/O");
        unsigned iStart = 0;
        for (auto& lrscf : _lrscf) {
          // prepare hdf5 file
          std::string fName = lrscf->getSys()->getSystemPath() + lrscf->getSys()->getSystemName() + "_lrscf.";
          if (_type == Options::LRSCF_TYPE::ISOLATED) {
            fName += "iso.";
          }
          else if (_type == Options::LRSCF_TYPE::UNCOUPLED) {
            fName += "fdeu.";
          }
          else {
            fName += "fdec.";
          }
          fName += (SCFMode == RESTRICTED) ? "res." : "unres.";
          fName += "h5.tmp";
          HDF5::H5File file(fName.c_str(), H5F_ACC_TRUNC);

          // initialize output eigenvectors
          unsigned nSets = eigenvectors.size();
          std::vector<Eigen::MatrixXd> eigenvectorsI(nSets);

          // get lrscf dimensions
          auto no = lrscf->getNOccupied();
          auto nv = lrscf->getNVirtual();
          unsigned nDimI = 0;
          for_spin(no, nv) {
            nDimI += no_spin * nv_spin;
          };
          // save eigenvectors of this system to disk
          for (unsigned iSet = 0; iSet < nSets; ++iSet) {
            eigenvectorsI[iSet] = eigenvectors[iSet].middleRows(iStart, nDimI);
            HDF5::save(file, iSet == 0 ? "RIGHT" : "LEFT", eigenvectorsI[iSet]);
          }
          HDF5::save(file, "EIGENVALUES", eigenvalues);
          iStart += nDimI;
          file.close();
        }
        Timings::timeTaken("LRSCF -              Disk I/O");
      };
  return func;
}

template class LRSCFRestart<Options::SCF_MODES::RESTRICTED>;
template class LRSCFRestart<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
