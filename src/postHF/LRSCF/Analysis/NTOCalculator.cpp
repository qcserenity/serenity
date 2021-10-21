/**
 * @file NTOCalculator.cpp
 *
 * @date Aug 26, 2019
 * @author Johannes Toelle
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

/* Include Class Header*/
#include "postHF/LRSCF/Analysis/NTOCalculator.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "io/HDF5.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/LRSCFOptions.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <stdlib.h>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
NTOCalculator<SCFMode>::NTOCalculator(std::vector<std::shared_ptr<SystemController>> activeSystems,
                                      std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                      const double plottingThreshold)
  : _activeSystems(activeSystems),
    _environmentSystems(environmentSystems),
    _plottingThreshold(plottingThreshold),
    _hasBeenCalculated(false),
    _state(999999) {
  // Each active subsystem has his own set of occ virt NTO matrices
  _occNTOs.resize(activeSystems.size());
  _virtNTOs.resize(activeSystems.size());
  // read information printed by the LRSCF
  LRSCFTaskSettings lrscfSettings;
  if (activeSystems.size() > 1) {
    lrscfSettings.loadType = Options::LRSCF_TYPE::COUPLED;
  }
  else if (environmentSystems.size() > 0) {
    lrscfSettings.loadType = Options::LRSCF_TYPE::UNCOUPLED;
  }
  else {
    lrscfSettings.loadType = Options::LRSCF_TYPE::ISOLATED;
  }
  // Read excitation vectors and data
  printBigCaption("NTO Calculation");
  for (auto sys : activeSystems) {
    auto temp = std::make_shared<LRSCFController<SCFMode>>(sys, lrscfSettings);
    _lrscf.push_back(temp);
  }

  if (lrscfSettings.loadType != Options::LRSCF_TYPE::COUPLED) {
    std::shared_ptr<std::vector<Eigen::MatrixXd>> excVec = _lrscf[0]->getExcitationVectors(lrscfSettings.loadType);
    _eigenvalues = *(_lrscf[0]->getExcitationEnergies(lrscfSettings.loadType));
    _XPY = (*excVec)[0] + (*excVec)[1];
  }
  else {
    // Determine Vector dimensions
    unsigned int dimension = 0;
    _eigenvalues = *(_lrscf[0]->getExcitationEnergies(lrscfSettings.loadType));
    unsigned int nExc = _eigenvalues.size();
    for (auto lrscf : _lrscf) {
      auto nOcc = lrscf->getNOccupied();
      auto nVirt = lrscf->getNVirtual();
      for_spin(nOcc, nVirt) {
        dimension += nOcc_spin * nVirt_spin;
      };
    }
    _XPY.conservativeResize(dimension, nExc);
    _XPY.setZero();
    // Fill excitation vector with coupled data
    unsigned int rowIndex = 0;
    for (auto lrscf : _lrscf) {
      auto excVec = lrscf->getExcitationVectors(lrscfSettings.loadType);
      if ((*excVec)[0].cols() != nExc) {
        throw SerenityError("The coupled excitation vectors have different number of excitations!");
      }
      _XPY.block(rowIndex, 0, (*excVec)[0].rows(), nExc) = (*excVec)[0] + (*excVec)[1];
      rowIndex += (*excVec)[0].rows();
    }
  }
}

template<Options::SCF_MODES SCFMode>
void NTOCalculator<SCFMode>::calcNTOs(int iState) {
  // Get number of occupied and virtual orbitals
  SpinPolarizedData<SCFMode, unsigned int> nOccSpinTot;
  SpinPolarizedData<SCFMode, unsigned int> nVirtSpinTot;
  unsigned nOcctot = 0;
  unsigned nVirttot = 0;
  for_spin(nOccSpinTot, nVirtSpinTot) {
    nOccSpinTot_spin = 0;
    nVirtSpinTot_spin = 0;
  };
  // Total number of occupied and virtual orbitals
  for (auto lrscf : _lrscf) {
    auto nOcc = lrscf->getNOccupied();
    auto nVirt = lrscf->getNVirtual();
    for_spin(nOccSpinTot, nVirtSpinTot, nOcc, nVirt) {
      nOccSpinTot_spin += nOcc_spin;
      nVirtSpinTot_spin += nVirt_spin;
      nOcctot += nOcc_spin;
      nVirttot += nVirt_spin;
    };
  }
  Eigen::MatrixXd g_total = Eigen::MatrixXd::Zero(nOcctot, nVirttot);
  // Vector to matrix transformation
  unsigned iStart_subsystem = 0;
  unsigned nOcc_start = 0;
  unsigned nVirt_start = 0;
  for (auto lrscf : _lrscf) {
    auto nOcc = lrscf->getNOccupied();
    auto nVirt = lrscf->getNVirtual();
    unsigned int iStart_spin = 0;
    for_spin(nOcc, nVirt) {
      for (unsigned int ia = 0; ia < nOcc_spin * nVirt_spin; ++ia) {
        unsigned int i = floor(ia / nVirt_spin);
        unsigned int a = ia - i * nVirt_spin;
        g_total(i + nOcc_start, a + nVirt_start) = _XPY(iStart_spin + ia + iStart_subsystem, iState);
      }
      iStart_spin += nOcc_spin * nVirt_spin;
      nOcc_start += nOcc_spin;
      nVirt_start += nVirt_spin;
    };
    iStart_subsystem += iStart_spin;
  }
  // Calculate g g^T and g^T g
  Eigen::MatrixXd ggT = g_total * g_total.transpose();
  Eigen::MatrixXd gTg = g_total.transpose() * g_total;
  // Determine eigenpairs
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> u;
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> v;
  nOcc_start = 0;
  // Particle
  for_spin(nOccSpinTot, _occEigenvalues, u) {
    es.compute(ggT.block(nOcc_start, nOcc_start, nOccSpinTot_spin, nOccSpinTot_spin));
    _occEigenvalues_spin = es.eigenvalues();
    u_spin = es.eigenvectors();
    nOcc_start += nOccSpinTot_spin;
  };
  // Electron
  nVirt_start = 0;
  for_spin(nVirtSpinTot, _virtEigenvalues, v) {
    es.compute(gTg.block(nVirt_start, nVirt_start, nVirtSpinTot_spin, nVirtSpinTot_spin));
    _virtEigenvalues_spin = es.eigenvalues();
    v_spin = es.eigenvectors();
    nVirt_start += nVirtSpinTot_spin;
  };
  // Calculate NTOs
  SpinPolarizedData<SCFMode, unsigned int> nOcc_counter;
  SpinPolarizedData<SCFMode, unsigned int> nVirt_counter;
  for_spin(nOcc_counter, nVirt_counter) {
    nOcc_counter_spin = 0;
    nVirt_counter_spin = 0;
  };
  unsigned ilrscf = 0;
  for (auto lrscf : _lrscf) {
    // Hole
    CoefficientMatrix<SCFMode> coeff = lrscf->getCoefficients();
    SpinPolarizedData<SCFMode, Eigen::MatrixXd> occNTOs;
    auto nOcc = lrscf->getNOccupied();
    for_spin(occNTOs, coeff, nOcc, u, nOcc_counter) {
      occNTOs_spin.resize(coeff_spin.rows(), u_spin.cols());
      occNTOs_spin.setZero();
      occNTOs_spin = coeff_spin.block(0, 0, coeff_spin.rows(), nOcc_spin) *
                     u_spin.block(nOcc_counter_spin, 0, nOcc_spin, u_spin.cols());
      nOcc_counter_spin += nOcc_spin;
    };
    _occNTOs[ilrscf] = occNTOs;
    // Particle
    SpinPolarizedData<SCFMode, Eigen::MatrixXd> virtNTOs;
    auto nVirt = lrscf->getNVirtual();
    for_spin(virtNTOs, coeff, nOcc, nVirt, v, nVirt_counter) {
      virtNTOs_spin.resize(coeff_spin.rows(), v_spin.cols());
      virtNTOs_spin.setZero();
      virtNTOs_spin = coeff_spin.block(0, nOcc_spin, coeff_spin.rows(), nVirt_spin) *
                      v_spin.block(nVirt_counter_spin, 0, nVirt_spin, v_spin.cols());
      nVirt_counter_spin += nVirt_spin;
    };
    _virtNTOs[ilrscf] = virtNTOs;
    ilrscf++;
  }
  // print the NTOs
  int spinCounter = 0;
  for_spin(_occEigenvalues, _virtEigenvalues, u, v) {
    std::string type = "";
    if (SCFMode == Options::SCF_MODES::RESTRICTED) {
      type = "RESTRICTED";
    }
    else if (spinCounter == 0 && SCFMode == Options::SCF_MODES::UNRESTRICTED) {
      type = "ALPHA";
    }
    else if (spinCounter == 0 && SCFMode == Options::SCF_MODES::UNRESTRICTED) {
      type = "BETA";
    }
    printNTOInfo(type, iState, _occEigenvalues_spin, _virtEigenvalues_spin, u_spin, v_spin);
    spinCounter++;
  };
  _hasBeenCalculated = true;
  _state = iState;
}

template<Options::SCF_MODES SCFMode>
void NTOCalculator<SCFMode>::printNTOInfo(std::string spin, const unsigned int iState, Eigen::VectorXd& oEigenvalues,
                                          Eigen::VectorXd& vEigenvalues, Eigen::MatrixXd& u, Eigen::MatrixXd& v) {
  for (auto sys : _activeSystems) {
    printSmallCaption((std::string) "\n NTOs for state " + (iState + 1) + " " + spin);
    // Make directory
    std::string dirName = sys->getSystemPath() + "/NTOS" + "/NTO" + std::to_string(iState + 1) + "/";
    std::string command = "mkdir -p " + dirName;
    auto stat = system(command.c_str());
    (void)stat;
    std::cout << " Plotting NTOs to " << dirName << "... \n" << std::endl;
    // create or update readme file
    FILE* readme;
    readme = fopen((dirName + "README").c_str(), "a");
    // Print occupied eigenvalues to readme
    fprintf(readme, "\n %s : \n", spin.c_str());
    fprintf(readme, "%s \n", " OCCUPIED EIGENVALUES");
    std::cout << " Occupied eigenvalues : " << std::endl;
    for (unsigned int i = 0; i < oEigenvalues.rows(); ++i) {
      if (oEigenvalues(i) > _plottingThreshold) {
        printf("   NTO %3i %3.2f \n ", i + 1, oEigenvalues(i));
        fprintf(readme, "NTO %d %f \n", i + 1, oEigenvalues(i));
      }
    }
    // Print occupied dominant contributions to readme
    fprintf(readme, "\n %s ", "DOMINANT CONTRIBUTIONS");
    print((std::string) "\n Occupied dominant contributions : ");
    for (unsigned int i = 0; i < u.cols(); ++i) {
      if (oEigenvalues(i) > _plottingThreshold) {
        printf("   NTO %3.i :", i + 1);
        fprintf(readme, "\n NTO %d : ", i + 1);
        for (unsigned int j = 0; j < u.rows(); ++j) {
          if (fabs(u(j, i)) > _plottingThreshold) {
            printf("%3.i (%5.3f) ", j + 1, u(j, i));
            fprintf(readme, "%d (%5.3f)  ", j + 1, u(j, i));
          }
        }
        printf("\n");
      }
    }
    // Print virtual eigenvalues to readme
    fprintf(readme, "\n \n %s \n", "VIRTUAL EIGENVALUES");
    print((std::string) "\n Virtual eigenvalues : ");
    for (unsigned int a = 0; a < vEigenvalues.rows(); ++a) {
      if (vEigenvalues((int)(vEigenvalues.rows() - a - 1)) > _plottingThreshold) {
        printf("   NTO %3i %3.2f \n ", (int)(oEigenvalues.rows() + a + 1), vEigenvalues((int)(vEigenvalues.rows() - a - 1)));
        fprintf(readme, "NTO %d %f \n", (int)(oEigenvalues.rows() + a + 1), vEigenvalues((int)(vEigenvalues.rows() - a - 1)));
      }
    }
    // Print virtual dominant contributions to readme
    fprintf(readme, "\n %s ", "DOMINANT CONTRIBUTIONS");
    print((std::string) "\n Virtual dominant contributions : ");
    for (unsigned int a = 0; a < vEigenvalues.rows(); ++a) {
      if (vEigenvalues((int)(vEigenvalues.rows() - a - 1)) > _plottingThreshold) {
        printf("   NTO %3.i :", (int)(oEigenvalues.rows() + a + 1));
        fprintf(readme, "\n NTO %d : ", (int)(oEigenvalues.rows() + a + 1));
        for (unsigned int j = 0; j < v.rows(); ++j) {
          if (fabs(v(j, v.cols() - a - 1)) > _plottingThreshold) {
            printf("%3.i (%5.3f) ", (int)(u.cols() + j + 1), v(j, (int)(v.cols() - a - 1)));
            fprintf(readme, "%d (%5.3f)  ", (int)(u.cols() + j + 1), v(j, (int)(v.cols() - a - 1)));
          }
        }
        printf("\n");
      }
    }
    fclose(readme);
  }
}

template<Options::SCF_MODES SCFMode>
unsigned int NTOCalculator<SCFMode>::getNumberOfStates() {
  return _eigenvalues.size();
}

template<Options::SCF_MODES SCFMode>
std::string NTOCalculator<SCFMode>::getDir(unsigned int i, unsigned int iSubsystem) {
  std::string dirName = _activeSystems[iSubsystem]->getSystemPath() + "/NTOS" + "/NTO" + std::to_string(i + 1) + "/";
  return dirName;
}

template class NTOCalculator<Options::SCF_MODES::RESTRICTED>;
template class NTOCalculator<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */