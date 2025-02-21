/**
 * @file NTOCalculator.cpp
 *
 * @date Aug 26, 2019
 * @author Johannes Toelle, Anton Rikus
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
#include "basis/BasisFunctionMapper.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/wrappers/Libint.h"
#include "io/HDF5.h"
#include "math/linearAlgebra/MatrixFunctions.h"
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
                                      const double ntoThreshold)
  : _activeSystems(activeSystems), _environmentSystems(environmentSystems), _ntoThreshold(ntoThreshold), _state(999999) {
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
    _lrscf.push_back(std::make_shared<LRSCFController<SCFMode>>(sys, lrscfSettings));
    _transitionDensity.push_back(sys->getBasisController());
    _holeDensity.push_back(sys->getBasisController());
    _holeDensityCorrection.push_back(sys->getBasisController());
    _particleDensity.push_back(sys->getBasisController());
  }

  if (lrscfSettings.loadType != Options::LRSCF_TYPE::COUPLED) {
    // in the isolated or uncoupled case, there is only one LRSCFController in _lrscf
    std::shared_ptr<std::vector<Eigen::MatrixXd>> excVec = _lrscf[0]->getExcitationVectors(lrscfSettings.loadType);
    _eigenvalues = *(_lrscf[0]->getExcitationEnergies(lrscfSettings.loadType));
    _XPY = (*excVec)[0];
    _XMY = (*excVec)[1];
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
    _XPY.resize(dimension, nExc);
    _XMY.resize(dimension, nExc);
    // Fill excitation vector with coupled data
    unsigned int rowIndex = 0;
    for (auto lrscf : _lrscf) {
      auto excVec = lrscf->getExcitationVectors(lrscfSettings.loadType);
      if ((*excVec)[0].cols() != nExc) {
        throw SerenityError("The coupled excitation vectors have different number of excitations!");
      }
      _XPY.block(rowIndex, 0, (*excVec)[0].rows(), nExc) = (*excVec)[0];
      _XMY.block(rowIndex, 0, (*excVec)[1].rows(), nExc) = (*excVec)[1];
      rowIndex += (*excVec)[0].rows();
    }
  }
}

template<Options::SCF_MODES SCFMode>
void NTOCalculator<SCFMode>::calcNTOs(unsigned int iState) {
  // Get number of occupied and virtual orbitals
  SpinPolarizedData<SCFMode, unsigned int> nOccSpinTot(0);
  SpinPolarizedData<SCFMode, unsigned int> nVirtSpinTot(0);
  // Total number of occupied and virtual orbitals
  for (auto lrscf : _lrscf) {
    auto nOcc = lrscf->getNOccupied();
    auto nVirt = lrscf->getNVirtual();
    for_spin(nOccSpinTot, nVirtSpinTot, nOcc, nVirt) {
      nOccSpinTot_spin += nOcc_spin;
      nVirtSpinTot_spin += nVirt_spin;
    };
  }
  Eigen::MatrixXd g_total = Eigen::MatrixXd::Zero(nOccSpinTot.total(), nVirtSpinTot.total());
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
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(g_total, Eigen::ComputeFullU | Eigen::ComputeFullV);
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
    else if (spinCounter == 1 && SCFMode == Options::SCF_MODES::UNRESTRICTED) {
      type = "BETA";
    }
    printNTOInfo(type, iState, _occEigenvalues_spin, _virtEigenvalues_spin, u_spin, v_spin);
    spinCounter++;
  };
}

template<Options::SCF_MODES SCFMode>
void NTOCalculator<SCFMode>::calcTransitionDensity(unsigned int iState) {
  unsigned iaStart = 0;
  for (unsigned iLRSCF = 0; iLRSCF < _lrscf.size(); iLRSCF++) {
    auto lrscf = _lrscf[iLRSCF];
    const SpinPolarizedData<SCFMode, unsigned int>& nOcc = lrscf->getNOccupied();
    const SpinPolarizedData<SCFMode, unsigned int>& nVirt = lrscf->getNVirtual();
    const CoefficientMatrix<SCFMode>& coeff = lrscf->getCoefficients();
    auto& tmp = _transitionDensity[iLRSCF];
    for_spin(nOcc, nVirt, coeff, tmp) {
      const Eigen::VectorXd& xpyBlock = _XPY.block(iaStart, iState, nOcc_spin * nVirt_spin, 1);
      tmp_spin = coeff_spin.leftCols(nOcc_spin) *
                 Eigen::Map<const Eigen::MatrixXd>(xpyBlock.data(), nVirt_spin, nOcc_spin).transpose() *
                 coeff_spin.middleCols(nOcc_spin, nVirt_spin).transpose();
      iaStart += nOcc_spin * nVirt_spin;
    };
  }
}

template<Options::SCF_MODES SCFMode>
void NTOCalculator<SCFMode>::calcParticleAndHoleDensity(unsigned int iState) {
  unsigned iaStart = 0;
  for (unsigned iLRSCF = 0; iLRSCF < _lrscf.size(); iLRSCF++) {
    auto lrscf = _lrscf[iLRSCF];
    const SpinPolarizedData<SCFMode, unsigned int>& nOcc = lrscf->getNOccupied();
    const SpinPolarizedData<SCFMode, unsigned int>& nVirt = lrscf->getNVirtual();
    const CoefficientMatrix<SCFMode>& coeff = lrscf->getCoefficients();
    auto& hole = _holeDensity[iLRSCF];
    auto& particle = _particleDensity[iLRSCF];
    for_spin(nOcc, nVirt, coeff, hole, particle) {
      const Eigen::VectorXd& xpyBlock = _XPY.block(iaStart, iState, nOcc_spin * nVirt_spin, 1);
      const Eigen::VectorXd& xmyBlock = _XMY.block(iaStart, iState, nOcc_spin * nVirt_spin, 1);
      hole_spin = coeff_spin.leftCols(nOcc_spin) * (-0.5) *
                  (Eigen::Map<const Eigen::MatrixXd>(xpyBlock.data(), nVirt_spin, nOcc_spin).transpose() *
                       Eigen::Map<const Eigen::MatrixXd>(xpyBlock.data(), nVirt_spin, nOcc_spin) +
                   Eigen::Map<const Eigen::MatrixXd>(xmyBlock.data(), nVirt_spin, nOcc_spin).transpose() *
                       Eigen::Map<const Eigen::MatrixXd>(xmyBlock.data(), nVirt_spin, nOcc_spin)) *
                  coeff_spin.leftCols(nOcc_spin).transpose();
      particle_spin = coeff_spin.middleCols(nOcc_spin, nVirt_spin) * 0.5 *
                      (Eigen::Map<const Eigen::MatrixXd>(xpyBlock.data(), nVirt_spin, nOcc_spin) *
                           Eigen::Map<const Eigen::MatrixXd>(xpyBlock.data(), nVirt_spin, nOcc_spin).transpose() +
                       Eigen::Map<const Eigen::MatrixXd>(xmyBlock.data(), nVirt_spin, nOcc_spin) *
                           Eigen::Map<const Eigen::MatrixXd>(xmyBlock.data(), nVirt_spin, nOcc_spin).transpose()) *
                      coeff_spin.middleCols(nOcc_spin, nVirt_spin).transpose();
    };
    unsigned iaStartJ = 0;
    for (unsigned jLRSCF = 0; jLRSCF < _lrscf.size(); jLRSCF++) {
      auto jlrscf = _lrscf[jLRSCF];
      const SpinPolarizedData<SCFMode, unsigned int>& nOccJ = jlrscf->getNOccupied();
      // remember that _lrscf only contains more than one LRSCFController in the coupled case
      if (jLRSCF != iLRSCF) {
        auto& libint = Libint::getInstance();
        Eigen::MatrixXd overlapIJ =
            libint.compute1eInts(LIBINT_OPERATOR::overlap, lrscf->getBasisController(), jlrscf->getBasisController());
        BasisFunctionMapper mapper(lrscf->getBasisController());
        auto sparseProj = mapper.getSparseProjection(jlrscf->getBasisController());
        const CoefficientMatrix<SCFMode>& coeffJ = jlrscf->getCoefficients();
        auto& holeCorr = _holeDensityCorrection[iLRSCF];
        std::shared_ptr<Eigen::SparseMatrix<double>> proj = mapper.getSparseProjection(jlrscf->getBasisController());
        for_spin(nOcc, nVirt, coeff, holeCorr, coeffJ, nOccJ) {
          Eigen::MatrixXd overlapIJMO = coeff_spin.middleCols(nOcc_spin, nVirt_spin).transpose() *
                                        overlapIJ.transpose() * coeffJ_spin.middleCols(nOccJ_spin, nVirt_spin);
          const Eigen::VectorXd& xpyBlockI = _XPY.block(iaStart, iState, nOcc_spin * nVirt_spin, 1);
          const Eigen::VectorXd& xmyBlockI = _XMY.block(iaStart, iState, nOcc_spin * nVirt_spin, 1);
          const Eigen::VectorXd& xpyBlockJ = _XPY.block(iaStartJ, iState, nOccJ_spin * nVirt_spin, 1);
          const Eigen::VectorXd& xmyBlockJ = _XMY.block(iaStartJ, iState, nOccJ_spin * nVirt_spin, 1);
          holeCorr_spin += coeff_spin.leftCols(nOcc_spin) * (-0.5) *
                           (Eigen::Map<const Eigen::MatrixXd>(xpyBlockI.data(), nVirt_spin, nOcc_spin).transpose() *
                                overlapIJMO * Eigen::Map<const Eigen::MatrixXd>(xpyBlockJ.data(), nVirt_spin, nOccJ_spin) +
                            Eigen::Map<const Eigen::MatrixXd>(xmyBlockI.data(), nVirt_spin, nOcc_spin).transpose() *
                                overlapIJMO * Eigen::Map<const Eigen::MatrixXd>(xmyBlockJ.data(), nVirt_spin, nOccJ_spin)) *
                           (*proj * coeffJ_spin.leftCols(nOccJ_spin)).transpose();
          iaStartJ += nOccJ_spin * nVirt_spin;
        };
      }
      else {
        for_spin(nOccJ, nVirt) {
          iaStartJ += nOccJ_spin * nVirt_spin;
        };
      }
    } // for jLRSCF
    for_spin(nOcc, nVirt) {
      iaStart += nOcc_spin * nVirt_spin;
    };
  } // for iLRSCF
}

template<Options::SCF_MODES SCFMode>
void NTOCalculator<SCFMode>::printNTOInfo(std::string spin, const unsigned int iState, Eigen::VectorXd& oEigenvalues,
                                          Eigen::VectorXd& vEigenvalues, Eigen::MatrixXd& u, Eigen::MatrixXd& v) {
  for (auto sys : _activeSystems) {
    printSmallCaption((std::string) "\n NTOs for state " + (iState + 1) + " " + spin);
    // Make directory
    std::string dirName = sys->getSystemPath() + "NTOS" + "/NTO" + std::to_string(iState + 1) + "/";
    std::string command = "mkdir -p " + dirName;
    // avoid unused result compile error
    auto stat = system(command.c_str());
    (void)stat;
    std::cout << " Plotting NTOs to " << dirName << "... \n" << std::endl;
    // create or update readme file
    std::remove((dirName + "README").c_str());
    FILE* readme;
    readme = fopen((dirName + "README").c_str(), "a");
    // Print occupied eigenvalues to readme
    fprintf(readme, "\n %s : \n", spin.c_str());
    fprintf(readme, "%s \n", " OCCUPIED EIGENVALUES");
    std::cout << " Occupied eigenvalues : " << std::endl;
    for (unsigned int i = 0; i < oEigenvalues.rows(); ++i) {
      if (oEigenvalues(i) > _ntoThreshold) {
        printf("   NTO %3i %3.2f \n ", i + 1, oEigenvalues(i));
        fprintf(readme, "NTO %d %f \n", i + 1, oEigenvalues(i));
      }
    }
    // Print occupied dominant contributions to readme
    fprintf(readme, "\n %s ", "DOMINANT CONTRIBUTIONS");
    print((std::string) "\n Occupied dominant contributions : ");
    for (unsigned int i = 0; i < u.cols(); ++i) {
      if (oEigenvalues(i) > _ntoThreshold) {
        printf("   NTO %3.i :", i + 1);
        fprintf(readme, "\n NTO %d : ", i + 1);
        for (unsigned int j = 0; j < u.rows(); ++j) {
          if (fabs(u(j, i)) > _ntoThreshold) {
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
      if (vEigenvalues((int)(vEigenvalues.rows() - a - 1)) > _ntoThreshold) {
        printf("   NTO %3i %3.2f \n ", (int)(oEigenvalues.rows() + a + 1), vEigenvalues((int)(vEigenvalues.rows() - a - 1)));
        fprintf(readme, "NTO %d %f \n", (int)(oEigenvalues.rows() + a + 1), vEigenvalues((int)(vEigenvalues.rows() - a - 1)));
      }
    }
    // Print virtual dominant contributions to readme
    fprintf(readme, "\n %s ", "DOMINANT CONTRIBUTIONS");
    print((std::string) "\n Virtual dominant contributions : ");
    for (unsigned int a = 0; a < vEigenvalues.rows(); ++a) {
      if (vEigenvalues((int)(vEigenvalues.rows() - a - 1)) > _ntoThreshold) {
        printf("   NTO %3.i :", (int)(oEigenvalues.rows() + a + 1));
        fprintf(readme, "\n NTO %d : ", (int)(oEigenvalues.rows() + a + 1));
        for (unsigned int j = 0; j < v.rows(); ++j) {
          if (fabs(v(j, v.cols() - a - 1)) > _ntoThreshold) {
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