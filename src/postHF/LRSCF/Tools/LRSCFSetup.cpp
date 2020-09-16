/**
 * @file LRSCFSetup.cpp
 * @date Jan. 10, 2019
 * @author Johannes Tölle
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
#include "postHF/LRSCF/Tools/LRSCFSetup.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "geometry/Point.h"
#include "integrals/wrappers/Libint.h"
#include "settings/LRSCFOptions.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void LRSCFSetup<SCFMode>::printInfo(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                    const LRSCFTaskSettings& settings,
                                    const std::vector<std::shared_ptr<SystemController>>& envSys,
                                    const Options::LRSCF_TYPE type) {
  // Some output (Information of the calculation)
  // Dummy just needed to compile for_spin macro

  SpinPolarizedData<SCFMode, int> dummy;

  printf(" System(s):\n");
  printf("%5s %20s %5s ", "", "Sys", "Stat");
  unsigned int iStart = 0;
  for_spin(dummy) {
    (void)dummy_spin; // no warning
    printf("%8s%1s ", "nOcc_", (iStart == 0) ? "a" : "b");
    iStart += 1;
  };
  iStart = 0;
  for_spin(dummy) {
    (void)dummy_spin; // no warning
    printf("%8s%1s ", "nVirt_", (iStart == 0) ? "a" : "b");
    iStart += 1;
  };
  printf("%8s\n", "nDim");
  printf(" -");
  for_spin(dummy) {
    (void)dummy_spin; // no warning
    printf("-------------------");
  };
  printf("--------------------------------------------\n");
  for (unsigned int I = 0; I < lrscf.size(); ++I) {
    printf("%5i ", I + 1);
    SpinPolarizedData<SCFMode, unsigned int> nOccupied = lrscf[I]->getNOccupied();
    SpinPolarizedData<SCFMode, unsigned int> nVirtual = lrscf[I]->getNVirtual();
    printf("%20s ", lrscf[I]->getSys()->getSystemName().c_str());
    printf("%5s ", "act");
    for_spin(nOccupied) {
      printf("%8i ", nOccupied_spin);
    };
    for_spin(nVirtual) {
      printf("%8i ", nVirtual_spin);
    };
    unsigned int nDimI = 0;
    for_spin(nOccupied, nVirtual) {
      nDimI += nOccupied_spin * nVirtual_spin;
    };
    printf("%8i \n", nDimI);
  }
  for (unsigned int I = 0; I < envSys.size(); ++I) {
    printf("%5i ", I + 1);
    auto nOccupied = envSys[I]->getNOccupiedOrbitals<SCFMode>();
    auto nVirtual = envSys[I]->getNVirtualOrbitals<SCFMode>();
    printf("%20s ", envSys[I]->getSystemPath().c_str());
    printf("%5s ", "env");
    for_spin(nOccupied) {
      printf("%8i ", nOccupied_spin);
    };
    for_spin(nVirtual) {
      printf("%8i ", nVirtual_spin);
    };
    unsigned int nDimI = 0;
    for_spin(nOccupied, nVirtual) {
      nDimI += nOccupied_spin * nVirtual_spin;
    };
    printf("%8i \n", nDimI);
  }

  printf("\n Type        : %22s \n",
         (type == Options::LRSCF_TYPE::ISOLATED) ? "Isolated" : (type == Options::LRSCF_TYPE::UNCOUPLED) ? "FDEu" : "FDEc");
  std::string naddXCFunc_string;
  std::string naddKinFunc_string;
  std::string func;
  auto naddxc = settings.embedding.naddXCFunc;
  auto naddkin = settings.embedding.naddKinFunc;
  auto xc = settings.func;
  Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(naddXCFunc_string, naddxc);
  Options::resolve<CompositeFunctionals::KINFUNCTIONALS>(naddKinFunc_string, naddkin);
  Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(func, xc);
  printf(" NaddKinFunc : %22s \n", naddKinFunc_string.c_str());
  printf(" NaddXCFunc  : %22s \n\n", naddXCFunc_string.c_str());
  if (settings.func != CompositeFunctionals::XCFUNCTIONALS::NONE)
    printf(" Func        : %22s \n\n", func.c_str());
}

template<Options::SCF_MODES SCFMode>
Point LRSCFSetup<SCFMode>::getGaugeOrigin(const LRSCFTaskSettings& settings,
                                          const std::vector<std::shared_ptr<SystemController>>& act,
                                          const std::vector<std::shared_ptr<SystemController>>& env) {
  // Set gauge-origin for dipole integrals
  // Create point object to pass to dipole integral constructor
  Point gaugeOrigin(settings.gaugeOrigin[0] / BOHR_TO_ANGSTROM, settings.gaugeOrigin[1] / BOHR_TO_ANGSTROM,
                    settings.gaugeOrigin[2] / BOHR_TO_ANGSTROM);
  // The default values are set to infinity to detect whether a custom input was given or not
  // Default: no custom input was given, gauge-origin is being set to center of mass
  if (std::find(settings.gaugeOrigin.begin(), settings.gaugeOrigin.end(), std::numeric_limits<double>::infinity()) !=
      settings.gaugeOrigin.end()) {
    printf("\n Set gauge-origin to center of mass.\n");
    Geometry supergeo;
    for (const auto& sys : act) {
      supergeo += (*sys->getGeometry());
    }
    for (const auto& sys : env) {
      supergeo += (*sys->getGeometry());
    }
    supergeo.deleteIdenticalAtoms();
    gaugeOrigin = supergeo.getCenterOfMass();
  }
  else {
    // Otherwise, do nothing and keep input gauge origin
    printf(" Found custom gauge-origin in the input.\n");
  }

  // Print gauge-origin
  printf(" Gauge-origin (Angstrom): %8.4f %8.4f %8.4f\n\n", gaugeOrigin.getX() * BOHR_TO_ANGSTROM,
         gaugeOrigin.getY() * BOHR_TO_ANGSTROM, gaugeOrigin.getZ() * BOHR_TO_ANGSTROM);

  return gaugeOrigin;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<std::vector<Eigen::MatrixXd>>
LRSCFSetup<SCFMode>::setupFDEcTransformation(const LRSCFTaskSettings& settings, const Eigen::MatrixXi& couplingPatternMatrix,
                                             const std::vector<Options::LRSCF_TYPE>& referenceLoadingType,
                                             const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                             const std::vector<std::shared_ptr<SystemController>>& act,
                                             const unsigned nDimension) {
  // Set Up Matrix Size:
  // Number of excitation from each previous calculation determines the number of cols for transformation
  std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors = nullptr;
  unsigned int iStart = 0;
  unsigned int nExcPrev = 0;
  if (settings.uncoupledSubspace.empty()) {
    assert(referenceLoadingType.size() == (unsigned int)couplingPatternMatrix.rows());

    for (unsigned int i = 0; i < lrscf.size(); i++) {
      if (i == 0) {
        if (lrscf[i]->getExcitationEnergies(referenceLoadingType[i])) {
          nExcPrev += (*lrscf[i]->getExcitationEnergies(referenceLoadingType[i])).rows();
        }
        else {
          throw SerenityError("You tried to perform a coupled sLRSCF calculation but Serenity could not find FDEu,"
                              " isolated or coupled solution vectors for all subsystems.");
        }
      }
      else {
        if (referenceLoadingType[i] == Options::LRSCF_TYPE::COUPLED && couplingPatternMatrix(i - 1, i) != 0) {
          continue;
        }
        if (lrscf[i]->getExcitationEnergies(referenceLoadingType[i])) {
          nExcPrev += (*lrscf[i]->getExcitationEnergies(referenceLoadingType[i])).rows();
        }
        else {
          throw SerenityError("You tried to perform a coupled sLRSCF calculation but Serenity could not find FDEu,"
                              " isolated or coupled solution vectors for all subsystems.");
        }
      }
    }
  }
  else {
    nExcPrev = settings.uncoupledSubspace.size() - act.size();
  }

  // Initialize new eigenvectors
  eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(2);
  (*eigenvectors)[0].resize(nDimension, nExcPrev);
  (*eigenvectors)[0].setZero();
  (*eigenvectors)[1].resize(nDimension, nExcPrev);
  (*eigenvectors)[1].setZero();
  unsigned int iCountRows = 0;
  unsigned int iCountCols = 0;
  iStart = 0;
  for (unsigned int I = 0; I < referenceLoadingType.size(); ++I) {
    auto type = referenceLoadingType[I];
    std::shared_ptr<std::vector<Eigen::MatrixXd>> vecI;
    if (lrscf[I]->getExcitationVectors(type)) {
      vecI = lrscf[I]->getExcitationVectors(type);
    }
    else {
      throw SerenityError("You tried to perform a coupled sLRSCF calculation but Serenity could not find FDEu"
                          " or isolated solution vectors for all subsystems.");
    }
    if (settings.uncoupledSubspace.empty()) {
      // In case of coupeld calculations the coloumn does not need to be increased because
      // Carefull when two coupled calculations are coupled then the col needs to be increased
      // e.g.: 1 1 0 0
      //      2 2 0 0
      //      0 0 3 3
      //      0 0 4 4
      // Therefore make sure that the entry in the coupling matrix under the acutal position is not equal to zero
      if (I != 0) {
        if (referenceLoadingType[I] == Options::LRSCF_TYPE::COUPLED && couplingPatternMatrix(I - 1, I) != 0) {
          iCountRows += (*vecI)[0].rows();
          iCountCols = (*vecI)[0].cols();
        }
        else {
          iCountRows += (*vecI)[0].rows();
          iCountCols += (*vecI)[0].cols();
        }
      }
      else {
        iCountRows += (*vecI)[0].rows();
        iCountCols += (*vecI)[0].cols();
      }

      (*eigenvectors)[0].block(iCountRows - (*vecI)[0].rows(), iCountCols - (*vecI)[0].cols(), (*vecI)[0].rows(),
                               (*vecI)[0].cols()) = (*vecI)[0] + (*vecI)[1];
      if (!(settings.tda)) {
        (*eigenvectors)[1].block(iCountRows - (*vecI)[1].rows(), iCountCols - (*vecI)[1].cols(), (*vecI)[1].rows(),
                                 (*vecI)[1].cols()) = (*vecI)[0] - (*vecI)[1];
        ;
      }
    }
    else {
      iCountRows += (*vecI)[0].rows();
      unsigned int nStates = settings.uncoupledSubspace[iStart];
      iStart += 1;
      for (unsigned int iState = 0; iState < nStates; ++iState) {
        (*eigenvectors)[0].block(iCountRows - (*vecI)[0].rows(), iCountCols, (*vecI)[0].rows(), 1) =
            (*vecI)[0].col(settings.uncoupledSubspace[iStart + iState] - 1) +
            (*vecI)[1].col(settings.uncoupledSubspace[iStart + iState] - 1);
        if (!(settings.tda)) {
          (*eigenvectors)[1].block(iCountRows - (*vecI)[1].rows(), iCountCols, (*vecI)[1].rows(), 1) =
              (*vecI)[0].col(settings.uncoupledSubspace[iStart + iState] - 1) -
              (*vecI)[1].col(settings.uncoupledSubspace[iStart + iState] - 1);
        }
        iCountCols += 1;
      }
      iStart += nStates;
    }
  }
  return eigenvectors;
}

template<Options::SCF_MODES SCFMode>
void LRSCFSetup<SCFMode>::setupLRSCFController(const LRSCFTaskSettings& settings, const Eigen::MatrixXi& couplingPatternMatrix,
                                               const std::vector<std::shared_ptr<SystemController>>& act,
                                               const std::vector<std::shared_ptr<SystemController>>& env,
                                               const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscfAll,
                                               const Options::LRSCF_TYPE type) {
  unsigned int iSystem = 0;
  for (auto lrscf : lrscfAll) {
    // AUTOMATICALLY: Read orbital space belonging to the individual subsystems from the uncoupled calculations
    if (type == Options::LRSCF_TYPE::COUPLED) {
      printf(" ------- Read Reference Orbitals -------\n");
      try {
        std::string mode = (SCFMode == RESTRICTED) ? "res" : "unres";
        std::string type_str = "";
        // If a regular coupled calculation is performed then read uncoupled reference orbitals
        // If the lrscfController is part of an coupled uncoupled calculation the coupled vectors and uncoupled vectors
        // are read in
        if ((unsigned int)couplingPatternMatrix.rows() > act.size() && iSystem < act.size()) {
          type_str = ".fdec";
        }
        else {
          type_str = "";
        }
        std::string fName =
            lrscf->getSys()->getSystemPath() + lrscf->getSys()->getSystemName() + type_str + ".lrscfSpace." + mode + ".h5";
        HDF5::Filepath name(fName);
        HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        SpinPolarizedData<SCFMode, std::vector<unsigned int>> indexWhiteList;
        unsigned int iCount = 0;
        for_spin(indexWhiteList) {
          std::string spin = "alpha";
          if (iCount > 0)
            spin = "beta";
          iCount += 1;
          HDF5::dataset_exists(file, spin);
          Eigen::VectorXi tmp;
          HDF5::load(file, spin, tmp);
          for (unsigned int i = 0; i < tmp.size(); ++i) {
            indexWhiteList_spin.push_back(tmp(i));
          }
        };
        file.close();
        lrscf->editReference(indexWhiteList, iSystem, type);
      }
      catch (...) {
        // nothing to be done here
      }
    }

    // Initialize index array
    auto indexWhiteList = lrscf->getReferenceOrbitals();
    // For reference
    auto oldIndexWhiteList = indexWhiteList;

    // Exclude Projection:
    // Exclude artificial projected orbitals from the occupied environment
    // orbitals into the virtual orbital space of the active subsystems
    if (settings.excludeProjection && type != Options::LRSCF_TYPE::ISOLATED) {
      printf("  ------- Exclude Projection -------\n");
      // calculate overlap
      auto& libint = Libint::getInstance();
      auto nOccupied = lrscf->getNOccupied();
      auto nVirtual = lrscf->getNVirtual();
      auto coefAct = lrscf->getCoefficients();

      for (unsigned int iSysAct = 0; iSysAct < lrscfAll.size() + env.size(); iSysAct++) {
        std::shared_ptr<SystemController> envSys = nullptr;
        if (iSysAct != iSystem && iSysAct < lrscfAll.size())
          envSys = lrscfAll[iSysAct]->getSys();
        if (iSysAct >= lrscfAll.size())
          envSys = env[iSysAct - lrscfAll.size()];
        if (envSys == nullptr) {
          // DO nothing
        }
        else {
          // First Basiscontroller --> Col entries in Overlap
          // Second Basiscontroller --> Row entries in Overlap
          // Order therefore here important
          auto overlapAB =
              libint.compute1eInts(libint2::Operator::overlap, lrscf->getBasisController(), envSys->getBasisController());
          std::vector<Eigen::MatrixXd> overlapMO;
          auto nOccEnv = envSys->getNOccupiedOrbitals<SCFMode>();
          auto coefEnv = envSys->getActiveOrbitalController<SCFMode>()->getCoefficients();

          // Calculate overlap Sia between occ env and virt act
          for_spin(nOccEnv, nVirtual, nOccupied, coefEnv, coefAct) {
            auto nBasisEnv = envSys->getBasisController()->getNBasisFunctions();
            auto nBasisAct = lrscf->getBasisController()->getNBasisFunctions();
            overlapMO.push_back(coefEnv_spin.block(0, 0, nBasisEnv, nOccEnv_spin).transpose() * overlapAB *
                                coefAct_spin.block(0, nOccupied_spin, nBasisAct, nVirtual_spin));
          };
          // Find orbitals of env system through overlap criteria
          unsigned int spinCounter = 0;
          for_spin(indexWhiteList, nOccupied, oldIndexWhiteList) {
            for (unsigned int col = 0; col < overlapMO[spinCounter].cols(); col++) {
              double sum = 0.0;
              for (unsigned int row = 0; row < overlapMO[spinCounter].rows(); row++) {
                sum += std::fabs(overlapMO[spinCounter](row, col));
              }
              // remove Orbital from Index
              if (sum > 0.95) {
                indexWhiteList_spin.erase(std::remove(indexWhiteList_spin.begin(), indexWhiteList_spin.end(),
                                                      oldIndexWhiteList_spin[col + nOccupied_spin]),
                                          indexWhiteList_spin.end());
              }
            }
            spinCounter += 1;
          };
          // edit Reference Orbitals
          lrscf->editReference(indexWhiteList, iSystem, type);
        }
      }
    }
    iSystem += 1;
  }
}

template class LRSCFSetup<Options::SCF_MODES::RESTRICTED>;
template class LRSCFSetup<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
