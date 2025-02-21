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
#include "io/FormattedOutputStream.h"
#include "io/HDF5.h"
#include "misc/WarningTracker.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Tools/Besley.h"
#include "settings/LRSCFOptions.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "tasks/VirtualOrbitalSpaceSelectionTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void LRSCFSetup<SCFMode>::printInfo(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                    const LRSCFTaskSettings& settings,
                                    const std::vector<std::shared_ptr<SystemController>>& envSys,
                                    const Options::LRSCF_TYPE type) {
  // Some output (Information of the calculation)
  // Dummy just needed to compile for_spin macro

  SpinPolarizedData<SCFMode, int> dummy;

  printBigCaption("LRSCF Info");

  printf(" System(s):\n");
  printf("%5s %20s %5s ", "", "Sys", "Stat");
  unsigned iStart = 0;
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
  for (unsigned I = 0; I < lrscf.size(); ++I) {
    printf("%5i ", I + 1);
    SpinPolarizedData<SCFMode, unsigned> nOccupied = lrscf[I]->getNOccupied();
    SpinPolarizedData<SCFMode, unsigned> nVirtual = lrscf[I]->getNVirtual();
    printf("%20s ", lrscf[I]->getSys()->getSystemName().c_str());
    printf("%5s ", "act");
    for_spin(nOccupied) {
      printf("%8i ", nOccupied_spin);
    };
    for_spin(nVirtual) {
      printf("%8i ", nVirtual_spin);
    };
    unsigned nDimI = 0;
    for_spin(nOccupied, nVirtual) {
      nDimI += nOccupied_spin * nVirtual_spin;
    };
    printf("%8i \n", nDimI);
  }
  for (unsigned I = 0; I < envSys.size(); ++I) {
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
    unsigned nDimI = 0;
    for_spin(nOccupied, nVirtual) {
      nDimI += nOccupied_spin * nVirtual_spin;
    };
    printf("%8i \n", nDimI);
  }

  printf("\n Type        : %22s \n", (type == Options::LRSCF_TYPE::ISOLATED)    ? "Isolated"
                                     : (type == Options::LRSCF_TYPE::UNCOUPLED) ? "FDEu"
                                                                                : "FDEc");
  std::string naddXCFunc_string;
  std::string naddKinFunc_string;
  std::string func;
  auto naddkin = settings.embedding.naddKinFunc;
  auto xc = settings.func;
  Options::resolve<CompositeFunctionals::KINFUNCTIONALS>(naddKinFunc_string, naddkin);
  Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(func, xc);
  printf(" NaddKinFunc : %22s \n", naddKinFunc_string.c_str());
  if (settings.embedding.embeddingModeList.size() > 0) {
    OutputControl::mOut << " -------------------------------------------- " << std::endl;
    OutputControl::mOut << " NaddXCFunc used for approx./exact embedding: " << std::endl;
    auto naddxc = settings.embedding.naddXCFuncList[0];
    Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(naddXCFunc_string, naddxc);
    printf("  NaddXCFunc for exact   : %22s \n", naddXCFunc_string.c_str());
    std::string naddXCFunc_string2;
    naddxc = settings.embedding.naddXCFuncList[1];
    Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(naddXCFunc_string2, naddxc);
    printf("  NaddXCFunc for approx. : %22s \n", naddXCFunc_string2.c_str());
  }
  else {
    auto naddxc = settings.embedding.naddXCFunc;
    Options::resolve<CompositeFunctionals::XCFUNCTIONALS>(naddXCFunc_string, naddxc);
    printf(" NaddXCFunc  : %22s \n\n", naddXCFunc_string.c_str());
  }
  if (settings.func != CompositeFunctionals::XCFUNCTIONALS::NONE)
    printf(" Func        : %22s \n\n", func.c_str());
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd LRSCFSetup<SCFMode>::getDiagonal(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf) {
  unsigned iStart = 0;
  Eigen::VectorXd diagonal(0);
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    auto e = ilrscf->getEigenvalues();
    for_spin(no, nv, e) {
      unsigned nvno = nv_spin * no_spin;
      diagonal.conservativeResize(diagonal.size() + nvno);
      for (unsigned ia = 0; ia < no_spin * nv_spin; ++ia) {
        unsigned i = ia / nv_spin;
        unsigned a = no_spin + ia % nv_spin;
        diagonal(iStart + ia) = e_spin(a) - e_spin(i);
      }
      iStart += nvno;
    };
  }
  return diagonal;
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
    printf("\n  Set gauge-origin to center of mass.\n");
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
    printf("  Found custom gauge-origin in the input.\n");
  }

  // Print gauge-origin
  printf("  Gauge-origin (Angstrom): %8.4f %8.4f %8.4f\n\n", gaugeOrigin.getX() * BOHR_TO_ANGSTROM,
         gaugeOrigin.getY() * BOHR_TO_ANGSTROM, gaugeOrigin.getZ() * BOHR_TO_ANGSTROM);

  return gaugeOrigin;
}

template<Options::SCF_MODES SCFMode>
void LRSCFSetup<SCFMode>::calculateMolecularWeightandNuclearDipoleMoment(
    const std::vector<std::shared_ptr<SystemController>>& act, const std::vector<std::shared_ptr<SystemController>>& env,
    Point gaugeOrigin, double& molWeight, Eigen::Vector3d& nucDipoleMoment) {
  molWeight = 0.0;
  nucDipoleMoment.setZero();
  for (const auto& sys : act) {
    auto atoms = sys->getGeometry()->getAtoms();
    for (const auto& atom : atoms) {
      if (!atom->isDummy()) {
        molWeight += atom->getAtomType()->getMass();
        auto mom = (Point(*atom) - gaugeOrigin) * atom->getEffectiveCharge();
        nucDipoleMoment(0) += mom.getX();
        nucDipoleMoment(1) += mom.getY();
        nucDipoleMoment(2) += mom.getZ();
      }
    }
  }
  printf("  Molecular weight: %10.4f\n\n", molWeight);
  printf("  Nuclear dipole moment: %10.4f %10.4f %10.4f\n\n", nucDipoleMoment(0), nucDipoleMoment(1), nucDipoleMoment(2));
  if (env.size() > 0) {
    for (const auto& sys : env) {
      auto atoms = sys->getGeometry()->getAtoms();
      for (const auto& atom : atoms) {
        if (!atom->isDummy()) {
          molWeight += atom->getAtomType()->getMass();
          auto mom = (Point(*atom) - gaugeOrigin) * atom->getEffectiveCharge();
          nucDipoleMoment(0) += mom.getX();
          nucDipoleMoment(1) += mom.getY();
          nucDipoleMoment(2) += mom.getZ();
        }
      }
    }
    printf("  Molecular weight with env. systems: %10.4f\n\n", molWeight);
    printf("  Nuclear dipole moment with env. systems: %10.4f %10.4f %10.4f\n\n", nucDipoleMoment(0),
           nucDipoleMoment(1), nucDipoleMoment(2));
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<std::vector<Eigen::MatrixXd>>
LRSCFSetup<SCFMode>::setupFDEcTransformation(Eigen::VectorXd& eigenvalues, const LRSCFTaskSettings& settings,
                                             const Eigen::MatrixXi& couplingPatternMatrix,
                                             const std::vector<Options::LRSCF_TYPE>& referenceLoadingType,
                                             const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                             const std::vector<std::shared_ptr<SystemController>>& act,
                                             const unsigned nDimension) {
  printBigCaption("Assembling FDEc Transformation Matrix");
  // Set Up Matrix Size:
  // Number of excitation from each previous calculation determines the number of cols for transformation
  std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors = nullptr;
  unsigned iStart = 0;
  unsigned nExcPrev = 0;
  if (settings.uncoupledSubspace.empty()) {
    assert(referenceLoadingType.size() == (unsigned)couplingPatternMatrix.rows());

    for (unsigned i = 0; i < lrscf.size(); i++) {
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
  eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(2, Eigen::MatrixXd::Zero(nDimension, nExcPrev));
  eigenvalues = Eigen::VectorXd::Zero(nExcPrev);

  unsigned iCountRows = 0;
  unsigned iCountCols = 0;
  iStart = 0;
  for (unsigned I = 0; I < referenceLoadingType.size(); ++I) {
    auto type = referenceLoadingType[I];
    std::shared_ptr<std::vector<Eigen::MatrixXd>> vecI;
    Eigen::VectorXd eigenvalues_I;
    if (lrscf[I]->getExcitationVectors(type)) {
      vecI = lrscf[I]->getExcitationVectors(type);
      eigenvalues_I = (*lrscf[I]->getExcitationEnergies(type));
    }
    else {
      throw SerenityError("You tried to perform a coupled sLRSCF calculation but Serenity could not find FDEu"
                          " or isolated solution vectors for all subsystems.");
    }
    Eigen::MatrixXd& R_I = (*vecI)[0];
    Eigen::MatrixXd& L_I = (*vecI)[1];
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
          iCountRows += R_I.rows();
          iCountCols = R_I.cols();
        }
        else {
          iCountRows += R_I.rows();
          iCountCols += R_I.cols();
        }
      }
      else {
        iCountRows += R_I.rows();
        iCountCols += R_I.cols();
      }

      (*eigenvectors)[0].block(iCountRows - R_I.rows(), iCountCols - R_I.cols(), R_I.rows(), R_I.cols()) = R_I;
      (*eigenvectors)[1].block(iCountRows - L_I.rows(), iCountCols - L_I.cols(), L_I.rows(), L_I.cols()) = L_I;
      // Todo: Add this for the uncoupled subspace.
      eigenvalues.segment(iCountCols - R_I.cols(), eigenvalues_I.size()) = eigenvalues_I;
    }
    else {
      iCountRows += R_I.rows();
      unsigned nStates = settings.uncoupledSubspace[iStart];
      iStart += 1;
      for (unsigned iState = 0; iState < nStates; ++iState) {
        unsigned index = settings.uncoupledSubspace[iStart + iState] - 1;
        (*eigenvectors)[0].block(iCountRows - R_I.rows(), iCountCols, R_I.rows(), 1) = R_I.col(index);
        (*eigenvectors)[1].block(iCountRows - L_I.rows(), iCountCols, L_I.rows(), 1) = L_I.col(index);
        iCountCols += 1;
      }
      iStart += nStates;
    }
  }

  return eigenvectors;
}

template<Options::SCF_MODES SCFMode>
void LRSCFSetup<SCFMode>::setupLRSCFController(const LRSCFTaskSettings& settings,
                                               const std::vector<std::shared_ptr<SystemController>>& act,
                                               const std::vector<std::shared_ptr<SystemController>>& env,
                                               const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscfAll) {
  for (auto lrscf : lrscfAll) {
    // Update everything
    auto sys = lrscf->getSys();
    if (settings.excludeProjection) {
      if (env.size() == 0 && act.size() == 1) {
        throw SerenityError("You need to specify an environment system for exclude projection!");
      }
      // System vector
      std::vector<std::shared_ptr<SystemController>> remainingSys;
      for (auto i : act) {
        if (i != sys) {
          remainingSys.push_back(i);
        }
      }
      for (auto i : env) {
        remainingSys.push_back(i);
      }
      VirtualOrbitalSpaceSelectionTask<SCFMode> voss({sys}, remainingSys);
      voss.settings.excludeProjection = true;
      voss.run();
      // Update everything
      auto coefs = sys->template getActiveOrbitalController<SCFMode>()->getCoefficients();
      auto orbitalEner = sys->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
      auto nOcc = sys->template getNOccupiedOrbitals<SCFMode>();
      auto nVirt = sys->template getNVirtualOrbitalsTruncated<SCFMode>();

      lrscf->setNOccupied(nOcc);
      lrscf->setNVirtual(nVirt);
      lrscf->setCoefficients(coefs);
      lrscf->setEigenvalues(orbitalEner);
    }
    // Restrict Orbitals according to Besley's criterion for occupied and virtual orbitals
    if (settings.besleyAtoms > 0) {
      if (settings.besleyCutoff.size() != 2)
        throw SerenityError("Keyword besleyCutoff needs two arguments!");
      Besley<SCFMode> besley(sys, settings.besleyAtoms, settings.besleyCutoff);
      auto indexWhiteList = besley.getWhiteList();
      lrscf->editReference(indexWhiteList);
    }
    if (settings.frozenCore) {
      lrscf->applyFrozenCore();
    }
    if (settings.coreOnly) {
      lrscf->applyCoreOnly();
    }
    if (settings.frozenVirtual != 0.0) {
      lrscf->applyFrozenVirtual();
    }
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFSetup<SCFMode>::prepareStabilityAnalysis(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                   LRSCFTaskSettings& settings) {
  // Must enforce symmetric algorithms and the TDA lambda function for calculating the matrix--vector products.
  settings.method = Options::LR_METHOD::TDA;

  if (lrscf.size() > 1) {
    WarningTracker::printWarning("Stability analysis only supported for isolated calculations. Use at your own risk.", true);
  }

  if (settings.scfstab == Options::STABILITY_ANALYSIS::REAL) {
    if (settings.triplet && SCFMode == Options::SCF_MODES::RESTRICTED) {
      printBigCaption("Performing RHF triplet (external) stability analysis: (A+B), triplet");
    }
    else {
      if (SCFMode == Options::SCF_MODES::RESTRICTED) {
        printBigCaption("Performing RHF singlet (internal) stability analysis: (A+B), singlet");
      }
      else {
        printBigCaption("Performing UHF singlet (internal) stability analysis: (A+B)");
      }
    }
  }
  else if (settings.scfstab == Options::STABILITY_ANALYSIS::NONREAL) {
    if (SCFMode == Options::SCF_MODES::RESTRICTED) {
      printBigCaption("Performing RHF non-real (external) stability analysis: (A-B)");
    }
    else {
      printBigCaption("Performing UHF non-real (external) stability analysis: (A-B)");
    }
  }
  else if (settings.scfstab == Options::STABILITY_ANALYSIS::SPINFLIP) {
    printBigCaption("Performing Spin-Flip TDDFT/TDHF within the TDA");
    for (auto& ilrscf : lrscf) {
      ilrscf->setupSpinFlipReference();
    }
  }
  else {
    throw SerenityError("Choose an SCF stability analysis.");
  }
}

template<Options::SCF_MODES SCFMode>
void LRSCFSetup<SCFMode>::printApproximateCoulombInfo(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                      LRSCFTaskSettings& settings) {
  if (!(settings.approxCoulomb.size() == 2 || settings.approxCoulomb.size() == 1)) {
    throw SerenityError("approxCoulomb expects one or two thresholds!");
  }
  if (settings.approxCoulomb.size() == 1) {
    settings.approxCoulomb.push_back(std::numeric_limits<double>::infinity());
  }
  if (settings.approxCoulomb[0] > settings.approxCoulomb[1]) {
    throw SerenityError("The second entry of approxCoulomb must be equal or greater than the first!");
  }

  if (settings.approxCoulomb[0] != std::numeric_limits<double>::infinity() ||
      settings.approxCoulomb[1] != std::numeric_limits<double>::infinity()) {
    printBigCaption("Approximate Coulomb Interaction");
    printf("  Distance Thresholds (Angstrom): %4.2f | %-4.2f\n\n", settings.approxCoulomb[0], settings.approxCoulomb[1]);

    // Work with Bohr internally.
    settings.approxCoulomb[0] *= ANGSTROM_TO_BOHR;
    settings.approxCoulomb[1] *= ANGSTROM_TO_BOHR;

    unsigned nInteractions = (lrscf.size() * (lrscf.size() - 1) / 2);

    unsigned nRICouplings = 0;
    unsigned nMRICouplings = 0;
    unsigned nGrimmeCouplings = 0;

    for (unsigned I = 0; I < lrscf.size(); ++I) {
      for (unsigned J = I + 1; J < lrscf.size(); ++J) {
        double distance = lrscf[I]->getSys()->getGeometry()->getMinimumDistance(*lrscf[J]->getSys()->getGeometry());
        if (distance < settings.approxCoulomb[0]) {
          nRICouplings = nRICouplings + 1;
        }
        if (distance >= settings.approxCoulomb[0] && distance < settings.approxCoulomb[1]) {
          nMRICouplings = nMRICouplings + 1;
        }
        if (distance >= settings.approxCoulomb[1]) {
          nGrimmeCouplings = nGrimmeCouplings + 1;
        }
      }
    }
    if (nRICouplings + nMRICouplings + nGrimmeCouplings != nInteractions) {
      throw SerenityError("Something is wrong with your approximate Coulomb interactions!");
    }
    printf("  Number of regular Coulomb interactions (RI or NORI) : %6i / %-6i\n", nRICouplings, nInteractions);
    printf("  Number of approximate Coulomb interactions (mon RI) : %6i / %-6i\n", nMRICouplings, nInteractions);
    printf("                                             (Grimme) : %6i / %-6i\n\n", nGrimmeCouplings, nInteractions);
  }
}

template class LRSCFSetup<Options::SCF_MODES::RESTRICTED>;
template class LRSCFSetup<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
