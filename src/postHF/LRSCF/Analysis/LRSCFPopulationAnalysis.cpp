/**
 * @file LRSCFPopulationAnalysis.cpp
 * @author Niklas Niemeyer, Anton Rikus
 *
 * @date October 10, 2021
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
#include "postHF/LRSCF/Analysis/LRSCFPopulationAnalysis.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/LoewdinPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutputStream.h"
#include "math/Matrix.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void LRSCFPopulationAnalysis<SCFMode>::calculateTransitionCharges(
    const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf, std::vector<Eigen::MatrixXd>& densities) {
  Timings::takeTime("LRSCF -    Transition Charges");

  unsigned nEigen = densities[0].cols();

  long iStartSys = 0;
  for (unsigned I = 0; I < lrscf.size(); ++I) {
    // Get system info.
    auto no = lrscf[I]->getNOccupied();
    auto nv = lrscf[I]->getNVirtual();
    auto P = lrscf[I]->getParticleCoefficients();
    auto H = lrscf[I]->getHoleCoefficients();
    auto cc2 = lrscf[I]->getCC2Controller();

    unsigned nAtoms = lrscf[I]->getSys()->getNAtoms();
    const auto& atoms = lrscf[I]->getSys()->getAtoms();
    const auto& ints = lrscf[I]->getSys()->getOneElectronIntegralController();
    const auto& basis = lrscf[I]->getSys()->getAtomCenteredBasisController();

    // Prepare populations for this subsystem (also, store coordinates in that file for easier access).
    Eigen::MatrixXd pops = Eigen::MatrixXd::Zero(nAtoms, 3 + nEigen);
    pops.leftCols(3) = lrscf[I]->getSys()->getGeometry()->getCoordinates();

    Eigen::MatrixXd particlepops = Eigen::MatrixXd::Zero(nAtoms, 3 + nEigen);
    Eigen::MatrixXd holepops = Eigen::MatrixXd::Zero(nAtoms, 3 + nEigen);

    // Perform MO/AO transformation and calculate transition charges.
    long iStartSpin = 0;
    printSmallCaption("Löwdin Transition Density Populations");
    printf("  State    Atom    Atomtype    Transition    Hole    Particle\n");
    for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
      DensityMatrix<RESTRICTED> densityAO(basis);
      DensityMatrix<RESTRICTED> holeAO(basis);
      DensityMatrix<RESTRICTED> particleAO(basis);

      iStartSpin = 0;
      for_spin(P, H, nv, no) {
        unsigned nb = nv_spin + no_spin;
        unsigned np = cc2 ? nb : nv_spin;
        unsigned nh = cc2 ? nb : no_spin;
        Eigen::Map<Eigen::MatrixXd> xpy(densities[0].col(iEigen).data() + iStartSys + iStartSpin, np, nh);
        Eigen::Map<Eigen::MatrixXd> xmy(densities[1].col(iEigen).data() + iStartSys + iStartSpin, np, nh);
        // Transition density
        densityAO += P_spin.middleCols(cc2 ? 0 : nh, np) * xpy * H_spin.leftCols(nh).transpose();
        // Hole density.
        holeAO -=
            0.5 * H_spin.leftCols(nh) * (xpy.transpose() * xpy + xmy.transpose() * xmy) * H_spin.leftCols(nh).transpose();
        // Particle density.
        particleAO += 0.5 * P_spin.middleCols(cc2 ? 0 : nh, np) * (xpy * xpy.transpose() + xmy * xmy.transpose()) *
                      P_spin.middleCols(cc2 ? 0 : nh, np).transpose();
        iStartSpin += np * nh;
      };
      pops.col(3 + iEigen) = LoewdinPopulationCalculator<RESTRICTED>::calculateAtomPopulations(
          densityAO, ints->getOverlapIntegrals(), basis->getBasisIndices());
      holepops.col(iEigen) = LoewdinPopulationCalculator<RESTRICTED>::calculateAtomPopulations(
          holeAO, ints->getOverlapIntegrals(), basis->getBasisIndices());
      particlepops.col(iEigen) = LoewdinPopulationCalculator<RESTRICTED>::calculateAtomPopulations(
          particleAO, ints->getOverlapIntegrals(), basis->getBasisIndices());

      printf("--------------------------------------------------------------\n");
      for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
        printf("%5i %8i %10s %14.4f %9.4f %9.4f\n", (iEigen + 1), (iAtom + 1), atoms[iAtom]->getAtomType()->getName().c_str(),
               pops(iAtom, 3 + iEigen), holepops(iAtom, iEigen), particlepops(iAtom, iEigen));
      }
      printf("--------------------------------------------------------------\n");
      printf("%40.4f %9.4f %9.4f\n", pops.col(3 + iEigen).sum(), holepops.col(iEigen).sum(), particlepops.col(iEigen).sum());

      // The correlation matrix is reversed column-wise so that the atom order is reversed.
      // This is purely a preference. The probability that a particle state is at atom 1 while the hole state is at atom
      // 1 as well is in the lower left corner. Particle states are on the y-axis while hole states are on the x-axis.
      std::ofstream file(lrscf[I]->getSys()->getSystemPath() + lrscf[I]->getSys()->getSystemName() + ".correlation" +
                         std::to_string(iEigen + 1) + ".txt");
      Eigen::MatrixXd correlation = -particlepops.col(iEigen) * holepops.col(iEigen).transpose();
      if (file.is_open()) {
        file << (correlation.colwise().reverse());
      }
    }
    iStartSys += iStartSpin;

    // Store on disk.
    if (SCFMode == Options::SCF_MODES::RESTRICTED) {
      OutputControl::nOut << "\n  Transition charges written to disk are multiplied with sqrt(2).\n" << std::endl;
      pops.rightCols(nEigen) *= std::sqrt(2);
    }
    std::string fileName = lrscf[I]->getSys()->getSystemPath() + lrscf[I]->getSys()->getSystemName() + ".transitioncharges.txt";
    std::ofstream file(fileName);
    if (file.is_open()) {
      file << std::scientific << std::setprecision(16) << pops << "\n";
    }
    file.close();
  }
  Timings::timeTaken("LRSCF -    Transition Charges");
}

template class LRSCFPopulationAnalysis<Options::SCF_MODES::RESTRICTED>;
template class LRSCFPopulationAnalysis<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
