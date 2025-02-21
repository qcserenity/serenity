/**
 * @file CouplingConstruction.cpp
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
#include "postHF/LRSCF/Tools/CouplingConstruction.h"
/* Include Serenity Internal Headers */
#include "math/linearAlgebra/Orthogonalization.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/LRSCFOptions.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void CouplingConstruction<SCFMode>::solve(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                          const LRSCFTaskSettings& settings,
                                          std::vector<Options::LRSCF_TYPE> referenceLoadingType, SigmaCalculator sigmaCalculator,
                                          std::shared_ptr<std::vector<Eigen::MatrixXd>> eigenvectors,
                                          Eigen::VectorXd& eigenvalues) {
  printBigCaption("Coupling-Matrix Construction");
  if (lrscf.size() < 2) {
    throw SerenityError("Partial response-matrix construction only for more than one subsystem");
  }

  unsigned nEigen = (*eigenvectors)[0].cols();

  Eigen::VectorXi nEigenSub = Eigen::VectorXi::Zero(lrscf.size());
  Eigen::VectorXi nToSubsystem = Eigen::VectorXi::Zero(nEigen);
  Eigen::VectorXi nToExcitation = Eigen::VectorXi::Zero(nEigen);
  for (unsigned I = 0, n = 0; I < lrscf.size(); ++I) {
    nEigenSub[I] = (*lrscf[I]->getExcitationVectors(referenceLoadingType[I]))[0].cols();
    for (int iEigenI = 0; iEigenI < nEigenSub[I]; ++iEigenI, ++n) {
      nToSubsystem[n] = I;
      nToExcitation[n] = iEigenI;
    }
  }

  if (settings.method == Options::LR_METHOD::TDA) {
    printf("  Using subsystem TDA in coupling-matrix construction!\n\n");
    std::vector<Eigen::MatrixXd> guessVectors = {(*eigenvectors)[0]};
    Eigen::MatrixXd partialSigmaVectors = (*sigmaCalculator(guessVectors))[0];
    Eigen::MatrixXd responseMatrix = guessVectors[0].transpose() * partialSigmaVectors;

    // Fill diagonal blocks with old excitation energies.
    unsigned iStart = 0;
    for (unsigned I = 0; I < lrscf.size(); ++I) {
      Eigen::VectorXd energies = *(lrscf[I]->getExcitationEnergies(referenceLoadingType[I]));
      responseMatrix.block(iStart, iStart, energies.size(), energies.size()) = energies.asDiagonal();
      iStart += energies.size();
    }

    // Copy upper triangular into lower triangular.
    responseMatrix.triangularView<Eigen::StrictlyLower>() = responseMatrix.triangularView<Eigen::StrictlyUpper>().transpose();

    // Write subspace response matrix to disk if requested.
    std::string name = lrscf[0]->getSys()->getSystemName();
    for (unsigned I = 1; I < lrscf.size(); ++I) {
      name += "_" + lrscf[I]->getSys()->getSystemName();
    }
    std::ofstream file(name + "_FDEcMatrix.txt");
    file << std::scientific << std::setprecision(16) << responseMatrix;
    file.close();

    printSmallCaption("Coupling Matrix / eV");
    for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
      printf(" ");
      for (unsigned jEigen = 0; jEigen < nEigen; ++jEigen) {
        printf("%12.3e", responseMatrix(iEigen, jEigen) * HARTREE_TO_EV);
      }
      printf("\n");
    }
    printf("\n");

    // Save coupling blocks to disk for later analysis.
    for (unsigned I = 0; I < lrscf.size(); ++I) {
      unsigned offI = nEigenSub.segment(0, I).sum();
      unsigned nEigenI = nEigenSub[I];
      std::string pathI = lrscf[I]->getSys()->getSystemPath();
      std::string nameI = lrscf[I]->getSys()->getSystemName();
      for (unsigned J = I + 1; J < lrscf.size(); ++J) {
        unsigned offJ = nEigenSub.segment(0, J).sum();
        unsigned nEigenJ = nEigenSub[J];
        std::string pathJ = lrscf[J]->getSys()->getSystemPath();
        std::string nameJ = lrscf[J]->getSys()->getSystemName();

        std::ofstream fileI(pathI + nameJ + ".TDACoupling.txt");
        std::ofstream fileJ(pathJ + nameI + ".TDACoupling.txt");
        Eigen::MatrixXd IJ_block = responseMatrix.block(offI, offJ, nEigenI, nEigenJ);
        Eigen::MatrixXd JI_block = IJ_block.transpose();

        fileI << std::scientific << std::setprecision(16) << IJ_block;
        fileJ << std::scientific << std::setprecision(16) << JI_block;
      }
    }

    // Subspace solution.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> subspaceSolver(responseMatrix);
    const Eigen::MatrixXd& eigv = subspaceSolver.eigenvectors();
    (*eigenvectors)[0] = guessVectors[0] * eigv;
    eigenvalues = subspaceSolver.eigenvalues();

    // Print subsystem contributions.
    printf("---------------------------------------------------------------------------------------\n");
    printf("                                Subsystem Contributions                                \n");
    printf("---------------------------------------------------------------------------------------\n");
    printf(" state       energy      wavelength      sys       excitation       contribution       \n");
    printf("              (eV)          (nm)                                      100*|c|^2        \n");
    printf("---------------------------------------------------------------------------------------\n");
    for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
      Eigen::VectorXd contributions = eigv.col(iEigen).cwiseProduct(eigv.col(iEigen));
      bool first = true;
      for (unsigned c = 0; c < nEigen; ++c) {
        // Print only when contribution dominant.
        if (contributions[c] > 0.1) {
          // Print contribution.
          if (first) {
            printf(" %3i %15.5f %12.1f %10i %14i %18.2f\n", iEigen + 1, eigenvalues[iEigen] * HARTREE_TO_EV,
                   HARTREE_TO_NM / eigenvalues[iEigen], nToSubsystem[c] + 1, nToExcitation[c] + 1, 100 * contributions[c]);
            first = false;
          }
          else {
            printf(" %32s %10i %14i %18.2f\n", "", nToSubsystem[c] + 1, nToExcitation[c] + 1, 100 * contributions[c]);
          }
        }
      }
    }
  }
  else if (settings.method == Options::LR_METHOD::TDDFT) {
    printf("  Using subsystem TDDFT in coupling-matrix construction!\n\n");

    std::vector<Eigen::MatrixXd> guessVectors = (*eigenvectors);
    Orthogonalization::modifiedGramSchmidtLinDep(guessVectors, 0.0);
    std::vector<Eigen::MatrixXd> partialSigmaVectors = (*sigmaCalculator(guessVectors));

    // Construct coupling matrix.
    Eigen::MatrixXd responseMatrix = Eigen::MatrixXd::Zero(2 * nEigen, 2 * nEigen);
    responseMatrix.topLeftCorner(nEigen, nEigen) = guessVectors[0].transpose() * partialSigmaVectors[0];
    responseMatrix.bottomRightCorner(nEigen, nEigen) = guessVectors[1].transpose() * partialSigmaVectors[1];

    // Copy upper triangular into lower triangular.
    responseMatrix.triangularView<Eigen::StrictlyLower>() = responseMatrix.triangularView<Eigen::StrictlyUpper>().transpose();

    Eigen::MatrixXd metric = Eigen::MatrixXd::Zero(2 * nEigen, 2 * nEigen);
    metric.topRightCorner(nEigen, nEigen) = guessVectors[0].transpose() * guessVectors[1];
    metric.bottomLeftCorner(nEigen, nEigen) = guessVectors[1].transpose() * guessVectors[0];

    Eigen::FullPivLU<Eigen::MatrixXd> lu(metric);
    if (lu.isInvertible()) {
      responseMatrix = (lu.inverse() * responseMatrix).eval();
    }
    else {
      throw SerenityError("Subspace metric cannot be inverted.");
    }

    Eigen::EigenSolver<Eigen::MatrixXd> subspaceSolver(responseMatrix);
    Eigen::MatrixXd seigenvectors = subspaceSolver.eigenvectors().real();
    Eigen::VectorXd seigenvalues = subspaceSolver.eigenvalues().real();

    // Sort in ascending order.
    unsigned iMin;
    for (unsigned i = 0; i < 2 * nEigen; ++i) {
      seigenvalues.tail(2 * nEigen - i).minCoeff(&iMin);
      seigenvalues.row(i).swap(seigenvalues.row(iMin + i));
      seigenvectors.col(i).swap(seigenvectors.col(iMin + i));
    }

    // Only take positive ones.
    eigenvalues = seigenvalues.tail(nEigen).eval();
    seigenvectors = seigenvectors.rightCols(nEigen).eval();

    // Ritz vectors.
    (*eigenvectors)[0] = guessVectors[0] * seigenvectors.topRows(nEigen);
    (*eigenvectors)[1] = guessVectors[1] * seigenvectors.bottomRows(nEigen);

    // Normalize <X+Y|X-Y> = 1.
    for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
      double clr = (*eigenvectors)[1].col(iEigen).dot((*eigenvectors)[0].col(iEigen));
      (*eigenvectors)[0].col(iEigen) *= 1.0 / std::sqrt(clr);
      (*eigenvectors)[1].col(iEigen) *= 1.0 / std::sqrt(clr);
    }
  }
}

template class CouplingConstruction<Options::SCF_MODES::RESTRICTED>;
template class CouplingConstruction<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
