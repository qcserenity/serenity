/**
 * @file   FXDTask.cpp
 *
 * @date   Mar 3, 2020
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
#include "tasks/FXDTask.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/HDF5.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FXDTask<SCFMode>::FXDTask(const std::shared_ptr<SystemController> systemController)
  : _systemController(systemController) {
}

template<Options::SCF_MODES SCFMode>
void FXDTask<SCFMode>::run() {
  printSectionTitle("FXD");
  // Check input
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "WARNING: Currently this task assumes that the previous calculation was a CIS/TDA calculation!" << std::endl;
  if (settings.donoratoms.size() != 2 || settings.acceptoratoms.size() != 2)
    throw SerenityError("For the definition of the donor/acceptor input two indices need to be specified. The number "
                        "of the first and the number of the last atom belonging to the donor/acceptor molecule");
  if (settings.donoratoms[0] > settings.donoratoms[1] || settings.acceptoratoms[0] > settings.acceptoratoms[1])
    throw SerenityError("Attention: Your first atom index of the donor/accpetor molecule is larger than your second");
  if (settings.loadType == Options::LRSCF_TYPE::COUPLED) {
    throw SerenityError("Coupled excitation vectors are not supported!");
  }
  // Setup LRSCFController
  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.loadType = settings.loadType;
  lrscfSettings.method = Options::LR_METHOD::TDA;
  _lrscf = std::make_shared<LRSCFController<SCFMode>>(_systemController, lrscfSettings);
  // load excitation energies from TDA/CIS calculation
  _excEner = _lrscf->getExcitationEnergies(settings.loadType);
  // load excitation vectors from TDA/CIS calculation
  std::shared_ptr<std::vector<Eigen::MatrixXd>> temp = _lrscf->getExcitationVectors(settings.loadType);
  _excVecs = std::make_shared<Eigen::MatrixXd>((*temp)[0]);

  // If only a specific amount of states should be coupled
  if (settings.states != 100) {
    if (settings.states > (*_excVecs).cols())
      throw SerenityError("You have less excitations calculated than states you want to analyse!");
    Eigen::MatrixXd temp = (*_excVecs).leftCols(settings.states);
    (*_excVecs) = temp;
    Eigen::VectorXd temp2 = (*_excEner).segment(0, settings.states);
    (*_excEner) = temp2;
  }

  // Print the donor and acceptor atoms
  auto geometry = _systemController->getGeometry();
  auto elements = geometry->getAtomSymbols();
  std::cout << "\n Donor atoms:" << std::endl;
  for (unsigned int iAtom = settings.donoratoms[0]; iAtom < settings.donoratoms[1] + 1; iAtom++) {
    std::cout << " " << iAtom << " " << elements[iAtom] << "\n";
  }
  std::cout << " Acceptor atoms:" << std::endl;
  for (unsigned int iAtom = settings.acceptoratoms[0]; iAtom < settings.acceptoratoms[1] + 1; iAtom++) {
    std::cout << " " << iAtom << " " << elements[iAtom] << "\n";
  }

  // Performs Fragment Excitation Difference (FED) approach [J. Chem. Phys. C, 2008, 12, 1204]
  if (settings.FED) {
    auto fed_matrix = calculateFXDMatrix(_excVecs, true);
    printf("\n --------------------- FED Coupling calculation ---------------------- \n");
    printf("     States      X(D)        X(A)         dX12      Coupling(eV)   \n");
    printf(" --------------------------------------------------------------------- \n");

    for (unsigned int n = 0; n < (*_excVecs).cols(); n++) {
      double delta_da_nn = fed_matrix[0](n, n) - fed_matrix[1](n, n);
      for (unsigned int m = n; m < (*_excVecs).cols(); m++) {
        double delta_da_nm = fed_matrix[0](n, m) - fed_matrix[1](n, m);
        double delta_da_mm = fed_matrix[0](m, m) - fed_matrix[1](m, m);
        double numerator = ((*_excEner)(n) - (*_excEner)(m)) * abs(delta_da_nm);
        double coupling = 0.0;
        if (numerator != 0.0) {
          double denominator = sqrt(pow(delta_da_nn - delta_da_mm, 2) + 4.0 * pow(delta_da_nm, 2));
          coupling = (numerator / denominator) * HARTREE_TO_EV;
        }
        printf("  %3i  %3i   %10.6f  %10.6f  %10.6f    %10.6f\n", n + 1, m + 1, fed_matrix[0](n, m),
               fed_matrix[1](n, m), delta_da_nm, coupling);
      }
    }
    _fed_matrix = fed_matrix;
    printf(" --------------------------------------------------------------------- \n");
  }
  // Performs Fragment Charge Difference (FCD) approach [J. Chem. Phys., 2002, 117, 5607]
  if (settings.FCD) {
    auto fcd_matrix = calculateFXDMatrix(_excVecs, false);
    printf("\n --------------------- FCD Coupling calculation ---------------------- \n");
    printf("     States      Q(D)        Q(A)         dQ12      Coupling(eV)   \n");
    printf(" --------------------------------------------------------------------- \n");

    for (unsigned int n = 0; n < (*_excVecs).cols(); n++) {
      double delta_da_nn = fcd_matrix[0](n, n) - fcd_matrix[1](n, n);
      for (unsigned int m = n; m < (*_excVecs).cols(); m++) {
        double delta_da_nm = fcd_matrix[0](n, m) - fcd_matrix[1](n, m);
        double delta_da_mm = fcd_matrix[0](m, m) - fcd_matrix[1](m, m);
        double numerator = ((*_excEner)(n) - (*_excEner)(m)) * abs(delta_da_nm);
        double coupling = 0.0;
        if (numerator != 0.0) {
          double denominator = sqrt(pow(delta_da_nn - delta_da_mm, 2) + 4.0 * pow(delta_da_nm, 2));
          coupling = (numerator / denominator) * HARTREE_TO_EV;
        }
        printf("  %3i  %3i   %10.6f  %10.6f  %10.6f    %10.6f\n", n + 1, m + 1, fcd_matrix[0](n, m),
               fcd_matrix[1](n, m), delta_da_nm, coupling);
      }
    }
    _fcd_matrix = fcd_matrix;
    printf(" --------------------------------------------------------------------- \n");
  }
  // Performs multistate FED-FCD approach [Photosynth. Res., 2018, 137, 215]
  if (settings.multistateFXD) {
    auto fed_matrix = calculateFXDMatrix(_excVecs, true);
    auto fcd_matrix = calculateFXDMatrix(_excVecs, false);
    // Calculate dq-dx matrix to separate CT and LE excitations
    Eigen::MatrixXd dx = 0.5 * (fed_matrix[0] - fed_matrix[1]);
    Eigen::MatrixXd dq = 0.5 * (fcd_matrix[0] - fcd_matrix[1]);
    Eigen::MatrixXd dq_dx = dq * dq - dx * dx;

    // Diagonalization of the matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_dq_dx(dq_dx);
    Eigen::MatrixXd u_1 = eigenSolver_dq_dx.eigenvectors();
    Eigen::VectorXd u_1_eigenvalues = eigenSolver_dq_dx.eigenvalues();
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(u_1.rows(), u_1.cols());
    for (unsigned int i = 0; i < u_1.cols(); i++) {
      temp.col(i) = u_1.col(u_1.cols() - i - 1);
    }
    u_1 = temp;
    temp.resize(0, 0);
    // Separate local and CT states
    unsigned int CT_counter = 0;
    unsigned int LE_counter = 0;
    for (unsigned int nEig = 0; nEig < u_1_eigenvalues.size(); nEig++) {
      if (u_1_eigenvalues[nEig] > 0.0) {
        CT_counter += 1;
        // This one is smaller equal 0.0 for no specific reason
      }
      else if (u_1_eigenvalues[nEig] <= 0.0) {
        LE_counter += 1;
      }
    }
    Eigen::MatrixXd transformed = u_1.transpose() * (*_excEner).asDiagonal() * u_1;

    dx = u_1.transpose() * dx * u_1;
    dq = u_1.transpose() * dq * u_1;

    // Second transformation for LE /CT subblocks
    // Large transformation matrix
    Eigen::MatrixXd u_2 = Eigen::MatrixXd::Zero(u_1.rows(), u_1.cols());
    // Subtransformation matrices for LE/CT space
    Eigen::MatrixXd u_2_LE;
    Eigen::VectorXd u_2_LE_eigenvalues;
    Eigen::MatrixXd u_2_CT;
    Eigen::VectorXd u_2_CT_eigenvalues;

    if (CT_counter != 0) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_dq(dq.block(0, 0, CT_counter, CT_counter));
      u_2_CT = eigenSolver_dq.eigenvectors();
      u_2_CT_eigenvalues = eigenSolver_dq.eigenvalues();
      u_2.block(0, 0, CT_counter, CT_counter) = u_2_CT;
    }

    if (LE_counter != 0) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_dx(
          dx.block(CT_counter, CT_counter, LE_counter, LE_counter).eval());
      u_2_LE = eigenSolver_dx.eigenvectors();
      u_2_LE_eigenvalues = eigenSolver_dx.eigenvalues();
      if (CT_counter != 0) {
        u_2.block(CT_counter, CT_counter, LE_counter, LE_counter) = u_2_LE;
      }
      else {
        u_2 = u_2_LE;
      }
    }

    // Second transformation
    transformed = (u_2.transpose() * transformed * u_2).eval();
    // third transformation to diagonalize subblock (LE/CT) for donor/acceptor subblock
    Eigen::MatrixXd u_3 = Eigen::MatrixXd::Zero(u_1.rows(), u_1.cols());
    // CT block
    unsigned int CT1_counter = 0;
    unsigned int CT2_counter = 0;
    if (CT_counter != 0) {
      for (unsigned int nEig = 0; nEig < u_2_CT_eigenvalues.size(); nEig++) {
        if (u_2_CT_eigenvalues[nEig] <= 0.0) {
          CT1_counter += 1;
        }
        else {
          CT2_counter += 1;
        }
      }
      if (CT1_counter != 0) {
        Eigen::MatrixXd CT1_block = transformed.block(0, 0, CT1_counter, CT1_counter);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_CT1(CT1_block);
        Eigen::MatrixXd u_3_CT_1 = eigenSolver_CT1.eigenvectors();
        u_3.block(0, 0, u_3_CT_1.rows(), u_3_CT_1.cols()) = u_3_CT_1;
      }
      if (CT2_counter != 0) {
        Eigen::MatrixXd CT2_block = transformed.block(CT1_counter, CT1_counter, CT2_counter, CT2_counter);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_CT2(CT2_block);
        Eigen::MatrixXd u_3_CT_2 = eigenSolver_CT2.eigenvectors();
        u_3.block(CT1_counter, CT1_counter, u_3_CT_2.rows(), u_3_CT_2.cols()) = u_3_CT_2;
      }
    }
    // LE block
    unsigned int combined_CT_block = CT1_counter + CT2_counter;
    unsigned int LE1_counter = 0;
    unsigned int LE2_counter = 0;
    if (LE_counter != 0) {
      for (unsigned int nEig = 0; nEig < u_2_LE_eigenvalues.size(); nEig++) {
        if (u_2_LE_eigenvalues[nEig] <= 0.0) {
          LE1_counter += 1;
        }
        else {
          LE2_counter += 1;
        }
      }
      if (LE1_counter != 0) {
        Eigen::MatrixXd LE1_block = transformed.block(combined_CT_block, combined_CT_block, LE1_counter, LE1_counter);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_LE1(LE1_block);
        Eigen::MatrixXd u_3_LE_1 = eigenSolver_LE1.eigenvectors();
        u_3.block(combined_CT_block, combined_CT_block, u_3_LE_1.rows(), u_3_LE_1.cols()) = u_3_LE_1;
      }
      if (LE2_counter != 0) {
        Eigen::MatrixXd LE2_block =
            transformed.block(combined_CT_block + LE1_counter, combined_CT_block + LE1_counter, LE2_counter, LE2_counter);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver_LE2(LE2_block);
        Eigen::MatrixXd u_3_LE_2 = eigenSolver_LE2.eigenvectors();
        u_3.block(combined_CT_block + LE1_counter, combined_CT_block + LE1_counter, u_3_LE_2.rows(), u_3_LE_2.cols()) = u_3_LE_2;
      }
    }
    // Third transformation
    transformed = u_3.transpose() * transformed * u_3;
    // Save transformed matrix and save transformed eigenvectors if requested
    _fxd_matrix = transformed;
    if (settings.writeTransformedExcitationVectors) {
      (*_excVecs) = (*_excVecs) * u_1 * u_2 * u_3;
      (*_excEner) = transformed.diagonal();
      std::cout << " Write the transformed excitation vectors and diabatic excitation energies to disk:" << std::endl;
      std::shared_ptr<std::vector<Eigen::MatrixXd>> excvec =
          std::make_shared<std::vector<Eigen::MatrixXd>>(std::vector<Eigen::MatrixXd>(1, *_excVecs));
      _lrscf->setSolution(excvec, _excEner, settings.loadType);
    }
    printf("\n ----------------------------------------------------------------\n");
    printf("                   Multistate-FED-FCD Calculation                 \n");
    printf(" ---------------------------------------------------------------- \n");
    // CT1 diabatic excitation energies and couplings
    if (CT1_counter != 0) {
      printf("\n ------------------ CT(A->D) --------------------- \n");
      printf("     State        energy(eV)       \n");
      printf(" --------------------------------------------------\n");
      for (unsigned int i = 0; i < CT1_counter; i++) {
        printf("     %3i        %10.6f \n", i + 1, transformed(i, i) * HARTREE_TO_EV);
      }
      printf(" --------------------------------------------\n");
    }
    if (CT2_counter != 0) {
      printf("\n ------------------ CT(D->A) --------------------- \n");
      printf("     State        energy(eV)       \n");
      printf(" -------------------------------------------------\n");
      for (unsigned int i = 0; i < CT2_counter; i++) {
        printf("     %3i        %10.6f \n", i + 1, transformed(CT1_counter + i, CT1_counter + i) * HARTREE_TO_EV);
      }
      printf(" -------------------------------------------------\n");
    }
    if (LE1_counter != 0) {
      printf("\n ------------------ LE(A*) --------------------- \n");
      printf("     State        energy(eV)       \n");
      printf(" -----------------------------------------------\n");
      for (unsigned int i = 0; i < LE1_counter; i++) {
        printf("     %3i        %10.6f \n", i + 1, transformed(combined_CT_block + i, combined_CT_block + i) * HARTREE_TO_EV);
      }
      printf(" --------------------------------------------\n");
    }
    if (LE2_counter != 0) {
      printf("\n ------------------ LE(D*) --------------------- \n");
      printf("     State        energy(eV)       \n");
      printf(" -----------------------------------------------\n");
      for (unsigned int i = 0; i < LE2_counter; i++) {
        printf("     %3i        %10.6f \n", i + 1,
               transformed(combined_CT_block + LE1_counter + i, combined_CT_block + LE1_counter + i) * HARTREE_TO_EV);
      }
      printf(" --------------------------------------------\n");
    }
    if (CT2_counter != 0 && CT1_counter != 0) {
      printf("\n ------------------ CT(A->D) x CT(D->A) coupling (eV) --------------------- \n");
      std::cout << transformed.block(0, CT1_counter, CT1_counter, CT2_counter).array().abs().matrix() * HARTREE_TO_EV
                << std::endl;
      printf(" --------------------------------------------------------------------------\n");
    }
    if (LE1_counter != 0 && LE2_counter != 0) {
      printf("\n ------------------ LE(A*) x LE(D*) coupling (eV) --------------------- \n");
      std::cout << transformed.block(combined_CT_block, combined_CT_block + LE1_counter, LE1_counter, LE2_counter)
                           .array()
                           .abs()
                           .matrix() *
                       HARTREE_TO_EV
                << std::endl;
      printf(" ----------------------------------------------------------------------\n");
    }
    if (LE1_counter != 0 && CT1_counter != 0) {
      printf("\n ------------------ CT(A->D) x LE(A*) coupling (eV) --------------------- \n");
      std::cout << transformed.block(0, combined_CT_block, CT1_counter, LE1_counter).array().abs().matrix() * HARTREE_TO_EV
                << std::endl;
      printf(" ------------------------------------------------------------------------\n");
    }
    if (LE1_counter != 0 && CT2_counter != 0) {
      printf("\n ------------------ CT(D->A) x LE(A*) coupling (eV) --------------------- \n");
      std::cout << transformed.block(CT1_counter, combined_CT_block, CT2_counter, LE1_counter).array().abs().matrix() * HARTREE_TO_EV
                << std::endl;
      printf(" ------------------------------------------------------------------------\n");
    }
    if (LE2_counter != 0 && CT1_counter != 0) {
      printf("\n ------------------ CT(A->D) x LE(D*) coupling (eV) --------------------- \n");
      std::cout << transformed.block(0, combined_CT_block + LE1_counter, CT1_counter, LE2_counter).array().abs().matrix() * HARTREE_TO_EV
                << std::endl;
      printf(" ------------------------------------------------------------------------\n");
    }
    if (LE2_counter != 0 && CT2_counter != 0) {
      printf("\n ------------------ CT(D->A) x LE(D*) coupling (eV) --------------------- \n");
      std::cout
          << transformed.block(CT1_counter, combined_CT_block + LE1_counter, CT2_counter, LE2_counter).array().abs().matrix() * HARTREE_TO_EV
          << std::endl;
      printf(" ------------------------------------------------------------------------\n");
    }
  }
  return;
}

template<Options::SCF_MODES SCFMode>
std::vector<Eigen::MatrixXd> FXDTask<SCFMode>::calculateFXDMatrix(std::shared_ptr<Eigen::MatrixXd> excVector, bool FED) {
  // Initialize some variables
  CoefficientMatrix<SCFMode> coefs =
      _systemController->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->getCoefficients();
  SpinPolarizedData<SCFMode, Eigen::VectorXd> orbitalEnergies =
      _systemController->getElectronicStructure<SCFMode>()->getMolecularOrbitals()->getEigenvalues();
  SpinPolarizedData<SCFMode, unsigned int> nOcc = _systemController->getNOccupiedOrbitals<SCFMode>();
  SpinPolarizedData<SCFMode, unsigned int> nVirt = _systemController->getNVirtualOrbitals<SCFMode>();

  Eigen::MatrixXd overlapMatrix = _systemController->getOneElectronIntegralController()->getOverlapIntegrals();
  unsigned int nBasisFunc = overlapMatrix.rows();
  auto densMat = _systemController->getElectronicStructure<SCFMode>()->getDensityMatrix();
  // Basis indices associated with donor and acceptor subsystems
  if (settings.donoratoms.size() > 2 || settings.acceptoratoms.size() > 2)
    throw SerenityError("More than two indices for donor/acceptor atoms given!");
  const auto basisIndices = _systemController->getAtomCenteredBasisController()->getBasisIndices();
  unsigned int indexStartDonor = basisIndices[settings.donoratoms[0]].first;
  unsigned int indexEndDonor = basisIndices[settings.donoratoms[1]].second;
  unsigned int indexStartAcceptor = basisIndices[settings.acceptoratoms[0]].first;
  unsigned int indexEndAcceptor = basisIndices[settings.acceptoratoms[1]].second;

  Eigen::MatrixXd fxd_values_d = Eigen::MatrixXd::Zero((*excVector).cols(), (*excVector).cols());
  Eigen::MatrixXd fxd_values_a = Eigen::MatrixXd::Zero((*excVector).cols(), (*excVector).cols());

  for (unsigned int n = 0; n < (*excVector).cols(); n++) {
    unsigned int iStartSpin = 0;

    for_spin(nOcc, nVirt, coefs) {
      Eigen::MatrixXd excMatrix_n(nOcc_spin, nVirt_spin);
      for (unsigned int i = 0, ia = iStartSpin; i < nOcc_spin; ++i) {
        for (unsigned int a = 0; a < nVirt_spin; ++a, ++ia) {
          excMatrix_n(i, a) = (*excVector)(ia, n);
        }
      }
      // First half transformation excitation vector matrix from AO to MO
      Eigen::MatrixXd transformed_Occ_n = coefs_spin.block(0, 0, nBasisFunc, nOcc_spin) * excMatrix_n;
      Eigen::MatrixXd transformed_Virt_n = excMatrix_n * coefs_spin.block(0, nOcc_spin, nBasisFunc, nVirt_spin).transpose();

      Eigen::MatrixXd overlap_sqrt;
      if (settings.loewdinpopulation)
        overlap_sqrt = mSqrt_Sym(overlapMatrix);

      for (unsigned int m = 0; m < (*excVector).cols(); m++) {
        Eigen::MatrixXd excMatrix_m(nOcc_spin, nVirt_spin);
        for (unsigned int i = 0, ia = iStartSpin; i < nOcc_spin; ++i) {
          for (unsigned int a = 0; a < nVirt_spin; ++a, ++ia) {
            excMatrix_m(i, a) = (*excVector)(ia, m);
          }
        }
        // Second half transformation excitation vector matrix from AO to MO
        Eigen::MatrixXd transformed_Occ_m = coefs_spin.block(0, 0, nBasisFunc, nOcc_spin) * excMatrix_m; // Full
                                                                                                         // transformation
                                                                                                         // in MO basis
        transformed_Occ_m = transformed_Occ_n * transformed_Occ_m.transpose();
        Eigen::MatrixXd transformed_Virt_m = excMatrix_m * coefs_spin.block(0, nOcc_spin, nBasisFunc, nVirt_spin).transpose();
        transformed_Virt_m = transformed_Virt_n.transpose() * transformed_Virt_m;
        Eigen::MatrixXd p_times_s;
        if (FED) {
          if (!settings.loewdinpopulation) {
            p_times_s = (transformed_Occ_m + transformed_Virt_m).transpose() * overlapMatrix;
          }
          else {
            p_times_s = overlap_sqrt * (transformed_Occ_m + transformed_Virt_m).transpose() * overlap_sqrt;
          }
        }
        else {
          if (!settings.loewdinpopulation) {
            p_times_s = (transformed_Occ_m - transformed_Virt_m).transpose() * overlapMatrix;
          }
          else {
            p_times_s = overlap_sqrt * (transformed_Occ_m - transformed_Virt_m).transpose() * overlap_sqrt;
          }
        }

        fxd_values_a(n, m) = p_times_s
                                 .block(indexStartAcceptor, indexStartAcceptor, indexEndAcceptor - indexStartAcceptor,
                                        indexEndAcceptor - indexStartAcceptor)
                                 .trace();
        fxd_values_d(n, m) =
            p_times_s
                .block(indexStartDonor, indexStartDonor, indexEndDonor - indexStartDonor, indexEndDonor - indexStartDonor)
                .trace();
      }
      iStartSpin += nOcc_spin * nVirt_spin;
    };
  }
  std::vector<Eigen::MatrixXd> fxd_matrix(2);
  // symmetrization
  fxd_matrix[0] = 0.5 * (fxd_values_d + fxd_values_d.transpose());
  fxd_matrix[1] = 0.5 * (fxd_values_a + fxd_values_a.transpose());
  return fxd_matrix;
}

template class FXDTask<Options::SCF_MODES::RESTRICTED>;
template class FXDTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */