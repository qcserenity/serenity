/**
 * @file HCorePotential.cpp
 *
 * @date Aug 16, 2017
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

/* Include Class Header*/
#include "potentials/HCorePotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/gradients/CoreCoreRepulsionDerivative.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libecpint.h"
#include "integrals/wrappers/Libint.h"
#include "misc/Timing.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HCorePotential<SCFMode>::HCorePotential(std::shared_ptr<SystemController> system)
  : Potential<SCFMode>(system->getBasisController()), _system(system), _potential(nullptr) {
  this->_basis->addSensitiveObject(_self);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& HCorePotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -     1e-Int Pot.");
  if (!_potential) {
    auto system = _system.lock();
    auto intC = system->getOneElectronIntegralController();
    _potential.reset(new FockMatrix<SCFMode>(intC->getOneElectronIntegrals()));
  }
  Timings::timeTaken("Active System -     1e-Int Pot.");
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
double HCorePotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (!_potential)
    this->getMatrix();
  Timings::takeTime("Active System -     1e-Int Pot.");
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -     1e-Int Pot.");
  return energy;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd HCorePotential<SCFMode>::getGeomGradients() {
  auto system = _system.lock();
  const auto& orbitalSet = system->template getActiveOrbitalController<SCFMode>();
  auto atoms = system->getAtoms();
  unsigned int nAtoms = atoms.size();
  unsigned int maxAtoms = 2;
  Eigen::MatrixXd gradientContr(nAtoms, 3);
  gradientContr.setZero();

  DensityMatrix<SCFMode> matrix(system->template getElectronicStructure<SCFMode>()->getDensityMatrix());
  matrix = calcEnergyWeightedDensityMatrix(system, orbitalSet);
  for_spin(matrix) {
    matrix_spin *= -1.0;
  };

  unsigned int nBasisFunctionsRed = orbitalSet->getBasisController()->getReducedNBasisFunctions();
  auto basisIndicesRed = system->getAtomCenteredBasisController()->getBasisIndicesRed();
  auto mapping = createBasisToAtomIndexMapping(basisIndicesRed, nBasisFunctionsRed);

  auto& basis = system->getAtomCenteredBasisController()->getBasis();
  Libint& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::overlap, 1, 2);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    bool significant;
    Eigen::MatrixXd gradientContrPriv(nAtoms, 3);
    gradientContrPriv.setZero();
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < basis.size(); ++i) {
      unsigned int offI = system->getBasisController()->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; ++j) {
        unsigned int offJ = system->getBasisController()->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        significant = libint.compute(libint2::Operator::overlap, 1, *basis[i], *basis[j], intDerivs);
        if (significant) {
          double perm = (i == j ? 1.0 : 2.0);
          for (unsigned int iAtom = 0; iAtom < maxAtoms; ++iAtom) {
            unsigned int nAtom;
            switch (iAtom) {
              case (0):
                nAtom = mapping[i];
                break;
              case (1):
                nAtom = mapping[j];
                break;
              default:
                nAtom = iAtom - 2;
                break;
            }
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom * 3 + iDirection).data(), nJ, nI);
              for_spin(matrix) {
                gradientContrPriv(nAtom, iDirection) += perm * tmp.cwiseProduct(matrix_spin.block(offJ, offI, nJ, nI)).sum();
              };
            }
          }
        }
      }
    }

#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */

  DensityMatrix<RESTRICTED> densMatrix(system->template getElectronicStructure<SCFMode>()->getDensityMatrix().total());

  libint.initialize(libint2::Operator::nuclear, 1, 2, atoms);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    bool significant;
    Eigen::MatrixXd gradientContrPriv(nAtoms, 3);
    gradientContrPriv.setZero();
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < basis.size(); ++i) {
      unsigned int offI = system->getBasisController()->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; ++j) {
        unsigned int offJ = system->getBasisController()->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        significant = libint.compute(libint2::Operator::nuclear, 1, *basis[i], *basis[j], intDerivs);
        maxAtoms = nAtoms + 2;
        if (significant) {
          double perm = (i == j ? 1.0 : 2.0);
          for (unsigned int iAtom = 0; iAtom < maxAtoms; ++iAtom) {
            unsigned int nAtom;
            switch (iAtom) {
              case (0):
                nAtom = mapping[i];
                break;
              case (1):
                nAtom = mapping[j];
                break;
              default:
                nAtom = iAtom - 2;
                break;
            }
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom * 3 + iDirection).data(), nJ, nI);

              gradientContrPriv(nAtom, iDirection) += perm * tmp.cwiseProduct(densMatrix.block(offJ, offI, nJ, nI)).sum();
            }
          }
        }
      }
    }

#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */

  libint.initialize(libint2::Operator::kinetic, 1, 2);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    bool significant;
    Eigen::MatrixXd gradientContrPriv(nAtoms, 3);
    gradientContrPriv.setZero();
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < basis.size(); ++i) {
      unsigned int offI = system->getBasisController()->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; ++j) {
        unsigned int offJ = system->getBasisController()->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        significant = libint.compute(libint2::Operator::kinetic, 1, *basis[i], *basis[j], intDerivs);
        maxAtoms = 2;
        if (significant) {
          double perm = (i == j ? 1.0 : 2.0);
          for (unsigned int iAtom = 0; iAtom < maxAtoms; ++iAtom) {
            unsigned int nAtom;
            switch (iAtom) {
              case (0):
                nAtom = mapping[i];
                break;
              case (1):
                nAtom = mapping[j];
                break;
              default:
                nAtom = iAtom - 2;
                break;
            }
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom * 3 + iDirection).data(), nJ, nI);
              gradientContrPriv(nAtom, iDirection) += perm * tmp.cwiseProduct(densMatrix.block(offJ, offI, nJ, nI)).sum();
            }
          }
        }
      }
    }

#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */
  libint.finalize(libint2::Operator::kinetic, 1, 2);
  libint.finalize(libint2::Operator::nuclear, 1, 2);
  libint.finalize(libint2::Operator::overlap, 1, 2);

  // Adding CC potential

  gradientContr += CoreCoreRepulsionDerivative::calculateDerivative(atoms);

  return gradientContr;
}

template<Options::SCF_MODES SCFMode>
std::vector<unsigned int>
HCorePotential<SCFMode>::createBasisToAtomIndexMapping(const std::vector<std::pair<unsigned int, unsigned int>>& basisIndicesRed,
                                                       unsigned int nBasisFunctionsRed) {
  std::vector<unsigned int> mapping(nBasisFunctionsRed);
  // Vector to check whether ALL basis function shells are assigned to an atom index
  std::vector<bool> hasElementBeenSet(nBasisFunctionsRed, false);
  for (unsigned int iAtom = 0; iAtom < basisIndicesRed.size(); ++iAtom) {
    const unsigned int firstIndex = basisIndicesRed[iAtom].first;
    const unsigned int endIndex = basisIndicesRed[iAtom].second;
    for (unsigned int iShell = firstIndex; iShell < endIndex; ++iShell) {
      mapping[iShell] = iAtom;
      hasElementBeenSet[iShell] = true;
    }
  }
  // Check
  for (bool x : hasElementBeenSet) {
    if (not x)
      throw SerenityError("HCorePotential: Missed gradient element in gradient evaluation.");
  }
  return mapping;
}

template<Options::SCF_MODES SCFMode>
DensityMatrix<SCFMode>
HCorePotential<SCFMode>::calcEnergyWeightedDensityMatrix(std::shared_ptr<SystemController> systemController,
                                                         const std::shared_ptr<OrbitalController<SCFMode>>& orbitalSet) {
  DensityMatrix<SCFMode> energyWeightedDensityMatrix(orbitalSet->getBasisController());
  const unsigned int nBasisFunctions = systemController->getBasisController()->getNBasisFunctions();
  const auto& nElectrons = systemController->getNElectrons<SCFMode>();
  auto eigenvalues(std::move(orbitalSet->getEigenvalues()));
  auto coefficients(std::move(orbitalSet->getCoefficients()));
  const double occ = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
  for_spin(nElectrons, energyWeightedDensityMatrix, eigenvalues, coefficients) {
    for (unsigned int i = 0; i < nBasisFunctions; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        // TODO this is wrong for a non-Aufbau occupation
        for (unsigned int k = 0; k < nElectrons_spin / occ; ++k) {
          energyWeightedDensityMatrix_spin(i, j) +=
              occ * eigenvalues_spin[k] * coefficients_spin(i, k) * coefficients_spin(j, k);
        }
        energyWeightedDensityMatrix_spin(j, i) = energyWeightedDensityMatrix_spin(i, j);
      }
    }
  };
  return energyWeightedDensityMatrix;
}

template class HCorePotential<Options::SCF_MODES::RESTRICTED>;
template class HCorePotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
