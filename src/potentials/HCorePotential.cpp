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
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Atom.h"
#include "geometry/efields/EFieldPlates.h"
#include "geometry/gradients/CoreCoreRepulsionDerivative.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libecpint.h"
#include "integrals/wrappers/Libint.h"
#include "misc/Timing.h"
#include "settings/Settings.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
HCorePotential<SCFMode>::HCorePotential(std::shared_ptr<SystemController> system)
  : Potential<SCFMode>(system->getBasisController()), _system(system), _potential(nullptr), _extCharges(0) {
  this->_basis->addSensitiveObject(_self);
  auto ef = _system.lock()->getSettings().efield;
  if (ef.use && !ef.analytical) {
    EFieldPlates plates(Eigen::Map<Eigen::Vector3d>(&ef.pos1[0]), Eigen::Map<Eigen::Vector3d>(&ef.pos2[0]), ef.distance,
                        ef.nRings, ef.radius, ef.fieldStrength, ef.nameOutput);
    _extCharges.insert(_extCharges.end(), plates.getPairList().begin(), plates.getPairList().end());
  }
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& HCorePotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -     1e-Int Pot.");
  if (!_potential) {
    auto intC = _system.lock()->getOneElectronIntegralController();
    _potential = std::make_unique<FockMatrix<SCFMode>>(intC->getOneElectronIntegrals());
    auto ef = _system.lock()->getSettings().efield;

    if (ef.use && ef.analytical) {
      // creates an HCore potential with dipole matrix and EEF.
      auto dipoleL = intC->getDipoleLengths();
      Eigen::Vector3d fieldVec(ef.pos2[0] - ef.pos1[0], ef.pos2[1] - ef.pos1[1], ef.pos2[2] - ef.pos1[2]);
      fieldVec.normalize();
      printSmallCaption("Homogeneous Electric Field");
      printf("  Field Strength (au): %4.2e\n\n", ef.fieldStrength);
      printf("  Electric Field Vector: %7.3f %7.3f %7.3f\n\n", fieldVec(0), fieldVec(1), fieldVec(2));
      auto& pot = *_potential;
      for_spin(pot) {
        for (unsigned int j = 0; j < fieldVec.size(); ++j) {
          pot_spin -= ef.fieldStrength * fieldVec(j) * dipoleL[j];
        }
      };
    }
    if (!_extCharges.empty()) {
      auto& libint = Libint::getInstance();
      auto& pot = *_potential;
      Eigen::MatrixXd ext = libint.compute1eInts(LIBINT_OPERATOR::nuclear, this->_basis, _extCharges);
      for_spin(pot) {
        pot_spin += ext;
      };
    }
  }
  Timings::timeTaken("Active System -     1e-Int Pot.");
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
double HCorePotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (!_potential) {
    this->getMatrix();
  }

  Timings::takeTime("Active System -     1e-Int Pot.");
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  auto ef = _system.lock()->getSettings().efield;
  if (ef.use && ef.analytical) {
    Eigen::Vector3d fieldVec(ef.pos2[0] - ef.pos1[0], ef.pos2[1] - ef.pos1[1], ef.pos2[2] - ef.pos1[2]);
    fieldVec.normalize();
    fieldVec = ef.fieldStrength * fieldVec;
    for (auto& atom : _system.lock()->getAtoms()) {
      energy -= atom->getEffectiveCharge() * fieldVec.dot(atom->coords());
    }
  }

  // add energy which is due to interaction with external charges
  for (auto& charge : _extCharges) {
    for (auto& atom : _system.lock()->getAtoms()) {
      energy += atom->getEffectiveCharge() * charge.first /
                (atom->coords() - Eigen::Map<Eigen::Vector3d>(&charge.second[0])).norm();
    }
  }
  Timings::timeTaken("Active System -     1e-Int Pot.");
  return energy;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd HCorePotential<SCFMode>::getGeomGradients() {
  auto system = _system.lock();
  // ToDo: Ask electronicstructure for new fock matrix?! orbital set could be localized and not canonical
  const auto& orbitalSet = system->template getActiveOrbitalController<SCFMode>();
  auto atoms = system->getAtoms();
  const unsigned int nAtoms = atoms.size();
  Eigen::MatrixXd gradientContr = Eigen::MatrixXd::Zero(nAtoms, 3);

  DensityMatrix<SCFMode> matrix(system->template getElectronicStructure<SCFMode>()->getDensityMatrix());
  matrix = calcEnergyWeightedDensityMatrix(system, orbitalSet);
  for_spin(matrix) {
    matrix_spin *= -1.0;
  };
  auto basisIndicesRed = system->getAtomCenteredBasisController()->getBasisIndicesRed();
  auto mapping = system->getAtomCenteredBasisController()->getAtomIndicesOfBasisShells();
  auto& basis = system->getAtomCenteredBasisController()->getBasis();
  Libint& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::overlap, 1, 2);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    bool significant;
    Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(nAtoms, 3);
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < basis.size(); ++i) {
      unsigned int offI = system->getBasisController()->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; ++j) {
        unsigned int offJ = system->getBasisController()->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        significant = libint.compute(LIBINT_OPERATOR::overlap, 1, *basis[i], *basis[j], intDerivs);
        if (significant) {
          double perm = (i == j ? 1.0 : 2.0);
          for (unsigned int iCol = 0; iCol < 2; ++iCol) {
            unsigned iAtom = iCol ? mapping[j] : mapping[i];
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iCol * 3 + iDirection).data(), nJ, nI);
              for_spin(matrix) {
                gradientContrPriv(iAtom, iDirection) += perm * tmp.cwiseProduct(matrix_spin.block(offJ, offI, nJ, nI)).sum();
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

  libint.initialize(LIBINT_OPERATOR::nuclear, 1, 2, atoms);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    bool significant;
    Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(nAtoms, 3);
#pragma omp for schedule(static, 1)
    for (unsigned int i = 0; i < basis.size(); ++i) {
      unsigned int offI = system->getBasisController()->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; ++j) {
        unsigned int offJ = system->getBasisController()->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        significant = libint.compute(LIBINT_OPERATOR::nuclear, 1, *basis[i], *basis[j], intDerivs);
        if (significant) {
          double perm = (i == j ? 1.0 : 2.0);
          // loop over integral centers for atom pairs, where iCol < 2
          for (unsigned int iCol = 0; iCol < 2; ++iCol) {
            unsigned iAtom = iCol ? mapping[j] : mapping[i];
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iCol * 3 + iDirection).data(), nJ, nI);
              gradientContrPriv(iAtom, iDirection) += perm * tmp.cwiseProduct(densMatrix.block(offJ, offI, nJ, nI)).sum();
            }
          }
          // loop over integral centers for atom pairs, where iCol > 2
          for (unsigned int iCol = 2; iCol < nAtoms + 2; ++iCol) {
            unsigned iAtom = iCol - 2;
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iCol * 3 + iDirection).data(), nJ, nI);
              gradientContrPriv(iAtom, iDirection) += perm * tmp.cwiseProduct(densMatrix.block(offJ, offI, nJ, nI)).sum();
            }
          }
        }
      }
    }

#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */

  libint.initialize(LIBINT_OPERATOR::kinetic, 1, 2);
#pragma omp parallel
  {
    Eigen::MatrixXd intDerivs;
    bool significant;
    Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(nAtoms, 3);
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < basis.size(); ++i) {
      unsigned int offI = system->getBasisController()->extendedIndex(i);
      const unsigned int nI = basis[i]->getNContracted();
      for (unsigned int j = 0; j <= i; ++j) {
        unsigned int offJ = system->getBasisController()->extendedIndex(j);
        const unsigned int nJ = basis[j]->getNContracted();
        significant = libint.compute(LIBINT_OPERATOR::kinetic, 1, *basis[i], *basis[j], intDerivs);
        if (significant) {
          double perm = (i == j ? 1.0 : 2.0);
          for (unsigned int iCol = 0; iCol < 2; ++iCol) {
            unsigned iAtom = iCol ? mapping[j] : mapping[i];
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iCol * 3 + iDirection).data(), nJ, nI);
              gradientContrPriv(iAtom, iDirection) += perm * tmp.cwiseProduct(densMatrix.block(offJ, offI, nJ, nI)).sum();
            }
          }
        }
      }
    }

#pragma omp critical
    { gradientContr += gradientContrPriv; }
  } /* END OpenMP parallel */

  /* if analytical external electric field is present */
  auto ef = system->getSettings().efield;
  if (ef.use && ef.analytical) {
    Eigen::VectorXd fieldVec = Eigen::VectorXd::Zero(3);
    for (unsigned i = 0; i < 3; ++i)
      fieldVec[i] = ef.pos2[i] - ef.pos1[i];
    fieldVec.normalize();
    fieldVec *= ef.fieldStrength;
    libint.initialize(LIBINT_OPERATOR::emultipole1, 1, 2, Point(0, 0, 0));
#pragma omp parallel
    {
      Eigen::MatrixXd intDerivs;
      Eigen::MatrixXd gradientContrPriv = Eigen::MatrixXd::Zero(nAtoms, 3);
/* electronic part */
#pragma omp for schedule(dynamic)
      for (unsigned i = 0; i < basis.size(); ++i) {
        unsigned offI = system->getBasisController()->extendedIndex(i);
        unsigned nI = basis[i]->getNContracted();
        for (unsigned j = 0; j <= i; ++j) {
          unsigned offJ = system->getBasisController()->extendedIndex(j);
          unsigned nJ = basis[j]->getNContracted();
          if (libint.compute(LIBINT_OPERATOR::emultipole1, 1, *basis[i], *basis[j], intDerivs)) {
            double perm = (i == j ? 1.0 : 2.0);
            for (unsigned iCol = 0; iCol < 2; ++iCol) {
              unsigned iAtom = iCol ? mapping[j] : mapping[i];
              for (unsigned iDirection = 0; iDirection < 3; ++iDirection) {
                Eigen::MatrixXd dipSum = Eigen::MatrixXd::Zero(nJ, nI);
                for (unsigned dipDirection = 0; dipDirection < 3; ++dipDirection) {
                  dipSum += fieldVec(dipDirection) *
                            Eigen::Map<Eigen::MatrixXd>(
                                intDerivs.col(iCol * 12 + 1 + iDirection * 4 + dipDirection).data(), nJ, nI);
                }
                gradientContrPriv(iAtom, iDirection) +=
                    perm * dipSum.cwiseProduct(densMatrix.block(offJ, offI, nJ, nI)).sum();
              }
            }
          }
        }
      }
/* nuclear part */
#pragma omp for schedule(dynamic)
      for (unsigned iAtom = 0; iAtom < nAtoms; ++iAtom) {
        for (unsigned iDirection = 0; iDirection < 3; ++iDirection) {
          gradientContrPriv(iAtom, iDirection) -= atoms[iAtom]->getEffectiveCharge() * fieldVec(iDirection);
        }
      }
#pragma omp critical
      { gradientContr += gradientContrPriv; }
    } /* END OpenMP parallel */
    libint.finalize(LIBINT_OPERATOR::emultipole1, 1, 2);
  }
  else if (ef.use && !ef.analytical) {
    throw SerenityError("Gradients for numerical electric field are not implmented, yet!");
  }
  libint.finalize(LIBINT_OPERATOR::kinetic, 1, 2);
  libint.finalize(LIBINT_OPERATOR::nuclear, 1, 2);
  libint.finalize(LIBINT_OPERATOR::overlap, 1, 2);

  // Add ECP contribution
  gradientContr += Libecpint::computeECPGradientContribution(system->getAtomCenteredBasisController(), atoms, densMatrix);

  // Adding CC potential
  gradientContr += CoreCoreRepulsionDerivative::calculateDerivative(atoms);

  return gradientContr;
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
