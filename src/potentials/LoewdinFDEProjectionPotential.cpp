/**
 * @file LoewdinFDEProjectionPotential.cpp
 *
 * @date Apr 22, 2024
 * @author Denis G. Artiukhin
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
#include "potentials/LoewdinFDEProjectionPotential.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "misc/SystemSplittingTools.h"
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundleFactory.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LoewdinFDEProjectionPotential<SCFMode>::LoewdinFDEProjectionPotential(std::shared_ptr<SystemController> activeSystem,
                                                                      std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                                      const EmbeddingSettings& settings)
  : Potential<SCFMode>(activeSystem->getBasisController()),
    _activeSystem(activeSystem),
    _loewdinOrder(settings.loewdinOrder),
    _loewdinWeights(settings.loewdinWeights) {
  // Sanity checks
  if (settings.longRangeNaddKinFunc != CompositeFunctionals::KINFUNCTIONALS::NONE)
    throw SerenityError("Loewdin embedding: Hybrid algorithm is not yet implemented!");
  if (environmentSystems.size() == 0)
    throw SerenityError("Loewdin embedding: There is no environment system given!");
  if (environmentSystems.size() > 1)
    throw SerenityError("Loewdin embedding: Implementation for multiple environment subsystems is not yet available!");
  if (settings.truncateProjector)
    throw SerenityError("Loewdin embedding: Projector truncation options are not yet implemented!");

  for (auto e : environmentSystems)
    _environmentSystems.push_back(e);

  // Expects that an electronic structure is already in place for the active system!
  activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->addSensitiveObject(
      ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  _envDensityCont = SystemSplittingTools<SCFMode>::getEnvironmentDensityControllers(environmentSystems, false);

  // Get the off-diagonal overlap matrix S_AB
  _overlapABs = SystemSplittingTools<SCFMode>::getProjectedSubsystems(activeSystem, environmentSystems, -10.0);
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& LoewdinFDEProjectionPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    takeTime("Evaluating Loewdin correction");

    // reset the fock matrix
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin.setZero();
    };

    switch (_loewdinOrder) {
      case 0: {
        // Nothing to be done here
        break;
      }
      case 1: {
        this->computeFirstOrderTruncSeries();
        break;
      }
      case 2: {
        this->computeSecondOrderTruncSeries();
        break;
      }
      case 3: {
        this->computeThirdOrderTruncSeries();
        break;
      }
      default: {
        throw SerenityError("Requested Loewdin expansion order is not implemented");
      }
    }

    timeTaken(2, "Evaluating Loewdin correction");
  } // if (!_potential)

  return *_potential;
}

template<Options::SCF_MODES SCFMode>
void LoewdinFDEProjectionPotential<SCFMode>::computeFirstOrderTruncSeries() {
  if (_overlapABs[0]) {
    auto environmentSystem = _environmentSystems[0].lock();
    auto activeSystem = _activeSystem.lock();

    // Get P_BB
    DensityMatrix<SCFMode> P_BB = _envDensityCont[0]->getDensityMatrix();

    // K_AB = S_AB * P_BB
    const Eigen::MatrixXd& S_AB = *_overlapABs[0];
    for_spin(_kAB, P_BB) {
      _kAB_spin = S_AB * P_BB_spin;
    };

    // Construct off-diagonal kinetic energy matrix T_AB
    auto& libint = Libint::getInstance();
    Eigen::MatrixXd T_AB = libint
                               .compute1eInts(LIBINT_OPERATOR::kinetic, activeSystem->getBasisController(),
                                              environmentSystem->getBasisController())
                               .transpose();

    // L_AA = K_AB * T_BA
    for_spin(_lAA, _kAB) {
      _lAA_spin = _kAB_spin * T_AB.transpose();
    };

    // A factor of one half for restricted to account for the factor of two in the density matrix.
    double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;

    // Compute the potential
    auto& pot = *_potential;
    for_spin(pot, _lAA) {
      pot_spin -= scfFactor * _loewdinWeights[0] * (_lAA_spin + _lAA_spin.transpose());
    };
  } // if (s_ABs[0])
  else {
    throw SerenityError("Off-diagonal overlaps are not available!");
  }
}

template<Options::SCF_MODES SCFMode>
void LoewdinFDEProjectionPotential<SCFMode>::computeSecondOrderTruncSeries() {
  // Compute lower-order terms and save intermediate results
  this->computeFirstOrderTruncSeries();

  auto environmentSystem = _environmentSystems[0].lock();
  auto activeSystem = _activeSystem.lock();

  // Construct diagonal kinetic energy matrices T_AA and T_BB
  auto& libint = Libint::getInstance();
  Eigen::MatrixXd t_AA = libint.compute1eInts(LIBINT_OPERATOR::kinetic, activeSystem->getBasisController());
  Eigen::MatrixXd t_BB = libint.compute1eInts(LIBINT_OPERATOR::kinetic, environmentSystem->getBasisController());

  // Get P_AA
  DensityMatrix<SCFMode> P_AA = activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrix();

  // M_AB = P_AA * S_AB
  const Eigen::MatrixXd& S_AB = *_overlapABs[0];
  for_spin(_mAB, P_AA) {
    _mAB_spin = P_AA_spin * S_AB;
  };

  // N_AA = K_AB * M_BA
  for_spin(_nAA, _kAB, _mAB) {
    _nAA_spin = _kAB_spin * _mAB_spin.transpose();
  };

  // A factor of one half for restricted to account for the factor of two in the density matrix.
  double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;

  // Compute the potential
  auto& pot = *_potential;
  for_spin(pot, _nAA, _kAB, _t2Corr) {
    auto mat = _nAA_spin * t_AA;
    _t2Corr_spin = scfFactor * scfFactor * 0.5 * _loewdinWeights[1] * (mat + mat.transpose());

    pot_spin += 2.0 * _t2Corr_spin;
    pot_spin += scfFactor * scfFactor * _loewdinWeights[1] * _kAB_spin * t_BB * _kAB_spin.transpose();
  };
}

template<Options::SCF_MODES SCFMode>
void LoewdinFDEProjectionPotential<SCFMode>::computeThirdOrderTruncSeries() {
  // Compute lower-order terms and save intermediate results
  this->computeSecondOrderTruncSeries();

  // A factor of one half for restricted to account for the factor of two in the density matrix.
  double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 0.5 : 1.0;

  // Compute the potential
  auto& pot = *_potential;
  for_spin(pot, _t3Corr, _lAA, _nAA) {
    Eigen::MatrixXd mat1 = _nAA_spin * _lAA_spin;
    Eigen::MatrixXd mat2 = _nAA_spin * _lAA_spin.transpose();

    _t3Corr_spin = -0.5 * scfFactor * scfFactor * scfFactor * _loewdinWeights[2] *
                   (mat1 + mat1.transpose() + mat2 + mat2.transpose());
    pot_spin += 2.0 * _t3Corr_spin;
  };
}

template<Options::SCF_MODES SCFMode>
double LoewdinFDEProjectionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  double energy = 0.0;
  double energy_corr = 0.0;
  if (!_potential)
    getMatrix();

  auto& pot = *_potential;
  for_spin(pot, P) {
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };

  if (_loewdinOrder > 1) {
    for_spin(_t2Corr, P) {
      energy_corr += _t2Corr_spin.cwiseProduct(P_spin).sum();
    };
  }

  if (_loewdinOrder > 2) {
    for_spin(_t3Corr, P) {
      energy_corr += _t3Corr_spin.cwiseProduct(P_spin).sum();
    };
  }

  energy -= energy_corr;
  return energy;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd LoewdinFDEProjectionPotential<SCFMode>::getGeomGradients() {
  Eigen::MatrixXd gradientContr(1, 3);
  gradientContr.setZero();
  return gradientContr;
}

template class LoewdinFDEProjectionPotential<Options::SCF_MODES::RESTRICTED>;
template class LoewdinFDEProjectionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
