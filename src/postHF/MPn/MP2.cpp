/**
 * @file   MP2.cpp
 *
 * @date   Jul 14, 2014
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
#include "postHF/MPn/MP2.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/transformer/Ao2MoTransformer.h"
#include "io/FormattedOutput.h"
#include "math/RegularRankFourTensor.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
MP2EnergyCorrector<SCFMode>::MP2EnergyCorrector(std::shared_ptr<SystemController> systemController,
                                                const double ssScaling, const double osScaling)
  : _systemController(systemController), _ssScaling(ssScaling), _osScaling(osScaling) {
  assert(_systemController);
}

template<>
double MP2EnergyCorrector<Options::SCF_MODES::RESTRICTED>::calculateElectronicEnergy() {
  const auto& basisController = _systemController->getBasisController();
  const unsigned int nBasisFunc = basisController->getNBasisFunctions();
  const auto& orbitals = _systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>();
  const auto& orbitalEnergies = orbitals->getEigenvalues();
  double MP2EnergyCorrection = 0.0;
  const unsigned int nocc = _systemController->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  const unsigned int nvirt = nBasisFunc - nocc;
  RegularRankFourTensor<double> eris(nBasisFunc, 0.0);

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, basisController->getPrescreeningThreshold());

  auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                 const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId; // no warnings, please
    eris(b, a, i, j) = integral(0);
    eris(b, a, j, i) = integral(0);
    eris(a, b, j, i) = integral(0);
    eris(a, b, i, j) = integral(0);
    eris(i, j, b, a) = integral(0);
    eris(i, j, a, b) = integral(0);
    eris(j, i, b, a) = integral(0);
    eris(j, i, a, b) = integral(0);
  };
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients =
      _systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  looper.loop(storeERIS, coefficients.lpNorm<Eigen::Infinity>());

  Ao2MoTransformer ao2mo(basisController);

  ao2mo.transformTwoElectronIntegrals(eris, eris, coefficients, nBasisFunc);

  Eigen::setNbThreads(1);
  Eigen::MatrixXd eigenvalues = Eigen::MatrixXd::Zero(nvirt, nvirt);
  eigenvalues.colwise() -= orbitalEnergies.segment(nocc, nvirt);
  eigenvalues.rowwise() -= orbitalEnergies.segment(nocc, nvirt).transpose();
#pragma omp parallel shared(eris)
  {
    /*
     * Summation
     */
#pragma omp for schedule(dynamic) reduction(+ : MP2EnergyCorrection)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = i; j < nocc; ++j) {
        for (unsigned int kappa = nocc; kappa < nBasisFunc; ++kappa) {
          for (unsigned int iota = nocc; iota < nBasisFunc; ++iota) {
            double t_ij_ab = eris(i, iota, j, kappa) /
                             (eigenvalues(kappa - nocc, iota - nocc) + orbitalEnergies(i) + orbitalEnergies(j));
            double t_ij_ba = eris(i, kappa, j, iota) /
                             (eigenvalues(kappa - nocc, iota - nocc) + orbitalEnergies(i) + orbitalEnergies(j));
            double ssEnergy = (i == j ? 1.0 : 2.0) * (t_ij_ab - t_ij_ba) * eris(i, iota, j, kappa);
            double osEnergy = (i == j ? 1.0 : 2.0) * t_ij_ab * eris(i, iota, j, kappa);
            MP2EnergyCorrection += _ssScaling * ssEnergy + _osScaling * osEnergy;
          }
        }
      }
    }
  }
  Eigen::setNbThreads(0);
  return MP2EnergyCorrection;
}

template class MP2EnergyCorrector<Options::SCF_MODES::RESTRICTED>;
// template class MP2EnergyCorrector<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
