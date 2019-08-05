/**
 * @file   MP2.cpp
 *
 * @date   Jul 14, 2014
 * @author Jan Unsleber
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
#include "postHF/MPn/MP2.h"
/* Include Serenity Internal Headers */
#include "integrals/transformer/Ao2MoTransformer.h"
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/ElectronicStructure.h"
#include "io/FormattedOutput.h"
#include "data/OrbitalController.h"
#include "math/RegularRankFourTensor.h"
#include "system/SystemController.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"


namespace Serenity {
using namespace std;

template <Options::SCF_MODES T>
MP2EnergyCorrector<T>::MP2EnergyCorrector(
      shared_ptr<SystemController> systemController,
      const double ssScaling,
      const double osScaling) :
    _systemController(systemController),
    _ssScaling(ssScaling),
    _osScaling(osScaling){
  assert(_systemController);
}

//template <Options::SCF_MODES T>
template<>
double MP2EnergyCorrector<Options::SCF_MODES::RESTRICTED>::calculateElectronicEnergy(){

  const auto& basisController = _systemController->getBasisController();
  const auto& nElectrons =_systemController->getNElectrons<Options::SCF_MODES::RESTRICTED>();
  const unsigned int nBasisFunc = basisController->getNBasisFunctions();
  const auto& orbitals = _systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>();
  const auto& orbitalEnergies = orbitals->getEigenvalues();
  double MP2EnergyCorrection = 0.0;
  RegularRankFourTensor<double> eris(nBasisFunc, 0.0);


  TwoElecFourCenterIntLooper looper(libint2::Operator::coulomb,0,basisController,1E-10);

  auto const storeERIS = [&eris]
                          (const unsigned int&  a,
                              const unsigned int&  b,
                              const unsigned int&  i,
                              const unsigned int&  j,
                              const Eigen::VectorXd&  integral,
                              const unsigned int threadId) {
    (void)threadId; //no warnings, please
    eris(b,a,i,j) = integral(0);
    eris(b,a,j,i) = integral(0);
    eris(a,b,j,i) = integral(0);
    eris(a,b,i,j) = integral(0);
    eris(i,j,b,a) = integral(0);
    eris(i,j,a,b) = integral(0);
    eris(j,i,b,a) = integral(0);
    eris(j,i,a,b) = integral(0);
  };
  looper.loop(storeERIS);

  Ao2MoTransformer ao2mo(basisController);

  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients=_systemController->
      getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();

  ao2mo.transformTwoElectronIntegrals(eris,eris,coefficients,nBasisFunc);

#pragma omp parallel shared(eris)
  {
  /*
   * Summation
   */
#pragma omp for schedule(dynamic) reduction(+:MP2EnergyCorrection)
  for (unsigned int beta = 0; beta < (nElectrons/2); ++beta){
    for (unsigned int alpha = 0; alpha < (nElectrons/2); ++alpha){
      for (unsigned int kappa = (nElectrons/2); kappa < nBasisFunc; ++kappa){
        for (unsigned int iota = (nElectrons/2); iota < nBasisFunc; ++iota){
          double tmp = eris(alpha,iota,beta,kappa)*(2.0*eris(iota,alpha,kappa,beta)-eris(iota,beta,kappa,alpha));
          tmp = tmp/(orbitalEnergies[alpha]+orbitalEnergies[beta]-orbitalEnergies[iota]-orbitalEnergies[kappa]);
          MP2EnergyCorrection += tmp;
        }
      }
    }
  }
  
  }

  MP2EnergyCorrection *= 0.5 * (_ssScaling + _osScaling);

  return MP2EnergyCorrection;
}

template class MP2EnergyCorrector<Options::SCF_MODES::RESTRICTED>;
//template class MP2EnergyCorrector<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
