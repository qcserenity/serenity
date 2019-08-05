/**
 * @file RIMP2.cpp
 *
 * @date Aug 14, 2017
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
#include "postHF/MPn/RIMP2.h"
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "integrals/RI_J_IntegralController.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <unsupported/Eigen/MatrixFunctions>


namespace Serenity {

template <Options::SCF_MODES SCFMode>
RIMP2<SCFMode>::RIMP2(
    std::shared_ptr<SystemController> systemController,
    const double ssScaling,
    const double osScaling) :
   _systemController(systemController),
   _ssScaling(ssScaling),
   _osScaling(osScaling){
  assert(_systemController);
}

template <Options::SCF_MODES SCFMode>
double RIMP2<SCFMode>::calculateCorrection(){
/*
 * RI-MP2
 */

/*
 *  Gather data objects
 */
auto es = _systemController->getElectronicStructure<SCFMode>();
auto basisController = _systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
auto auxBasisController = _systemController->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
CoefficientMatrix<SCFMode> coeff = es->getMolecularOrbitals()->getCoefficients();

if (_systemController->getSettings().basis.auxCLabel == ""){
  printf("%4s Default Auxiliary Basis Set (Correlation):             %5"
      ""
      "s\n","",(_systemController->getSettings().basis.label + "-RI-C").c_str());
  print("");
}

/*
 * Get orbitals counts
 *  (n = defaut basis, m = aux. basis)
 */
auto nocc = _systemController->getNOccupiedOrbitals<SCFMode>();
const unsigned int n = basisController->getNBasisFunctions();
const unsigned int m = auxBasisController->getNBasisFunctions();

RI_J_IntegralController riints(basisController,auxBasisController);
Eigen::MatrixXd minv = riints.getInverseM().sqrt();
Eigen::MatrixXd ints = Eigen::MatrixXd::Zero(m, n*n);

// generate pointer for faster access
auto invM_ptr =  minv.data();
auto ints_ptr =  ints.data();

/*
 * Generate (ia|K)
 * Scaling [O(n^2 * m) = O(x^3)]
 */
auto const loopEvalFunction = [&ints_ptr,&m,&n]
                               (const unsigned int& i,
                                   const unsigned int& j,
                                   const unsigned int& K,
                                   Eigen::VectorXd& integral,
                                   const unsigned int threadId){
  (void)threadId;
  ints_ptr[(unsigned long long)K*n*n+i*n+j] = integral.data()[0];
  ints_ptr[(unsigned long long)K*n*n+j*n+i] = integral.data()[0];
};

TwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,0,
                                   basisController,auxBasisController,1E-10);
looper.loop(loopEvalFunction);

/*
 *  Transform (ia|K) to (ia|K)(K|Q)^{-1/2}
 *  Scaling [O(n^2 * m^2) = O(x^4)]
 */
#pragma omp parallel for schedule(static, 1)
for(unsigned int i=0;i<n;i++){
  for(unsigned int j=0;j<n;j++){
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(m);
    auto tmp_ptr = tmp.data();
    for(unsigned int Q=0;Q<m;Q++){
      for(unsigned int K=0;K<m;K++){
        tmp_ptr[Q] += invM_ptr[K*m+Q] * ints_ptr[(unsigned long long)K*n*n+i*n+j];
      }
    }
    for(unsigned int Q=0;Q<m;Q++){
      ints_ptr[(unsigned long long)Q*n*n+i*n+j] = tmp_ptr[Q];
    }
  }
}

/*
 * Transform (ia|K)(K|Q)^{-1/2} to B_{ia,Q} = (mu nu|K)(K|Q)^{-1/2}
 * Scaling [O(n^2 * m * nocc * nvirt) = O(x^5)]
 */
for_spin(coeff,_BiaQ,nocc){
  unsigned int nvirt = n-nocc_spin;
  auto coeff_ptr = coeff_spin.data();
  _BiaQ_spin.resize(m,nocc_spin*nvirt);
  auto BiaQ_ptr =  _BiaQ_spin.data();
  for(unsigned int mu=0;mu<nocc_spin;mu++){
    for(unsigned int nu=0;nu<nvirt;nu++){
#pragma omp parallel for schedule(static, 1)
      for(unsigned int Q=0;Q<m;Q++){
        double tmp = 0.0;
        for(unsigned int i=0;i<n;i++){
          for(unsigned int j=0;j<n;j++){
            tmp += coeff_ptr[i+n*mu] * coeff_ptr[j+n*(nu+nocc_spin)] * ints_ptr[(unsigned long long)Q*n*n+i*n+j];
          }
        }
        BiaQ_ptr[(unsigned long long)Q*nvirt*nocc_spin+mu*nvirt+nu] = tmp;
      }
    }
  }
};
// free some memory
ints.resize(0,0);

return this->calculateEnergy();
}

template <>
double RIMP2<UNRESTRICTED>::calculateEnergy(){

  /*
   * Get orbitals counts
   *  (n = defaut basis, m = aux. basis)
   */
  auto nocc = _systemController->getNOccupiedOrbitals<UNRESTRICTED>();
  const unsigned int n = _systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT)->getNBasisFunctions();
  const unsigned int m = _systemController->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL)->getNBasisFunctions();

  SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd > orbitalEnergies = _systemController->getElectronicStructure<UNRESTRICTED>()->getMolecularOrbitals()->getEigenvalues();

  /*
   * Calculate energy using B_{ia,Q}
   * Scaling [O(m * nocc^2 * nvirt^2) = O(x^5)]
   */

  // same spin integrals
  double MP2EnergyCorrection = 0.0;
  for_spin(_BiaQ,nocc,orbitalEnergies){
    auto BiaQ_ptr = _BiaQ_spin.data();
    const unsigned int nvirt = n-nocc_spin;
    double tmpEnergy = 0.0;
#pragma omp parallel for schedule(dynamic) reduction(+:tmpEnergy)
    for (unsigned int j = 0; j < nocc_spin; ++j){
      for (unsigned int i = 0; i < nocc_spin; ++i){
        for (unsigned int b = 0; b < nvirt; ++b){
          for (unsigned int a = 0; a < nvirt; ++a){
            double iajb = 0.0;
            double ibja = 0.0;
            for(unsigned int Q=0;Q<m;Q++){
              iajb += BiaQ_ptr[(unsigned long long)Q*nvirt*nocc_spin+i*nvirt+a]*BiaQ_ptr[(unsigned long long)Q*nvirt*nocc_spin+j*nvirt+b];
              ibja += BiaQ_ptr[(unsigned long long)Q*nvirt*nocc_spin+i*nvirt+b]*BiaQ_ptr[(unsigned long long)Q*nvirt*nocc_spin+j*nvirt+a];
            }
            double tmp = iajb *(iajb-ibja);
            tmp = tmp/(orbitalEnergies_spin[i]+orbitalEnergies_spin[j]-orbitalEnergies_spin[a+nocc_spin]-orbitalEnergies_spin[b+nocc_spin]);
            tmpEnergy += tmp;
          }
        }
      }
    }
    MP2EnergyCorrection += 0.5*tmpEnergy;
  };
  MP2EnergyCorrection *= _ssScaling;

  // mixed spin integrals
  const unsigned int nvirtA = n-nocc.alpha;
  const unsigned int nvirtB = n-nocc.beta;
  auto BiaQ_ptrA = _BiaQ.alpha.data();
  auto BiaQ_ptrB = _BiaQ.beta.data();
  double tmpEnergy = 0.0;
#pragma omp parallel for schedule(dynamic) reduction(+:tmpEnergy)
    for (unsigned int j = 0; j < nocc.alpha; ++j){
      for (unsigned int i = 0; i < nocc.beta; ++i){
        for (unsigned int b = 0; b < nvirtA; ++b){
          for (unsigned int a = 0; a < nvirtB; ++a){
            double iajb = 0.0;
            for(unsigned int Q=0;Q<m;Q++){
              iajb += BiaQ_ptrB[(unsigned long long)Q*nvirtB*nocc.beta+i*nvirtB+a]*BiaQ_ptrA[(unsigned long long)Q*nvirtA*nocc.alpha+j*nvirtA+b];
            }
            double tmp = iajb *iajb;
            tmpEnergy += tmp/(orbitalEnergies.beta[i]+orbitalEnergies.alpha[j]-orbitalEnergies.beta[a+nocc.beta]-orbitalEnergies.alpha[b+nocc.alpha]);
          }
        }
      }
    }
  MP2EnergyCorrection += _osScaling * tmpEnergy;

  return MP2EnergyCorrection;
}

template <>
double RIMP2<RESTRICTED>::calculateEnergy(){

  /*
   * Get orbitals counts
   *  (n = defaut basis, m = aux. basis)
   */
  auto nocc = _systemController->getNOccupiedOrbitals<RESTRICTED>();
  const unsigned int n = _systemController->getBasisController(Options::BASIS_PURPOSES::DEFAULT)->getNBasisFunctions();
  const unsigned int nvirt = n-nocc;
  const unsigned int m = _systemController->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL)->getNBasisFunctions();

  SpinPolarizedData<RESTRICTED, Eigen::VectorXd > orbitalEnergies = _systemController->getElectronicStructure<RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
  auto BiaQ_ptr = _BiaQ.data();

  /*
   * Calculate energy using B_{ia,Q}
   * Scaling [O(m * nocc^2 * nvirt^2) = O(x^5)]
   */

  double MP2EnergyCorrection = 0.0;

  #pragma omp parallel for schedule(dynamic) reduction(+:MP2EnergyCorrection)
  for (unsigned int j = 0; j < nocc; ++j){
    for (unsigned int i = 0; i < nocc; ++i){
      for (unsigned int b = 0; b < nvirt; ++b){
        for (unsigned int a = 0; a < nvirt; ++a){
          double iajb = 0.0;
          double ibja = 0.0;
          for(unsigned int Q=0;Q<m;Q++){
            iajb += BiaQ_ptr[(unsigned long long)Q*nvirt*nocc+i*nvirt+a]*BiaQ_ptr[(unsigned long long)Q*nvirt*nocc+j*nvirt+b];
            ibja += BiaQ_ptr[(unsigned long long)Q*nvirt*nocc+i*nvirt+b]*BiaQ_ptr[(unsigned long long)Q*nvirt*nocc+j*nvirt+a];
          }
          double tmp = iajb *(2.0*iajb-ibja);
          tmp = tmp/(orbitalEnergies[i]+orbitalEnergies[j]-orbitalEnergies[a+nocc]-orbitalEnergies[b+nocc]);
          MP2EnergyCorrection += tmp;
        }
      }
    }
  }

  MP2EnergyCorrection *= 0.5 * (_ssScaling + _osScaling);
  return MP2EnergyCorrection;
}



template class RIMP2<Options::SCF_MODES::RESTRICTED>;
template class RIMP2<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
