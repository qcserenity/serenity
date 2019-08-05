/**
 * @file ExchangeInteractionPotential.cpp
 * @author: Kevin Klahr
 *
 * @date 30. November 2016
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
#include "potentials/ExchangeInteractionPotential.h"


namespace Serenity {

template <Options::SCF_MODES SCFMode> double
ExchangeInteractionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  if (!_potential) this->getMatrix();
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot,P){
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  return energy;
};

template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
ExchangeInteractionPotential<SCFMode>::getMatrix(){
  if (!_potential){

    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& F = *_potential;
    for_spin(F){
      F_spin.setZero();
    };
    for (unsigned int i=0;i<_dMatControllers.size();++i){
#ifdef _OPENMP
      // create a vector of matrices for each thread
      std::vector<MatrixInBasis<SCFMode> > eriContr;
      for (int j=0; j<omp_get_max_threads(); ++j) {
        eriContr.push_back(MatrixInBasis<SCFMode>(this->_basis));
        auto& reference = eriContr[j];
        for_spin(reference){
        (reference_spin).setZero();
        };
      }
#else
      // or just one
      std::vector<MatrixInBasis<SCFMode> > eriContr(1,MatrixInBasis<SCFMode>(this->_basis));
      auto& reference = eriContr[0];
      for_spin(reference){
      (reference_spin).setZero();
      };
#endif
      const auto& envMat = _dMatControllers[i]->getDensityMatrix();

              ExchangeInteractionIntLooper excLooper(libint2::Operator::coulomb, 0, _basis, envMat.getBasisController(),_prescreeningThreshold);

              auto const excLooperFunction = [&]
                                           (const unsigned int&  i,
                                               const unsigned int&  a,
                                               const unsigned int&  j,
                                               const unsigned int&  b,
                                               Eigen::VectorXd& intValues,
                                               const unsigned int& threadID) {
                double perm = 1.0;
                // res has to be 0.5 in the restricted case to compensate double counting
                double res = 1.0;
                if (SCFMode == Options::SCF_MODES::RESTRICTED) res *= 0.5;
//                if (a!=b && i==j) perm *= 2.0;
                auto& reference = eriContr[threadID];
                for_spin(reference, envMat){

                /*
                 * Exchange contribution
                 */

                  const double exc = res * perm * envMat_spin(a,b) * intValues(0) * _exc;
                  (reference_spin)(i,j) -= exc;
//                  if (i!=j) eriContr[threadID](j,i) -= exc;
                };
              };
              excLooper.loop(excLooperFunction);


#ifdef _OPENMP
      // sum over all threads
      for (int j=0; j<omp_get_max_threads(); ++j) {
        auto& reference = eriContr[j];
        for_spin(F,reference){
          F_spin += reference_spin;
        };
      }
#else
      auto& reference = eriContr[0];
      for_spin(F,reference){
        F_spin += reference_spin;
      };
#endif
    }
  }
  return *_potential;
}


template <Options::SCF_MODES SCFMode> Eigen::MatrixXd
ExchangeInteractionPotential<SCFMode>::getGeomGradients(){
  Eigen::MatrixXd gradientContr(1,3);
  gradientContr.setZero();
  assert(false && "No gradients available for the experimental ExchangeInteractionPotential.");
  return gradientContr;
}

template class ExchangeInteractionPotential<Options::SCF_MODES::RESTRICTED>;
template class ExchangeInteractionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
