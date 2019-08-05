/**
 * @file ABCoreHamiltonian.cpp
 *
 * @date May 15, 2018
 * @author Moritz Bensberg
 *
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

#include "potentials/ABFockMatrixConstruction/ABCoreHamiltonian.h"
#include "potentials/ABFockMatrixConstruction/ABEffectiveCorePotential.h"
#include "geometry/Geometry.h"

namespace Serenity {

template <Options::SCF_MODES SCFMode>
ABCoreHamiltonian<SCFMode>::ABCoreHamiltonian(
    std::shared_ptr<BasisController> basisA,
    std::shared_ptr<BasisController> basisB,
    std::shared_ptr<Geometry> geometry) :
    ABPotential<SCFMode>(basisA,basisB),
    _geom(geometry){
  // Setting the notifying system up.
  this->_basisA->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_basisB->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  _abEffectiveCorePotential = std::make_shared<ABEffectiveCorePotential<SCFMode> > (
      this->_basisA,
      this->_basisB,
      _geom->getAtoms());
}

template <Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>& ABCoreHamiltonian<SCFMode>::getMatrix() {
  if (!_abPotential) {
    unsigned int nBasisA = this->_basisA->getNBasisFunctions();
    unsigned int nBasisB = this->_basisB->getNBasisFunctions();
    const auto& basisA = this->_basisA->getBasis();
    const auto& basisB = this->_basisB->getBasis();
    // initialize fock matrix
    _abPotential.reset(new SPMatrix<SCFMode>(nBasisA,nBasisB));
    SPMatrix<SCFMode>& f_AB = *_abPotential;
    // add ECP contribution. (Zero if no ECPs are used)
    f_AB += _abEffectiveCorePotential->getMatrix();
    // Libint operators
    libint2::Operator opKin = libint2::Operator::kinetic;
    libint2::Operator opNuc = libint2::Operator::nuclear;
    /* ======================= */
    /* 1. kinetic contribution */
    /* ======================= */
    _libint->initialize(opKin,0,2);
  #pragma omp parallel
    {
      Eigen::MatrixXd ints;
  #pragma omp for schedule(static, 1)
      for(unsigned int a = 0; a < basisA.size(); ++a){
        unsigned int offA = this->_basisA->extendedIndex(a);
        const unsigned int nA = basisA[a]->getNContracted();
        for (unsigned int b = 0; b < basisB.size(); ++b){
          unsigned int offB = this->_basisB->extendedIndex(b);
          const unsigned int nB = basisB[b]->getNContracted();
          if(_libint->compute(opKin,0,*basisA[a],*basisB[b],ints)){
            Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(),nB,nA); // nA x nB
            for_spin(f_AB) {
              f_AB_spin.block(offA,offB,nA,nB)+=tmp.transpose();
            };
          }
        }
      }
    } /* END OpenMP parallel */
    _libint->finalize(opKin,0,2);

    /* ======================= */
    /* 2. nuclear contribution */
    /* ======================= */
    _libint->initialize(opNuc,0,2,_geom->getAtoms());
  #pragma omp parallel
    {
      Eigen::MatrixXd ints;
  #pragma omp for schedule(static, 1)
      for(unsigned int a = 0; a < basisA.size(); ++a){
        unsigned int offA = this->_basisA->extendedIndex(a);
        const unsigned int nA = basisA[a]->getNContracted();
        for (unsigned int b = 0; b < basisB.size(); ++b){
          unsigned int offB = this->_basisB->extendedIndex(b);
          const unsigned int nB = basisB[b]->getNContracted();
          if(_libint->compute(opNuc,0,*basisA[a],*basisB[b],ints)){
            Eigen::Map<Eigen::MatrixXd> tmp(ints.col(0).data(),nB,nA); // nA x nB
            for_spin(f_AB) {
              f_AB_spin.block(offA,offB,nA,nB)+=tmp.transpose();
            };
          }
        }
      }
    } /* END OpenMP parallel */
    _libint->finalize(opNuc,0,2);
  } /* if !_abPotential */
  return *_abPotential;
}


template class ABCoreHamiltonian<Options::SCF_MODES::RESTRICTED>;
template class ABCoreHamiltonian<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
