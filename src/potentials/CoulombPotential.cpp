/**
 * @file   CoulombPotential.cpp
 *
 * @date   Nov 24, 2016
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
#include "potentials/CoulombPotential.h"
/* Include Serenity Internal Headers */
#include "integrals/RI_J_IntegralController.h"
#include "misc/Timing.h"


namespace Serenity {

template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
CoulombPotential<SCFMode>::getMatrix(){
  Timings::takeTime("Active System -    Coulomb Pot.");
  if (!_potential){
    const auto& densityMatrix = _dMatController->getDensityMatrix();
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot){
      pot_spin.setZero();
    };
    this->addToMatrix(*_potential,densityMatrix);
  }
  Timings::timeTaken("Active System -    Coulomb Pot.");
  return *_potential;
}

template <Options::SCF_MODES SCFMode> double
CoulombPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  if (!_potential) this->getMatrix();
  Timings::takeTime("Active System -    Coulomb Pot.");
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot,P){
    energy += 0.5*pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -    Coulomb Pot.");
  return energy;
};

template <Options::SCF_MODES SCFMode>  void
CoulombPotential<SCFMode>::addToMatrix(
    FockMatrix<SCFMode>& F,
    const DensityMatrix<SCFMode>& deltaP){

  Eigen::MatrixXd totalDensity = deltaP.total();
  const unsigned int nAuxFunctions = _ri_j_IntController->getAuxBasisController()->getNBasisFunctions();
  /*
   * Calculation of the rightmost double sum in the referenced equation ([1].(4)) which is over
   * (Q|kappa,lamda)*D_{kappa,lamda}.
   * Calculate sum_{mu,nu} (P|mu,nu)*Density(mu,nu) for each auxiliary basis function P
   */

  // Maximum absolute value in densityMatrix
  const double maxDens = totalDensity.lpNorm<Eigen::Infinity>();
  /*
   * Detailed prescreening function
   */
  auto prescreeningFunc = [&](
      unsigned int i, unsigned int j,
      unsigned int nI, unsigned int nJ, double schwartz) {
    (void)i;
    (void)j;
    (void)nI;
    (void)nJ;
    /*
     * Shell quadruple is significant.
     */
    if (fabs(maxDens*schwartz) < _prescreeningThreshold) return true;
    /*
     * Shell quadruple is insignificant.
     */
    return false;
  };

#ifdef _OPENMP
  Eigen::MatrixXd sumMat(nAuxFunctions,omp_get_max_threads());
#else
  Eigen::MatrixXd sumMat(nAuxFunctions,1);
#endif
  sumMat.setZero();
  auto densptr = totalDensity.data();
  auto sumptr = sumMat.data();
  const unsigned int nBFs = totalDensity.cols();
  auto const loopEvalFunction = [&densptr,&sumptr,&nAuxFunctions,&nBFs]
                                 (const unsigned int& i,
                                     const unsigned int& j,
                                     const unsigned int& K,
                                     const double& integral,
                                     const unsigned int threadId){
    sumptr[K+nAuxFunctions*threadId] +=  (i == j ? 1.0 : 2.0) * integral * densptr[i*nBFs+j];
  };
  _ri_j_IntController->loopOver3CInts(loopEvalFunction,prescreeningFunc);
  Eigen::VectorXd sumPMuNu_DMuNu = sumMat.rowwise().sum();
  sumMat.resize(1,1);
  /*
   *  Done calculating sumPMuNu_DMuNu
   */

  /*
   * Calculate product of sumPMuNu_DMuNu and inverse of 2c2e matrix M
   */

  Eigen::VectorXd coefficients = _ri_j_IntController->getInverseM()*sumPMuNu_DMuNu.eval();

  /*
   * The actual calculation of the resulting matrix elements, which is now basically a scalar
   * product of two vectors for each {mu, nu}.
   * potential(mu,nu) = sum_P coefficients(P)*(P|mu,nu)
   */
#ifdef _OPENMP
  std::vector<MatrixInBasis<RESTRICTED> > potential(omp_get_max_threads(),MatrixInBasis<RESTRICTED>(this->_basis));
#else
  std::vector<MatrixInBasis<RESTRICTED> > potential(1,MatrixInBasis<RESTRICTED>(this->_basis));
#endif
  takeTime("Coulomb 4");
  auto coeff = coefficients.data();
  auto const loopEvalFunction2 = [&coeff,&potential,&nBFs]
                                  (const unsigned int& i,
                                      const unsigned int& j,
                                      const unsigned int& K,
                                      const double& integral,
                                      const unsigned int threadId){
    (void)threadId; // no warning
    if (i==j){
      potential[threadId].data()[i*nBFs+j] += integral * coeff[K];
    } else {
      potential[threadId].data()[i*nBFs+j] += integral * coeff[K];
      potential[threadId].data()[j*nBFs+i] += integral * coeff[K];
    }
  };
  _ri_j_IntController->loopOver3CInts(loopEvalFunction2,prescreeningFunc);
  timeTaken(3,"Coulomb 4");
  for_spin(F){
#ifdef _OPENMP
  for (unsigned int i=0; i<(unsigned int)omp_get_max_threads(); i++){
    F_spin += potential[i];
  }
#else
  F_spin += potential[0];
#endif
  };
}

template<Options::SCF_MODES SCFMode>
std::vector<unsigned int> CoulombPotential<SCFMode>::createBasisToAtomIndexMapping(
      const std::vector<std::pair<unsigned int, unsigned int> >& basisIndicesRed,
      unsigned int nBasisFunctionsRed) {
  std::vector<unsigned int> mapping(nBasisFunctionsRed);
  // Vector to check whether ALL basis function shells are assigned to an atom index
  std::vector<bool> hasElementBeenSet(nBasisFunctionsRed, false);
  for (unsigned int iAtom=0; iAtom < basisIndicesRed.size(); ++iAtom) {
    const unsigned int firstIndex = basisIndicesRed[iAtom].first;
    const unsigned int endIndex = basisIndicesRed[iAtom].second;
    for (unsigned int iShell=firstIndex; iShell < endIndex; ++iShell) {
      mapping[iShell] = iAtom;
      hasElementBeenSet[iShell] = true;
    }
  }
  // Check
  for (bool x : hasElementBeenSet){
    if (not x)
      throw SerenityError("CoulombPotential: Missed basis function in mapping procedure.");
  }
  return mapping;
}

template <Options::SCF_MODES SCFMode> Eigen::MatrixXd
CoulombPotential<SCFMode>::getGeomGradients(){

  const auto& orbitalSet = _actSystem->getActiveOrbitalController<SCFMode>();
  auto atoms = _actSystem->getAtoms();
  unsigned int nAtoms = atoms.size();
  Matrix<double> eriContr(nAtoms,3);
  Matrix<double> exchangeContr(nAtoms,3);
  eriContr.setZero();
  exchangeContr.setZero();

  // Collecting all the necessities, Mapping Basis and AuxBasis to atoms

  auto basisController = _actSystem->getAtomCenteredBasisController();
  auto auxBasisController = _actSystem->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);

  unsigned int nBasisFunctionsRed = orbitalSet->getBasisController()->getReducedNBasisFunctions();
  auto basisIndicesRed = _actSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
  auto mapping = createBasisToAtomIndexMapping(basisIndicesRed,nBasisFunctionsRed);

  auto& auxBasis = auxBasisController->getBasis();
  auto nAuxBasFunc = auxBasisController->getNBasisFunctions();
  auto nAuxBasFuncRed = auxBasisController->getReducedNBasisFunctions();
  auto auxBasisIndicesRed = auxBasisController->getBasisIndicesRed();
  auto auxMapping = createBasisToAtomIndexMapping(auxBasisIndicesRed,nAuxBasFuncRed);

  DensityMatrix<RESTRICTED> densityMatrix(_actSystem->getElectronicStructure<SCFMode>()->getDensityMatrix().total());


  // Get the existing RI_J Controller from Config (there should be one from the first SCF)


  auto ri_j_IntController = RI_J_IntegralControllerFactory::getInstance()
      .produce(basisController,auxBasisController);

  auto inverseM = ri_j_IntController->getInverseM();

  /*
   * Constructing the different parts of eq. ([2].[6]).
   *
   * Starting with the C vectors. The first step combines the
   * Density Matrix elements ij with the respecting 3c Integral <J|ij>.
   * This is technically a Matrix-Vector multiplication.
   *
   */
#ifdef _OPENMP
  // create a vector of vectors for each thread
  std::vector<Eigen::VectorXd> CvecPriv(omp_get_max_threads(),Eigen::VectorXd::Zero(nAuxBasFunc));
#else
  // or just one
  std::vector<Eigen::VectorXd> CvecPriv(1,Eigen::VectorXd::Zero(nAuxBasFunc));
#endif


  TwoElecThreeCenterIntLooper looper(libint2::Operator::coulomb,
      0,
      basisController,
      auxBasisController,
      1E-10);

  Eigen::VectorXd Cvec(nAuxBasFunc);
  Cvec.setZero();
  auto const calcCVector = [&CvecPriv,&densityMatrix]
                            (const unsigned int& i,
                                const unsigned int& j,
                                const unsigned int& K,
                                Eigen::VectorXd& intValues,
                                const unsigned int threadId) {

    const double perm = (i == j ? 1.0 : 2.0);
    CvecPriv[threadId][K] += perm * intValues(0) * densityMatrix(i,j);

  };
  looper.loop(calcCVector);

#ifdef _OPENMP
for (unsigned int i =0;i<(unsigned int)omp_get_max_threads();++i){
  Cvec += CvecPriv[i];
}
#else
  Cvec += CvecPriv[0];
#endif

  /*
   * The vector resulting from this Matrix * Vector multiplication is then
   * part of another Matrix-Vector Multiplication with the inverted 2c Integral Matrix (inverseM)
   */

  Eigen::VectorXd Dvec = inverseM*Cvec;

  /*
   * The resulting Vector Dvec is the equivalent of the C^J vector in the paper.
   * The following set of loops calculates the complete 3rd part of eq ([2].[6])
   * (easily recognizable by the sign).
   *
   * This is done analogously to the calculation of the inverseM in the RI_J_IntegralController,
   * but with the added cartesian directions. Note that to get a certain part of an Eigen Matrix,
   * the method .block() is invoked, the equivalent for a vector is .segment().
   *
   */
  Libint& libint = Libint::getInstance();
  libint.initialize(libint2::Operator::coulomb,1,2);
#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> eriContrPriv(omp_get_max_threads(),Eigen::MatrixXd::Zero(nAtoms,3));
#else
  // or just one
  std::vector<Eigen::MatrixXd> eriContrPriv(1,Eigen::MatrixXd::Zero(nAtoms,3));
#endif

  for (unsigned int auxJ = 0; auxJ < auxBasis.size(); auxJ++){
    const unsigned int nAuxJ = auxBasis[auxJ]->getNContracted();
    unsigned int offauxJ = auxBasisController->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
    for (unsigned int auxK = 0; auxK <= auxJ; auxK++){
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      Eigen::MatrixXd intDerivs;
      const unsigned int nAuxK = auxBasis[auxK]->getNContracted();
      unsigned int offauxK = auxBasisController->extendedIndex(auxK);

      if (libint.compute(libint2::Operator::coulomb,
          1,
          *auxBasis[auxJ],
          *auxBasis[auxK],
          intDerivs)){

        double perm = (auxJ == auxK ? 1.0 : 2.0);

        Eigen::MatrixXd prefac = 0.5 * perm * Dvec.segment(offauxJ,nAuxJ) * Dvec.segment(offauxK,nAuxK).transpose();

        for (unsigned int iAtom = 0; iAtom < 2; ++iAtom) {

          unsigned int nAtom;
          switch(iAtom){
            case(0): nAtom = auxMapping[auxJ]; break;
            case(1): nAtom = auxMapping[auxK]; break;
          }

          for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
            Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom*3+iDirection).data(),nAuxK,nAuxJ);
            eriContrPriv[threadId](nAtom,iDirection) -= tmp.transpose().cwiseProduct(prefac).sum();
          }
        }
      }

    }
  }
  libint.finalize(libint2::Operator::coulomb,1,2);
  /*
   *
   * The following looper calculates the complete 1st and 2nd part of eq ([2].[6])
   * (easily recognizable by the sign).
   *
   */

  TwoElecThreeCenterIntLooper derivLooper(libint2::Operator::coulomb,
      1,
      basisController,
      auxBasisController,
      1E-10);

  auto const add3cDerivs = [&eriContrPriv,&Dvec,&densityMatrix,&mapping,&auxMapping,
              &basisController,&auxBasisController]
                            (const unsigned int& i,
                                const unsigned int& j,
                                const unsigned int& K,
                                Eigen::VectorXd& intValues,
                                const unsigned int threadId) {
    double perm = (i == j ? 1.0 : 2.0);
    for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
      for (unsigned int iAtom = 0; iAtom < 3; ++iAtom) {
        unsigned int nAtom;
        switch(iAtom){
          case(0): nAtom = auxMapping[auxBasisController->reducedIndex(K)]; break;
          case(1): nAtom = mapping[basisController->reducedIndex(i)]; break;
          case(2): nAtom = mapping[basisController->reducedIndex(j)]; break;
        }
          eriContrPriv[threadId](nAtom,iDirection) += perm * Dvec[K] * densityMatrix(i,j)
                                       * intValues(iAtom*3+iDirection);
      }
    }

  };
  derivLooper.loop(add3cDerivs);

#ifdef _OPENMP
  for (unsigned int i =0;i<(unsigned int)omp_get_max_threads();++i){
    eriContr += eriContrPriv[i];
  }
#else
      eriContr = eriContrPriv[0];
#endif
  return eriContr;

}


template class CoulombPotential<Options::SCF_MODES::RESTRICTED>;
template class CoulombPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
