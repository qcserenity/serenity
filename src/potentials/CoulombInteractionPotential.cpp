/**
 * @file CoulombInteractionPotential.cpp
 * @author: Kevin Klahr
 *
 * @date 29. November 2016
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
#include "potentials/CoulombInteractionPotential.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "integrals/looper/CoulombInteractionIntLooper.h"
#include "basis/CustomBasisController.h"
#include "data/matrices/FockMatrix.h"
#include "misc/Timing.h"

namespace Serenity {

template <Options::SCF_MODES SCFMode>
CoulombInteractionPotential<SCFMode>::CoulombInteractionPotential(
    std::shared_ptr<SystemController> actSystem,
    std::vector<std::shared_ptr<SystemController> > envSystems,
    const std::shared_ptr<BasisController> actBasis,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > envDensityMatrixController,
    const std::shared_ptr<BasisController> actAuxBasis,
    std::vector<std::shared_ptr<BasisController> > envAuxBasis):
    Potential<SCFMode>(actBasis),
    _actSystem(actSystem),
    _envSystems(envSystems),
    _envDMatController(envDensityMatrixController),
    _actAuxBasis(actAuxBasis),
    _envAuxBasis(envAuxBasis){
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_actAuxBasis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  for (auto& mat  : envDensityMatrixController){
    mat->getDensityMatrix().getBasisController()->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  }
  for (auto& bas  : envAuxBasis){
    bas->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  }
  _mode = Options::DENS_FITS::RI;
};

template <Options::SCF_MODES SCFMode>
CoulombInteractionPotential<SCFMode>::CoulombInteractionPotential(
    std::shared_ptr<SystemController> actSystem,
    std::vector<std::shared_ptr<SystemController> > envSystems,
    const std::shared_ptr<BasisController> actBasis,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > envDensityMatrixController):
    Potential<SCFMode>(actBasis),
    _actSystem(actSystem),
    _envSystems(envSystems),
    _envDMatController(envDensityMatrixController),
    _actAuxBasis(nullptr),
    _envAuxBasis(0,nullptr){
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  for (auto& mat  : envDensityMatrixController){
    mat->getDensityMatrix().getBasisController()->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  }
  _mode = Options::DENS_FITS::NONE;
};

template <Options::SCF_MODES SCFMode> double
CoulombInteractionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  if (!_potential) this->getMatrix();
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot,P){
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  return energy;
};

template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
CoulombInteractionPotential<SCFMode>::getMatrix(){
  Timings::takeTime("FDE -        Coulomb Pot.");
  if (!_potential){

    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& F = *_potential;
    for_spin(F){
      F_spin.setZero();
    };

    if (_mode == Options::DENS_FITS::RI){
      for (unsigned int i=0;i<_envDMatController.size();++i){

        // active
        auto nActAuxBasFunc = _actAuxBasis->getNBasisFunctions();
        auto envBasis = _envDMatController[i]->getDensityMatrix().getBasisController();
        auto envAuxBasis = _envAuxBasis[i];
        // env
        auto nEnvAuxBasFunc = envAuxBasis->getNBasisFunctions();
        auto envDensityMatrix = _envDMatController[i]->getDensityMatrix().total();
        // super system
        Basis superBas;
        Basis superAuxBas;
        superBas.reserve( this->_basis->getBasis().size() + envBasis->getBasis().size() );
        superAuxBas.reserve( _actAuxBasis->getBasis().size() + envAuxBasis->getBasis().size() );
        superBas.insert( superBas.end(), this->_basis->getBasis().begin(), this->_basis->getBasis().end() );
        superBas.insert( superBas.end(), envBasis->getBasis().begin(), envBasis->getBasis().end() );
        superAuxBas.insert( superAuxBas.end(), _actAuxBasis->getBasis().begin(), _actAuxBasis->getBasis().end() );
        superAuxBas.insert( superAuxBas.end(), envAuxBasis->getBasis().begin(), envAuxBasis->getBasis().end() );

        auto superSystemBasisController = std::shared_ptr<BasisController>(new CustomBasisController(superBas,"tmpSuperBas"));
        auto superSystemAuxBasisController = std::shared_ptr<BasisController>(new CustomBasisController(superAuxBas,"tmpSuperAuxBas"));

#ifdef _OPENMP
        // create a vector of matrices for each thread
        Eigen::MatrixXd sumMat(nActAuxBasFunc,omp_get_max_threads());
        std::vector<MatrixInBasis<RESTRICTED> > eriContr(omp_get_max_threads(),MatrixInBasis<RESTRICTED> (this->_basis));
#else
        // or just one
        Eigen::MatrixXd sumMat(nActAuxBasFunc,1);
        std::vector<MatrixInBasis<RESTRICTED> > eriContr(1,MatrixInBasis<RESTRICTED> (this->_basis));
        eriContr[0].setZero();
#endif
        sumMat.setZero();

#ifdef _OPENMP
        Eigen::MatrixXd sumMatEnv(nEnvAuxBasFunc,omp_get_max_threads());
#else
        Eigen::MatrixXd sumMatEnv(nEnvAuxBasFunc,1);
#endif
        sumMatEnv.setZero();
        auto ri_j_IntController = RI_J_IntegralControllerFactory::getInstance()
        .produce(superSystemBasisController,superSystemAuxBasisController);
        const auto inverseM = ri_j_IntController->getInverseM();
        TwoElecThreeCenterIntLooper looper1(libint2::Operator::coulomb,
            0,
            envBasis,
            _actAuxBasis,
            1E-10);
        auto const loopEvalFunction1 = [&]
                                        (const unsigned int& i,
                                            const unsigned int& j,
                                            const unsigned int& K,
                                            Eigen::VectorXd& integral,
                                            const unsigned int threadId){
          const double perm = (i == j ? 1.0 : 2.0);
          sumMat(K,threadId) +=  perm * integral[0] * envDensityMatrix(i,j);
        };
        looper1.loop(loopEvalFunction1);
        TwoElecThreeCenterIntLooper looper2(libint2::Operator::coulomb,
            0,
            envBasis,
            envAuxBasis,
            1E-10);
        auto const loopEvalFunction2 = [&]
                                        (const unsigned int& i,
                                            const unsigned int& j,
                                            const unsigned int& K,
                                            Eigen::VectorXd& integral,
                                            const unsigned int threadId){
          const double perm = (i == j ? 1.0 : 2.0);
          sumMatEnv(K,threadId) +=  perm * integral[0] * envDensityMatrix(i,j);
        };
        looper2.loop(loopEvalFunction2);
        Eigen::VectorXd sumPMuNu_DMuNuAct = sumMat.rowwise().sum();
        sumMat.resize(1,1);
        Eigen::VectorXd sumPMuNu_DMuNuEnv = sumMatEnv.rowwise().sum();
        sumMatEnv.resize(1,1);
        Eigen::VectorXd sumPMuNu_DMuNu(sumPMuNu_DMuNuAct.rows()+sumPMuNu_DMuNuEnv.rows());
        sumPMuNu_DMuNu << sumPMuNu_DMuNuAct,sumPMuNu_DMuNuEnv;
        Eigen::VectorXd coefficients = inverseM*sumPMuNu_DMuNu.eval();
        TwoElecThreeCenterIntLooper looper3(libint2::Operator::coulomb,
            0,
            this->_basis,
            _actAuxBasis,
            1E-10);
        auto const loopEvalFunction3 = [&]
                                        (const unsigned int& i,
                                            const unsigned int& j,
                                            const unsigned int& K,
                                            Eigen::VectorXd& integral,
                                            const unsigned int threadID){
          eriContr[threadID](i,j) += integral[0] * coefficients(K);
          if (i!=j){
            eriContr[threadID](j,i) +=  integral[0] * coefficients(K);
          }
        };
        looper3.loop(loopEvalFunction3);
        TwoElecThreeCenterIntLooper looper4(libint2::Operator::coulomb,
            0,
            this->_basis,
            envAuxBasis,
            1E-10);
        auto const loopEvalFunction4 = [&]
                                        (const unsigned int& i,
                                            const unsigned int& j,
                                            const unsigned int& K,
                                            Eigen::VectorXd& integral,
                                            const unsigned int threadID){
          eriContr[threadID](i,j) += integral[0] * coefficients(K+nActAuxBasFunc);
          if (i!=j){
            eriContr[threadID](j,i) +=  integral[0] * coefficients(K+nActAuxBasFunc);
          }
        };
        looper4.loop(loopEvalFunction4);
#ifdef _OPENMP
        // sum over all threads
        for (int t=0; t<omp_get_max_threads(); ++t) {
          for_spin(F){
            F_spin += eriContr[t];
          };
        }
#endif
      }
    } else {
      for (unsigned int i=0;i<_envDMatController.size();++i){
#ifdef _OPENMP
        // create a vector of matrices for each thread
        std::vector<MatrixInBasis<RESTRICTED> > eriContr;
        for (int t=0; t<omp_get_max_threads(); ++t) {
          eriContr.push_back(MatrixInBasis<RESTRICTED>(this->_basis));
          eriContr[t].setZero();
        }
#else
        // or just one
        std::vector<MatrixInBasis> eriContr(1,MatrixInBasis(this->_basis));
        eriContr[0].setZero();
#endif
        auto envMat = _envDMatController[i]->getDensityMatrix().total();
        CoulombInteractionIntLooper coulLooper(libint2::Operator::coulomb, 0, this->_basis, envMat.getBasisController(), 1E-10);

        auto const coulLooperFunction = [&eriContr,&envMat]
                                         (const unsigned int&  i,
                                             const unsigned int&  j,
                                             const unsigned int&  a,
                                             const unsigned int&  b,
                                             Eigen::VectorXd& intValues,
                                             const unsigned int& threadID) {
          double perm = 1.0;
          if (a!=b) perm *= 2.0;
          /*
           * Coulomb contribution
           */
          const double coul = perm * envMat(a,b) * intValues(0);
          eriContr[threadID](i,j) += coul;
          if (i!=j) eriContr[threadID](j,i) += coul;
        };
        coulLooper.loop(coulLooperFunction);

#ifdef _OPENMP
        // sum over all threads
        for (int t=0; t<omp_get_max_threads(); ++t) {
          for_spin(F){
            F_spin += eriContr[t];
          };
        }
#else
        for_spin(F){
          F_spin += eriContr[0];
        };
#endif
      }
    }
  }
  Timings::timeTaken("FDE -        Coulomb Pot.");
  return *_potential;
}

template<Options::SCF_MODES T>
std::vector<unsigned int> CoulombInteractionPotential<T>::createBasisToAtomIndexMapping(
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
      throw SerenityError("CoulombInteractionPotential: Missed basis function in mapping procedure.");
  }
  return mapping;
}

template <Options::SCF_MODES SCFMode> Eigen::MatrixXd
CoulombInteractionPotential<SCFMode>::getGeomGradients(){

  auto atomsAct = _actSystem->getAtoms();
  unsigned int nAtomsAct = atomsAct.size();
  Matrix<double> activeSystemGradientContr(atomsAct.size(),3);
  activeSystemGradientContr.setZero();
  DensityMatrix<RESTRICTED> densityMatrix(_actSystem->getElectronicStructure<SCFMode>()->getDensityMatrix().total());

  for (unsigned int i=0;i<_envSystems.size();++i){

    // Creating supersystem Geometry and Basis, necessary still for inverseM

    auto atomsEnv = _envSystems[i]->getAtoms();
    std::vector<std::shared_ptr<Atom> > ActiveFrozenAtoms;
    for (auto& atom : atomsAct){
      ActiveFrozenAtoms.push_back(atom);
    }
    for (auto& atom : atomsEnv){
      ActiveFrozenAtoms.push_back(atom);
    }
    unsigned int nAtoms = ActiveFrozenAtoms.size();
    Matrix<double> eriContr(nAtoms,3);
    eriContr.setZero();
    double test = 0.0;

    auto superSystemGeometry = std::make_shared<Geometry>(ActiveFrozenAtoms);
    auto superSystemBasisController = AtomCenteredBasisControllerFactory::produce(
        superSystemGeometry,
        _actSystem->getSettings().basis.basisLibPath,
        _actSystem->getSettings().basis.makeSphericalBasis,
        true,
        _actSystem->getSettings().basis.firstECP,
        _actSystem->getSettings().basis.label);
    auto superSystemAuxBasisController = AtomCenteredBasisControllerFactory::produce(
        superSystemGeometry,
        _actSystem->getSettings().basis.basisLibPath,
        _actSystem->getSettings().basis.makeSphericalBasis,
        false,
        999999999,
        _actSystem->getSettings().basis.auxJLabel);

    // Getting necessary information from the active and environment system
    // Active
    auto basisController = _actSystem->getAtomCenteredBasisController();
    auto auxBasisController = _actSystem->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
    auto& auxBasisFunctions = auxBasisController->getBasis();
    unsigned int nBasisFunctionsRed = basisController->getReducedNBasisFunctions();
    auto basisIndicesRed = _actSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
    auto mapping = createBasisToAtomIndexMapping(basisIndicesRed,nBasisFunctionsRed);
    auto nAuxBasFunc = auxBasisController->getNBasisFunctions();
    auto nAuxBasFuncRed = auxBasisController->getReducedNBasisFunctions();
    auto auxBasisIndicesRed = auxBasisController->getBasisIndicesRed();
    auto auxMapping = createBasisToAtomIndexMapping(auxBasisIndicesRed,nAuxBasFuncRed);

    // Environment
    auto basisControllerEnv = _envSystems[i]->getAtomCenteredBasisController();
    auto auxBasisControllerEnv = _envSystems[i]->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
    auto& auxBasisFunctionsEnv = auxBasisControllerEnv->getBasis();
    unsigned int nBasisFunctionsRedEnv = basisControllerEnv->getReducedNBasisFunctions();
    auto basisIndicesRedEnv = _envSystems[i]->getAtomCenteredBasisController()->getBasisIndicesRed();
    auto mappingEnv = createBasisToAtomIndexMapping(basisIndicesRedEnv,nBasisFunctionsRedEnv);
    auto nAuxBasFuncEnv = auxBasisControllerEnv->getNBasisFunctions();
    auto nAuxBasFuncRedEnv = auxBasisControllerEnv->getReducedNBasisFunctions();
    auto auxBasisIndicesRedEnv = auxBasisControllerEnv->getBasisIndicesRed();
    auto auxMappingEnv = createBasisToAtomIndexMapping(auxBasisIndicesRedEnv,nAuxBasFuncRedEnv);
    DensityMatrix<RESTRICTED> densityMatrixEnv = _envSystems[i]->template getElectronicStructure<SCFMode>()->getDensityMatrix().total();

    auto ri_j_IntController = RI_J_IntegralControllerFactory::getInstance()
    .produce(superSystemBasisController,superSystemAuxBasisController);

    // Get 2c Integrals (inverseM)
    auto inverseM = ri_j_IntController->getInverseM();


    /*
     * Constructing the different parts of eq. ([1].[6]).
     *
     * Starting with the C vectors. The first step combines the
     * Density Matrix elements ij with the respecting 3c Integral <J|ij>.
     * This is technically a Matrix-Vector multiplication.
     * Has to be done 4 times in order to get all combinations of
     * active and environment systems basis and auxBasis.
     *
     */
#ifdef _OPENMP
    // create a vector of matrices for each thread
    std::vector<Eigen::VectorXd> CvecActPriv(omp_get_max_threads(),Eigen::VectorXd::Zero(nAuxBasFunc));
    std::vector<Eigen::VectorXd> CvecEnvPriv(omp_get_max_threads(),Eigen::VectorXd::Zero(nAuxBasFuncEnv));
#else
    // or just one
    std::vector<Eigen::VectorXd> CvecActPriv(1,Eigen::VectorXd::Zero(nAuxBasFunc));
    std::vector<Eigen::VectorXd> CvecEnvPriv(1,Eigen::VectorXd::Zero(nAuxBasFuncEnv));
#endif

    TwoElecThreeCenterIntLooper looper1(libint2::Operator::coulomb,
        0,
        basisControllerEnv,
        auxBasisController,
        1E-10);

    auto const calcCVector1 = [&CvecActPriv,&densityMatrixEnv]
                               (const unsigned int& i,
                                   const unsigned int& j,
                                   const unsigned int& K,
                                   Eigen::VectorXd& intValues,
                                   const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecActPriv[threadId][K] += perm * intValues(0) * densityMatrixEnv(i,j);
    };
    looper1.loop(calcCVector1);

    TwoElecThreeCenterIntLooper looper2(libint2::Operator::coulomb,
        0,
        basisController,
        auxBasisController,
        1E-10);
    auto const calcCVector2 = [&CvecActPriv,&densityMatrix]
                               (const unsigned int& i,
                                   const unsigned int& j,
                                   const unsigned int& K,
                                   Eigen::VectorXd& intValues,
                                   const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecActPriv[threadId][K] += perm * intValues(0) * densityMatrix(i,j);
    };
    looper2.loop(calcCVector2);


    TwoElecThreeCenterIntLooper looper3(libint2::Operator::coulomb,
        0,
        basisController,
        auxBasisControllerEnv,
        1E-10);

    auto const calcCVector3 = [&CvecEnvPriv,&densityMatrix]
                               (const unsigned int& i,
                                   const unsigned int& j,
                                   const unsigned int& K,
                                   Eigen::VectorXd& intValues,
                                   const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecEnvPriv[threadId][K] += perm * intValues(0) * densityMatrix(i,j);
    };
    looper3.loop(calcCVector3);

    TwoElecThreeCenterIntLooper looper4(libint2::Operator::coulomb,
        0,
        basisControllerEnv,
        auxBasisControllerEnv,
        1E-10);
    auto const calcCVector4 = [&CvecEnvPriv,&densityMatrixEnv]
                               (const unsigned int& i,
                                   const unsigned int& j,
                                   const unsigned int& K,
                                   Eigen::VectorXd& intValues,
                                   const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecEnvPriv[threadId][K] += perm * intValues(0) * densityMatrixEnv(i,j);
    };
    looper4.loop(calcCVector4);
    /*
     * Has to be merged to a supersystem C^J because calculating the inverse of the
     * 2c-integral matrix in parts gives nonsensical results.
     */
    Eigen::VectorXd CvecAct(nAuxBasFunc);
    CvecAct.setZero();
    Eigen::VectorXd CvecEnv(nAuxBasFuncEnv);
    CvecEnv.setZero();
#ifdef _OPENMP
    for (unsigned int i =0;i<(unsigned int)omp_get_max_threads();++i){
      CvecAct += CvecActPriv[i];
      CvecEnv += CvecEnvPriv[i];
    }
#else
    CvecAct += CvecActPriv[0];
    CvecEnv += CvecEnvPriv[0];
#endif
    Eigen::VectorXd Cvec(CvecAct.rows()+CvecEnv.rows());
    Cvec << CvecAct,CvecEnv;
    Eigen::VectorXd Dvec = inverseM * Cvec;

    /*
     * The resulting Vector Dvec is the equivalent of the C^J vector in the paper.
     * The following set of loops calculates the complete 3rd part of eq ([1].[6])
     * (easily recognizable by the sign).
     *
     * Check CoulombPotential::getGeomGradients()!
     *
     * Has to be done twice in order to get the combinations
     * auxBasisEnv/Basis and auxBasis/BasisEnv.
     *
     */
    Libint& libint = Libint::getInstance();
    libint.initialize(libint2::Operator::coulomb,1,2);
#ifdef _OPENMP
    // create a vector of matrices for each thread
    std::vector<Eigen::MatrixXd> eriContrPriv(omp_get_max_threads(),Eigen::MatrixXd::Zero(nAtoms,3));
#else
    // or just one
    std::vector<Eigen::MatrixXd> eriContrPriv(1,Eigen::MatrixXd::Zero(nAtoms,3));
#endif

    for (unsigned int auxJ = 0; auxJ < auxBasisFunctions.size(); auxJ++){
      const unsigned int nAuxJ = auxBasisFunctions[auxJ]->getNContracted();
      unsigned int offauxJ = auxBasisController->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
      for (unsigned int auxK = 0; auxK < auxBasisFunctionsEnv.size(); auxK++){
        Eigen::MatrixXd intDerivs;
#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int nAuxK = auxBasisFunctionsEnv[auxK]->getNContracted();
        unsigned int offauxK = auxBasisControllerEnv->extendedIndex(auxK);
        if (libint.compute(libint2::Operator::coulomb,
            1,
            *auxBasisFunctions[auxJ],
            *auxBasisFunctionsEnv[auxK],
            intDerivs)){
          double perm = 1.0;
          Eigen::MatrixXd prefac = 0.5 * perm * Dvec.segment(offauxJ,nAuxJ) * Dvec.segment(offauxK+nAuxBasFunc,nAuxK).transpose();
          for (unsigned int iAtom = 0; iAtom < 2; ++iAtom) {
            unsigned int nAtom;
            switch(iAtom){
              case(0): nAtom = auxMapping[auxJ]; break;
              case(1): nAtom = auxMappingEnv[auxK]+atomsAct.size(); break;
            }
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom*3+iDirection).data(),nAuxK,nAuxJ);
              eriContrPriv[threadId](nAtom,iDirection) -= tmp.transpose().cwiseProduct(prefac).sum();
            }
          }
        }
      }
    }

    for (unsigned int auxJ = 0; auxJ < auxBasisFunctionsEnv.size(); auxJ++){
      const unsigned int nAuxJ = auxBasisFunctionsEnv[auxJ]->getNContracted();
      unsigned int offauxJ = auxBasisControllerEnv->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
      for (unsigned int auxK = 0; auxK < auxBasisFunctions.size(); auxK++){
        Eigen::MatrixXd intDerivs;
#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int nAuxK = auxBasisFunctions[auxK]->getNContracted();
        unsigned int offauxK = auxBasisController->extendedIndex(auxK);
        if (libint.compute(libint2::Operator::coulomb,
            1,
            *auxBasisFunctionsEnv[auxJ],
            *auxBasisFunctions[auxK],
            intDerivs)){
          double perm = 1.0;

          Eigen::MatrixXd prefac = 0.5 * perm * Dvec.segment(offauxJ+nAuxBasFunc,nAuxJ) * Dvec.segment(offauxK,nAuxK).transpose();
          for (unsigned int iAtom = 0; iAtom < 2; ++iAtom) {
            unsigned int nAtom;
            switch(iAtom){
              case(0): nAtom = auxMappingEnv[auxJ]+atomsAct.size(); break;
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

    for (unsigned int auxJ = 0; auxJ < auxBasisFunctions.size(); auxJ++){
      const unsigned int nAuxJ = auxBasisFunctions[auxJ]->getNContracted();
      unsigned int offauxJ = auxBasisController->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
      for (unsigned int auxK = 0; auxK < auxBasisFunctions.size(); auxK++){
        Eigen::MatrixXd intDerivs;
#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int nAuxK = auxBasisFunctions[auxK]->getNContracted();
        unsigned int offauxK = auxBasisController->extendedIndex(auxK);
        if (libint.compute(libint2::Operator::coulomb,
            1,
            *auxBasisFunctions[auxJ],
            *auxBasisFunctions[auxK],
            intDerivs)){
          double perm = 1.0;

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
     * The following set of loops calculates the complete 1st and 2nd part of eq ([1].[6])
     * (easily recognizable by the sign).
     *
     * Check CoulombPotential::getGeomGradients()!
     *
     * Has to be done twice in order to get the combinations
     * auxBasisEnv/Basis and auxBasis/BasisEnv.
     */

    TwoElecThreeCenterIntLooper derivLooper2(libint2::Operator::coulomb,
        1,
        basisControllerEnv,
        auxBasisController,
        1E-10);
    auto const add3cDerivs2 = [&eriContrPriv,&Dvec,&densityMatrixEnv,&mappingEnv,&auxMapping,
                               &basisControllerEnv,&auxBasisController,&nAtomsAct,&test]
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
            case(1): nAtom = mappingEnv[basisControllerEnv->reducedIndex(i)]+nAtomsAct; break;
            case(2): nAtom = mappingEnv[basisControllerEnv->reducedIndex(j)]+nAtomsAct; break;
          }
          eriContrPriv[threadId](nAtom,iDirection) += perm * Dvec[K] * densityMatrixEnv(i,j) * intValues(iAtom*3+iDirection);
        }
      }
    };
    derivLooper2.loop(add3cDerivs2);


    TwoElecThreeCenterIntLooper derivLooper1(libint2::Operator::coulomb,
        1,
        basisController,
        auxBasisControllerEnv,
        1E-10);
    auto const add3cDerivs1 = [&eriContrPriv,&Dvec,&densityMatrix,&mapping,&auxMappingEnv,
                               &basisController,&auxBasisControllerEnv,&nAtomsAct,&nAuxBasFunc,&test]
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
            case(0): nAtom = auxMappingEnv[auxBasisControllerEnv->reducedIndex(K)]+nAtomsAct; break;
            case(1): nAtom = mapping[basisController->reducedIndex(i)]; break;
            case(2): nAtom = mapping[basisController->reducedIndex(j)]; break;
          }
          eriContrPriv[threadId](nAtom,iDirection) += perm * Dvec[K+nAuxBasFunc] * densityMatrix(i,j) * intValues(iAtom*3+iDirection);
        }
      }
    };
    derivLooper1.loop(add3cDerivs1);


    TwoElecThreeCenterIntLooper derivLooper3(libint2::Operator::coulomb,
        1,
        basisController,
        auxBasisController,
        1E-10);
    auto const add3cDerivs3 = [&eriContrPriv,&Dvec,&densityMatrix,&mapping,&auxMapping,
                               &basisController,&auxBasisController,&nAtomsAct,&nAuxBasFunc,&test]
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
          eriContrPriv[threadId](nAtom,iDirection) += perm * Dvec[K] * densityMatrix(i,j) * intValues(iAtom*3+iDirection);
        }
      }
    };
    derivLooper3.loop(add3cDerivs3);


#ifdef _OPENMP
for (unsigned int i =0;i<(unsigned int)omp_get_max_threads();++i){
  eriContr += eriContrPriv[i];
}
#else
eriContr = eriContrPriv[0];
#endif

for (unsigned int i=0; i < atomsAct.size(); ++i){
  for (unsigned int j=0; j < 3; ++j){

    activeSystemGradientContr(i,j) += eriContr(i,j);
  }
}

  }
  // TODO: The above scheme also includes innersystem Coulomb gradient contributions once per environment system,
  // hence this has to be explicitely calculated and subtracted in the end. This is not incredibly expensive,
  // but unneccessary and, alas, very hard to extract from the above algorithm.
  CoulombPotential<SCFMode> coulPot(_actSystem,
      _actSystem->getElectronicStructure<SCFMode>()
      ->getDensityMatrixController(),
      RI_J_IntegralControllerFactory::getInstance().produce(
          _actSystem->getAtomCenteredBasisController(),
          _actSystem->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)),
          1E-10);
  auto ElEl = coulPot.getGeomGradients();
  return activeSystemGradientContr-(_envSystems.size()*ElEl);
}

template class CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED>;
template class CoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
