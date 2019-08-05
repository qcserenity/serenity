/**
 * @file NEInteractionPotential.cpp
 *
 * @date Nov 24, 2016
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
#include "potentials/NEInteractionPotential.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "integrals/wrappers/Libint.h"
#include "misc/Timing.h"


namespace Serenity {

template <Options::SCF_MODES SCFMode>
NEInteractionPotential<SCFMode>::NEInteractionPotential(
    std::shared_ptr<SystemController> activeSystem,
    std::vector<std::shared_ptr<SystemController> > environmentSystems,
    std::shared_ptr<BasisController> basis,
    std::vector<std::shared_ptr<const Geometry> > geometries):
    Potential<SCFMode>(basis),
    _actSystem(activeSystem),
    _envSystems(environmentSystems),
    _geometries(geometries){
  for (auto& geo : _geometries){
    for (auto& atom : geo->getAtoms()){
      atom->addSensitiveObject(ObjectSensitiveClass<Atom>::_self);
    }
  }
}


template <Options::SCF_MODES SCFMode> FockMatrix<SCFMode>&
NEInteractionPotential<SCFMode>::getMatrix(){
  Timings::takeTime("FDE -         1e-Int Pot.");
  if (!_potential){
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot){
      pot_spin.setZero();
    };

    auto& libint = Libint::getInstance();

    Geometry totalEnvGeometry;
    for (auto subsystemGeometry : _geometries) {
      totalEnvGeometry +=  *subsystemGeometry;
    }
    auto nucInts = libint.compute1eInts(libint2::Operator::nuclear,this->_basis,totalEnvGeometry.getAtoms());

    for_spin(pot){
      pot_spin += nucInts;
    };

  }
  Timings::timeTaken("FDE -         1e-Int Pot.");
  return *_potential;
}

template <Options::SCF_MODES SCFMode> double
NEInteractionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P){
  if (!_potential) this->getMatrix();
  Timings::takeTime("FDE -         1e-Int Pot.");
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot,P){
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("FDE -         1e-Int Pot.");
  return energy;
};

template<Options::SCF_MODES T> std::vector<unsigned int>
NEInteractionPotential<T>::createBasisToAtomIndexMapping(
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
      throw SerenityError("NEInteractionPotential: Missed gradient element in gradient evaluation.");
  } 
  return mapping;
}

template <Options::SCF_MODES SCFMode> Eigen::MatrixXd
NEInteractionPotential<SCFMode>::getGeomGradients(){
  /*
     * Supersystem geometry is needed within a single vector for libint to calculate all Nuc-El terms
     */
  auto atomsAct = _actSystem->getAtoms();
  Matrix<double> activeSystemGradientContr(atomsAct.size(),3);
  activeSystemGradientContr.setZero();

  for (unsigned int env=0;env<_envSystems.size();++env){
    auto atomsEnv = _envSystems[env]->getAtoms();
    std::vector<std::shared_ptr<Atom> > ActiveFrozenAtoms;
    for (unsigned int i = 0; i < atomsAct.size(); i++){
      ActiveFrozenAtoms.push_back(atomsAct[i]);
    }
    for (unsigned int i = 0; i < atomsEnv.size(); i++){
      ActiveFrozenAtoms.push_back(atomsEnv[i]);
    }
    unsigned int nAtoms = ActiveFrozenAtoms.size();
    Matrix<double> gradientContr(nAtoms,3);
    gradientContr.setZero();
    /*
     * Both density matrices are needed for the prefactors. Cross-system density matrix elements evaluate to 0,
     * which means there is no need for a supersystem density matrix.
     */
    DensityMatrix<RESTRICTED> matrix(_actSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix().total());
    DensityMatrix<RESTRICTED> matrixEnv(_envSystems[env]->template getElectronicStructure<SCFMode>()->getDensityMatrix().total());

    /*
     * Mapping basis functions of both active and environment system to their respecting atoms.
     * Both mapping vectors are handled seperately, which means they both start at atom 0 of their
     * respecting system.
     */
    auto atomCenteredBasisController = _actSystem->getAtomCenteredBasisController();
    auto atomCenteredBasisControllerEnv = _envSystems[env]->getAtomCenteredBasisController();
    unsigned int nBasisFunctionsRed = atomCenteredBasisController->getReducedNBasisFunctions();
    unsigned int nBasisFunctionsRedEnv = atomCenteredBasisControllerEnv->getReducedNBasisFunctions();
    auto basisIndicesRed = _actSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
    auto basisIndicesRedEnv = _envSystems[env]->getAtomCenteredBasisController()->getBasisIndicesRed();
    auto mapping = createBasisToAtomIndexMapping(basisIndicesRed,nBasisFunctionsRed);
    auto mappingEnv = createBasisToAtomIndexMapping(basisIndicesRedEnv,nBasisFunctionsRedEnv);

    auto& basis = atomCenteredBasisController->getBasis();
    auto& basisEnv = atomCenteredBasisControllerEnv->getBasis();

    /*
     * Since the following first loop will also include the Nuc-El contributions of basis functions of the active system
     * with atoms of the active system (a contribution that is already handled in the active system gradient), it has to
     * be subtracted in the end. This is extra work and there are certainly ways around having to do this.
     */
    /*
     * First loop calculates all Nuc-El terms involving the active systems basis functions
     */
    Libint& libint = Libint::getInstance();
    libint.initialize(libint2::Operator::nuclear,1,2,atomsEnv);

  #pragma omp parallel
    {
      //  std::vector<double> intDerivs;
      Eigen::MatrixXd intDerivs;
      Matrix<double> gradientContrPriv(nAtoms,3);
      gradientContrPriv.setZero();
      bool significant;
  #pragma omp for schedule(static,1)
      for (unsigned int i = 0; i < basis.size(); ++i) {
        unsigned int offI = _actSystem->getBasisController()->extendedIndex(i);
        const unsigned int nI = basis[i]->getNContracted();

        for (unsigned int j = 0; j <= i; ++j) {
          unsigned int offJ = _actSystem->getBasisController()->extendedIndex(j);
          const unsigned int nJ = basis[j]->getNContracted();

          /*
           * Libint needs an extra list of all environment system atoms for Nuc-El terms
           */
          significant = libint.compute(libint2::Operator::nuclear,
              1,
              *basis[i],
              *basis[j],
              intDerivs);

          if ( significant){
            double perm = (i==j ? 1.0 : 2.0);
            /*
             * Mapping contributions to each respecting atom. Since Active-Active contributions are unnecessary here,
             * simply mapping the Active-Environment contributions to the two atoms that carry the respecitve basis functions i
             * and j is sufficient.
             */
            for (unsigned int iAtom = 0; iAtom < 2; ++iAtom) {
              unsigned int nAtom;
              switch (iAtom){
                case(0): nAtom=mapping[i]; break;
                case(1): nAtom=mapping[j]; break;
              }
              for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
                Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom*3+iDirection).data(),nJ,nI);

                gradientContrPriv(nAtom,iDirection) += perm * tmp.cwiseProduct(matrix.block(offJ,offI,nJ,nI)).sum();
              }
            }
          }
        }
      }
  #pragma omp critical
      {
        gradientContr += gradientContrPriv;
      }
    } /* END OpenMP parallel */
    libint.finalize(libint2::Operator::nuclear,1,2);
    /*
     * Second loop calculates all Nuc-El terms involving the environment systems basis functions
     */
    libint.initialize(libint2::Operator::nuclear,1,2,atomsAct);
  #pragma omp parallel
    {
      //  std::vector<double> intDerivs;
      Eigen::MatrixXd intDerivs;
      Matrix<double> gradientContrPriv(nAtoms,3);
      gradientContrPriv.setZero();
      bool significant;
  #pragma omp for schedule(static,1)
      for (unsigned int i = 0; i < basisEnv.size(); ++i) {
        unsigned int offI = _envSystems[env]->getBasisController()->extendedIndex(i);
        const unsigned int nI = basisEnv[i]->getNContracted();
        for (unsigned int j = 0; j <= i; ++j) {
          unsigned int offJ = _envSystems[env]->getBasisController()->extendedIndex(j);
          const unsigned int nJ = basisEnv[j]->getNContracted();

          significant = libint.compute(libint2::Operator::nuclear,
              1,
              *basisEnv[i],
              *basisEnv[j],
              intDerivs);

          if ( significant){
            double perm = (i==j ? 1.0 : 2.0);

            /*
             * Only terms between basisfunctions of Env and atoms of Act are needed here. The loop below is
             * altered to not add the contributions that arise at the coordinates of the basis functions (case 0 and 1)
             * to anything since they are purely located in the environment and thus not needed.
             */

            for (unsigned int iAtom = 0; iAtom < atomsAct.size()+2; ++iAtom) {
              unsigned int nAtom;
              switch (iAtom){
                case(0): break;
                case(1): break;
                default: nAtom=iAtom-2;
                for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {

                  Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom*3+iDirection).data(),nJ,nI);
                    gradientContrPriv(nAtom,iDirection) += perm * tmp.cwiseProduct(matrixEnv.block(offJ,offI,nJ,nI)).sum();
                } break; }
            }
          }
        }
      }
  #pragma omp critical
      {
        gradientContr += gradientContrPriv;
      }
    } /* END OpenMP parallel */
    libint.finalize(libint2::Operator::nuclear,1,2);
    for (unsigned int i=0; i < atomsAct.size(); ++i){
      for (unsigned int j=0; j < 3; ++j){
        activeSystemGradientContr(i,j) += gradientContr(i,j);
      }
    }
  }
    return activeSystemGradientContr;
}

template class NEInteractionPotential<Options::SCF_MODES::RESTRICTED>;
template class NEInteractionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
