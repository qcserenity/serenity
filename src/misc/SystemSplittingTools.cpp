/**
 * @file SystemSplittingTools.cpp
 *
 * @date Sep 25, 2018
 * @author Moritz Bensberg
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
#include "misc/SystemSplittingTools.h"
/* Include Serenity Internal Headers */
#include "system/SystemController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "integrals/OneElectronIntegralController.h"
#include "basis/AtomCenteredBasisController.h"
#include "integrals/wrappers/Libint.h"

namespace Serenity {
template<Options::SCF_MODES SCFMode>
MatrixInBasis<SCFMode> SystemSplittingTools<SCFMode>::projectMatrixIntoNewBasis(
    const MatrixInBasis<SCFMode>& oldMatrix,
    std::shared_ptr<BasisController> newBasis,
    std::shared_ptr<MatrixInBasis<Options::SCF_MODES::RESTRICTED> > newOverlap) {
  auto& libint = Libint::getInstance();
  Eigen::MatrixXd overlapB;
  if(newOverlap) {
    overlapB = *newOverlap;
  } else {
    overlapB = libint.compute1eInts(libint2::Operator::overlap,oldMatrix.getBasisController(),oldMatrix.getBasisController());
  }
  Eigen::MatrixXd overlapAB = libint.compute1eInts(libint2::Operator::overlap,oldMatrix.getBasisController(),newBasis);
  //Calculate inverse of AO overlap integrals in basis B. Use SVD,
  //i.e. S_B = U * D * V^T, since overlapB could be ill-conditioned
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(overlapB,Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.setThreshold(1e-6);
  //Calculate projection operator P_{BA}
  Eigen::MatrixXd projectionOperator = svd.solve(overlapAB);
  MatrixInBasis<SCFMode> newMatrix(newBasis);
  for_spin(newMatrix,oldMatrix){
    newMatrix_spin = (projectionOperator*oldMatrix_spin*projectionOperator.transpose()).eval();
  };
  return newMatrix;
}
template<Options::SCF_MODES SCFMode>
unsigned int SystemSplittingTools<SCFMode>::matchAtom(
      std::shared_ptr<Geometry> geometry,
      std::shared_ptr<Atom> atom) {
  unsigned int actIndex = 0;
  for (const auto& actAtom : geometry->getAtoms()) {
    if(*actAtom == *makeAtomFromDummy(atom)){
      return actIndex;
    }
    ++actIndex;
  }
  ++actIndex;
  // if not found actIndex > nActAtoms
  return actIndex;
}

template<Options::SCF_MODES SCFMode>std::shared_ptr<DensityMatrix<SCFMode> >
SystemSplittingTools<SCFMode>::buildNonOrthogonalDensityMatrix(
    std::shared_ptr<SystemController> environmentSystem,
    SpinPolarizedData<SCFMode,std::vector<bool> > distantOrbitals) {
  const auto& envCoeff = environmentSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  double nOcc  = (SCFMode == Options::SCF_MODES::UNRESTRICTED) ? 1.0 : 2.0;
  auto nonOrthoDensMatPtr = std::make_shared<DensityMatrix<SCFMode> > (environmentSystem->getBasisController());
  auto& nonOrthoDensMat = *nonOrthoDensMatPtr;
  for_spin(envCoeff,distantOrbitals,nonOrthoDensMat) {
    nonOrthoDensMat_spin.setZero();
    unsigned int nNonOrtho = 0;
    for(const auto& nonOrtho : distantOrbitals_spin) if(nonOrtho) ++nNonOrtho;
    if(nNonOrtho > 0) {
      Eigen::MatrixXd nonOrthoCoeff(environmentSystem->getBasisController()->getNBasisFunctions(),nNonOrtho);
      unsigned int nAssigned = 0;
      for(unsigned int coeffIndex = 0; coeffIndex < distantOrbitals_spin.size(); ++coeffIndex) {
        if(distantOrbitals_spin[coeffIndex]) {
          nonOrthoCoeff.block(0,nAssigned,nonOrthoCoeff.rows(),1) = envCoeff_spin.col(coeffIndex);
          ++nAssigned;
        }//if distantOrbitals_spin[coeffIndex]
      }//for coeffIndex
      //Build density matrix as C*C^T * 2.0 or 1.0 (Restricted/Unrestricted)
      nonOrthoDensMat_spin = nOcc * (nonOrthoCoeff * nonOrthoCoeff.transpose());
    }//if nNonOrtho > 0
  };
  return nonOrthoDensMatPtr;
}

template<Options::SCF_MODES SCFMode>SpinPolarizedData<SCFMode,std::vector<bool> >
SystemSplittingTools<SCFMode>::selectDistantOrbitals(
    SPMatrix<SCFMode>& orbitalPopulations,
    std::shared_ptr<SystemController> activeSystem,
    std::shared_ptr<SystemController> environmentSystem,
    double basisFunctionRatio,
    double borderAtomThreshold) {
  SpinPolarizedData<SCFMode,std::vector<bool> >  distantOrbitals;
  unsigned int nEnvAtoms = environmentSystem->getGeometry()->getNAtoms();
  /*
   * 1. Match atoms and select distant atoms which have their number of shells reduced significantly.
   */
  auto actAtoms = activeSystem->getGeometry()->getAtoms();
  const auto& actSystemBasisLoadingData = activeSystem->getAtomCenteredBasisController()->getBasisLoadingData();

  std::vector<bool> distantAtoms(nEnvAtoms,true);
  unsigned int envAtomIndex = 0;
  for(const auto& envAtom : environmentSystem->getGeometry()->getAtoms()) {
    unsigned int matchedIndex = matchAtom(activeSystem->getGeometry(),envAtom);
    if(matchedIndex < actAtoms.size()) {
      //Check remaining basis function shells in the active system.
      double remainingBasisFunctionRatio = (double)actSystemBasisLoadingData[matchedIndex].shells.sum()
          /(double)actSystemBasisLoadingData[matchedIndex].shells.size();
      if(remainingBasisFunctionRatio >= basisFunctionRatio)distantAtoms[envAtomIndex] = false;
    }// if
    ++envAtomIndex;
  }// for envAtom
  /*
   * 2. Select the orbitals which are localized on not distant atom
   */
  // search for orbitals located on "distant" atoms
  // by getting the sum of Mulliken populations on not "distant" atoms.
  // If this sum is lower than the given threshold, the complete orbital is considered to be distant.
  const auto& nOccEnv = environmentSystem->getNOccupiedOrbitals<SCFMode>();
  unsigned int numberOfDistantOrbitals = 0;
  for_spin(distantOrbitals,nOccEnv,orbitalPopulations) {
    distantOrbitals_spin.resize(nOccEnv_spin,false);
    //Loop over all environment orbitals and atoms and evaluate the population on the non-distant atoms.
    //Sanity check for dimension of populations and environment atoms.
    assert((int)orbitalPopulations_spin.rows() == (int)distantAtoms.size()
        &&"Dimensions of atoms and their populations are not fitting!");
    for(unsigned int iEnvOrb = 0; iEnvOrb < nOccEnv_spin; ++iEnvOrb) {
      double popOnNotDistantAtoms = 0.0;
      for(unsigned int atomIndex = 0; atomIndex < distantAtoms.size(); ++atomIndex) {
        if(!distantAtoms[atomIndex]){
          popOnNotDistantAtoms += std::fabs(orbitalPopulations_spin(atomIndex,iEnvOrb));
        }//if !distantAtoms[atomIndex]
      }//for atomIndex
      if(popOnNotDistantAtoms < borderAtomThreshold) {
        distantOrbitals_spin[iEnvOrb] = true;
        ++numberOfDistantOrbitals;
      }
    }// for iEnvOrb
  };
  std::cout<<std::endl;
  std::cout << "Number of orbitals considered to be distant (alpha+beta): " << numberOfDistantOrbitals << std::endl;
  std::cout<<std::endl;
  return distantOrbitals;
}

template<Options::SCF_MODES SCFMode>
std::pair<std::shared_ptr<ElectronicStructure<SCFMode> >,std::shared_ptr<ElectronicStructure<SCFMode> > >
SystemSplittingTools<SCFMode>::splitElectronicStructure(
    std::shared_ptr<SystemController> system,
    SpinPolarizedData<SCFMode,std::vector<bool> >& activeOrbitals) {
  // get data
  const auto basisController = system->getBasisController();
  const unsigned int nBasisFunctions = basisController->getNBasisFunctions();
  const auto& superCoefficients = system->getActiveOrbitalController<SCFMode>()->getCoefficients();
  const auto& superEigenvalues  = system->getActiveOrbitalController<SCFMode>()->getEigenvalues();
  // Build coefficient matrix and eigenvalue vectors
  auto actCoeffPtr = std::unique_ptr<CoefficientMatrix<SCFMode> >(
      new CoefficientMatrix<SCFMode> (basisController));
  auto actEigenPtr = std::unique_ptr<SpinPolarizedData<SCFMode,Eigen::VectorXd> >(
      new SpinPolarizedData<SCFMode,Eigen::VectorXd> (nBasisFunctions));
  auto envCoeffPtr = std::unique_ptr<CoefficientMatrix<SCFMode> >(
      new CoefficientMatrix<SCFMode> (basisController));
  auto envEigenPtr = std::unique_ptr<SpinPolarizedData<SCFMode,Eigen::VectorXd> >(
      new SpinPolarizedData<SCFMode,Eigen::VectorXd> (nBasisFunctions));

  // Split coefficients
  auto& actCoefficients = *actCoeffPtr;
  auto& envCoefficients = *envCoeffPtr;
  for_spin(activeOrbitals,actCoefficients,superCoefficients,envCoefficients) {
    actCoefficients_spin.setZero();
    envCoefficients_spin.setZero();
    unsigned int aIndex = 0;
    unsigned int eIndex = 0;
    // occupied orbitals
    for(unsigned int i=0; i < activeOrbitals_spin.size(); ++i ) {
      if(activeOrbitals_spin[i]) {
        actCoefficients_spin.col(aIndex) = superCoefficients_spin.col(i);
        ++aIndex;
      } else {
        envCoefficients_spin.col(eIndex) = superCoefficients_spin.col(i);
        ++eIndex;
      }
    }
    // virtual orbital
    for(unsigned int i=activeOrbitals_spin.size(); i < nBasisFunctions; ++i) {
      actCoefficients_spin.col(aIndex) = superCoefficients_spin.col(i);
      envCoefficients_spin.col(eIndex) = superCoefficients_spin.col(i);
      ++aIndex;
      ++eIndex;
    }
  };
  // Split eigenvalues
  auto& actEigenvalues =  *actEigenPtr;
  auto& envEigenvalues =  *envEigenPtr;
  for_spin(activeOrbitals,superEigenvalues,actEigenvalues,envEigenvalues) {
    actEigenvalues_spin = Eigen::VectorXd::Constant(nBasisFunctions,std::numeric_limits<double>::infinity());
    envEigenvalues_spin = Eigen::VectorXd::Constant(nBasisFunctions,std::numeric_limits<double>::infinity());
    unsigned int aIndex = 0;
    for(const auto& act : activeOrbitals_spin) {
      if(act) aIndex++;
    }
    unsigned int eIndex = activeOrbitals_spin.size()-aIndex;
    // virtual orbital
    for(unsigned int i=activeOrbitals_spin.size(); i < nBasisFunctions; ++i) {
      actEigenvalues_spin(aIndex)= superEigenvalues_spin(i);
      envEigenvalues_spin(eIndex)= superEigenvalues_spin(i);
      ++aIndex;
      ++eIndex;
    }
  };
  // Number of occupied orbitals
  auto occupations=getNOccupiedOrbitals(activeOrbitals);
  // Build orbital controller
  auto activeOrbitalSet = std::make_shared<OrbitalController<SCFMode> >(
      std::move(actCoeffPtr),
      basisController,
      std::move(actEigenPtr));
  auto actES = std::make_shared<ElectronicStructure<SCFMode> >(
      activeOrbitalSet,
      system->getOneElectronIntegralController(),
      occupations.first);
  auto environmentOrbitalSet = std::make_shared<OrbitalController<SCFMode> >(
      std::move(envCoeffPtr),
      basisController,
      std::move(envEigenPtr));
  auto envES = std::make_shared<ElectronicStructure<SCFMode> >(
      environmentOrbitalSet,
      system->getOneElectronIntegralController(),
      occupations.second);
  //update density matrices
  actES->getDensityMatrixController()->updateDensityMatrix();
  envES->getDensityMatrixController()->updateDensityMatrix();
  assert(actES->getDensityMatrix().getBasisController()==basisController);
  assert(envES->getDensityMatrix().getBasisController()==basisController);
  return std::pair<std::shared_ptr<ElectronicStructure<SCFMode> >,std::shared_ptr<ElectronicStructure<SCFMode> > >(actES,envES);
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Atom> SystemSplittingTools<SCFMode>::makeAtomFromDummy(std::shared_ptr<Atom> atom) {
  auto atomType = atom->getAtomType();
  std::string nonDummyName =atomType->getName();
  if (atom->isDummy()) nonDummyName = nonDummyName.substr(0,nonDummyName.size()-1);
  auto newAtomType =std::make_shared<AtomType>(nonDummyName,
      atomType->getPSEPosition(),
      atomType->getMass(),
      atomType->getBraggSlaterRadius(),
      atomType->getVanDerWaalsRadius(),
      atomType->getOccupations());
  std::string basisLabel=atom->getPrimaryBasisLabel();
  auto basisFunctions=atom->getBasisFunctions();
  std::pair<std::string, std::vector<std::shared_ptr<Shell> >> atomBasisFunctions(basisLabel,basisFunctions);
  auto newAtom=std::make_shared<Atom>(
      newAtomType,
      atom->getX(),
      atom->getY(),
      atom->getZ(),
      atomBasisFunctions);
  return newAtom;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode,std::vector<bool> > SystemSplittingTools<SCFMode>::partitionOrbitals(
    SPMatrix<SCFMode>& orbitalPopulations,
    std::shared_ptr<SystemController> supersystem,
    std::vector<bool> activeAtoms,
    double orbitalThreshold,
    SpinPolarizedData<SCFMode,unsigned int> nOccOrbsActive,
    bool enforceCharges) {
  SpinPolarizedData<SCFMode,std::vector<bool> > listOfActiveOrbs;
  const unsigned int nAtoms = supersystem->getGeometry()->getAtoms().size();
  auto nOccOrbsSuper = supersystem->getNOccupiedOrbitals<SCFMode>();

  if (enforceCharges){
    // Pick n orbitals most localized on the active system.
    //   n is the number of orbitals defined by the number of electrons
    //   as generated from charge and spin in the original user input
    //   for the active system.
    for_spin(listOfActiveOrbs,nOccOrbsActive,nOccOrbsSuper,orbitalPopulations) {
      listOfActiveOrbs_spin.resize(nOccOrbsSuper_spin,false);
      Eigen::VectorXd popOnAct = Eigen::VectorXd::Zero(nOccOrbsSuper_spin);
      for (unsigned int i=0; i<nOccOrbsSuper_spin; ++i) {
        for (unsigned int a=0; a<nAtoms; ++a) {
          if (activeAtoms[a]) popOnAct[i] += fabs(orbitalPopulations_spin(a, i));
        }
      }
      for (unsigned int i=0; i<nOccOrbsActive_spin; ++i) {
        unsigned int m;
        popOnAct.maxCoeff(&m);
        listOfActiveOrbs_spin[m] = true;
        popOnAct[m] = 0.0;
      }
    };
  }else {
    // Pick all orbitals with a localization above the threshold
    //   given in the task settings.
    //   The number of electrons is then adjusted according to the
    //   number of electron that have been picked.
    for_spin(listOfActiveOrbs,nOccOrbsSuper,orbitalPopulations) {
      listOfActiveOrbs_spin.resize(nOccOrbsSuper_spin,false);

      for (unsigned int i=0; i<nOccOrbsSuper_spin; ++i) {
        double populationOnAllAtoms = 0.0;
        for (unsigned int a=0; a<nAtoms; ++a) {
          if (activeAtoms[a]) populationOnAllAtoms += fabs(orbitalPopulations_spin(a, i));
        }
        if (populationOnAllAtoms >= orbitalThreshold) {
          listOfActiveOrbs_spin[i] = true;
        }
      }
    };
  }
  return listOfActiveOrbs;
}

template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int>>
SystemSplittingTools<SCFMode>::getNOccupiedOrbitals(
    SpinPolarizedData<SCFMode,std::vector<bool> >& activeOrbitals) {
  std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> > occupations(
      SpinPolarizedData<SCFMode,unsigned int>(0),SpinPolarizedData<SCFMode,unsigned int>(0));
  auto& actOcc = occupations.first;
  auto& envOcc = occupations.second;
  for_spin(activeOrbitals,actOcc,envOcc) {
    for(const auto& orbBool : activeOrbitals_spin) {
      if (orbBool) {
        ++actOcc_spin;
      } else {
        ++envOcc_spin;
      }
    }
  };
  return occupations;
}

template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> >
SystemSplittingTools<SCFMode>::getNElectrons(
    std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> > occ) {
  auto nOrbsAct = occ.first;
  auto nOrbsEnv = occ.second;
  SpinPolarizedData<SCFMode,unsigned int> actElectrons = 0;
  SpinPolarizedData<SCFMode,unsigned int> envElectrons = 0;
  unsigned int nElectronsPerOrbital = (SCFMode == Options::SCF_MODES::RESTRICTED)? 2 : 1;
  for_spin(nOrbsAct,nOrbsEnv,actElectrons,envElectrons) {
    actElectrons_spin += nElectronsPerOrbital * nOrbsAct_spin;
    envElectrons_spin += nElectronsPerOrbital * nOrbsEnv_spin;
  };
  return std::pair<SpinPolarizedData<SCFMode,unsigned int>,SpinPolarizedData<SCFMode,unsigned int> > (actElectrons,envElectrons);
}

template<>
int SystemSplittingTools<Options::SCF_MODES::RESTRICTED>::getSpin(
    SpinPolarizedData<Options::SCF_MODES::RESTRICTED,unsigned int> nOcc) {
  (void) nOcc;
  return 0;
}
template<>
int SystemSplittingTools<Options::SCF_MODES::UNRESTRICTED>::getSpin(
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED,unsigned int> nOcc) {
  return (int)nOcc.alpha - (int)nOcc.beta;
}

template<Options::SCF_MODES SCFMode> std::vector<std::shared_ptr<Eigen::MatrixXd> >
SystemSplittingTools<SCFMode>::getProjectedSubsystems(
    std::shared_ptr<SystemController> activeSystem,
    std::vector<std::shared_ptr<SystemController> > environmentSystems,
    double truncationThreshold) {
  std::vector<std::shared_ptr<Eigen::MatrixXd> > overlapMatrices;
  auto basisContA = activeSystem->getBasisController();
  for(auto env : environmentSystems) {
    auto& libint = Libint::getInstance();
    auto s_AB = std::make_shared<Eigen::MatrixXd> (
        libint.compute1eInts(libint2::Operator::overlap, env->getBasisController(),basisContA));
    double totalInterSystemBasisOverlap = 0.0;
    for(unsigned int k = 0; k < s_AB->rows(); ++k) {
      for (unsigned int l = 0; l < s_AB->cols(); ++l) {
        totalInterSystemBasisOverlap += std::abs((*s_AB)(k,l));
      }// for l
    }// for k
    if (totalInterSystemBasisOverlap <= truncationThreshold) {
      // If the system does not get projected, forget the basis overlap.
      s_AB = nullptr;
    }
    overlapMatrices.push_back(s_AB);
  }// for env
  return overlapMatrices;
}

template<Options::SCF_MODES SCFMode> std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > >
SystemSplittingTools<SCFMode>::getEnvironmentDensityControllers(
    std::vector<std::shared_ptr<SystemController> > environmentSystems,
    bool topDown) {
  // Build the density matrix controllers of the environment systems in the same
  // spin-polarization as the active system.
  // These are needed for the Coulomb and exchange contribution to the supersystem fock operator of systems,
  // which are not in the system pair.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > environmentDensityControllers;
  for (auto sys : environmentSystems){
    if (sys->getSettings().scfMode == SCFMode || topDown) {
      assert(sys->hasElectronicStructure<SCFMode>());
      environmentDensityControllers.push_back(sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());
    } else {
      if (sys->getSettings().scfMode == Options::SCF_MODES::RESTRICTED) {
        //Build unrestricted DensityMatrixController
        DensityMatrix<SCFMode> uDensMat(sys->getBasisController());
        for_spin(uDensMat) {
          uDensMat_spin = 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        };
        environmentDensityControllers.push_back(std::make_shared<DensityMatrixController<SCFMode>>(uDensMat));
      } else if (sys->getSettings().scfMode == Options::SCF_MODES::UNRESTRICTED) {
        //Build restricted DensityMatrixController
        DensityMatrix<SCFMode> rDensMat(sys->getBasisController());
        for_spin(rDensMat) {
          rDensMat_spin = sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
        };
        environmentDensityControllers.push_back(std::make_shared<DensityMatrixController<SCFMode>>(rDensMat));
      } else {
        assert(false);
      }
    }
  }// for sys
  return environmentDensityControllers;
}


template class SystemSplittingTools<Options::SCF_MODES::RESTRICTED>;
template class SystemSplittingTools<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
