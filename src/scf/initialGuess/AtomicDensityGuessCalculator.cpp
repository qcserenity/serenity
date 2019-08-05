/**
 * @file AtomicDensityGuessCalculator.cpp
 * 
 * @date Jul 12, 2014
 * @author Thomas Dresselhaus, David Schnieders
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
#include "scf/initialGuess/AtomicDensityGuessCalculator.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/Basis.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrix.h"
#include "geometry/Geometry.h"
#include "integrals/wrappers/Libint.h"
#include "data/matrices/MatrixInBasis.h"
#include "integrals/OneElectronIntegralController.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "system/SystemController.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <vector>


namespace Serenity{
using namespace std;

std::unique_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED> > AtomicDensityGuessCalculator::calculateInitialDensity(
    shared_ptr<SystemController> systemController,
    bool keepMinimalBasis,
    bool scale) {
  assert(systemController);

  std::string pathToGuessDensities;
  if(const char* env_p = std::getenv("SERENITY_RESOURCES")){
    pathToGuessDensities = env_p;
  }else{
    std::cout << "ERROR Environment variable SERENITY_RESOURCES not set."<< std::endl;
    assert(false);
  }
  if(systemController->getSettings().basis.makeSphericalBasis){
    pathToGuessDensities += "initialGuess/spherical/";
  }else{
    pathToGuessDensities += "initialGuess/cartesian/";
  }
  if (_scf==GUESSMODES::SCF){
    if(systemController->getSettings().method==Options::ELECTRONIC_STRUCTURE_THEORIES::DFT){
      pathToGuessDensities += "dft_scf/";
    }else{
      pathToGuessDensities += "hf_scf/";
    }
  } else {
    pathToGuessDensities += "ao_occupation/";
  }

  /*
   * Get in some data locally
   */
  const auto& atoms = systemController->getAtoms();
  auto minimalBasisController =
      systemController->getAtomCenteredBasisController(
          (_scf==GUESSMODES::SCF)?Options::BASIS_PURPOSES::SCF_DENS_GUESS : Options::BASIS_PURPOSES::MINBAS);
  auto oneEIntC = systemController->getOneElectronIntegralController(
      (_scf==GUESSMODES::SCF)?Options::BASIS_PURPOSES::SCF_DENS_GUESS : Options::BASIS_PURPOSES::MINBAS);

  auto basisController=systemController->getBasisController();

  auto& libint = Libint::getInstance();
  libint.keepEngines(libint2::Operator::overlap,0,2);
  // these two variables are needed in the outer scope
  Eigen::MatrixXd targetBasisOverlap = systemController->getOneElectronIntegralController()->getOverlapIntegrals();

  auto newGuessDensity = unique_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED> >(new DensityMatrix<Options::SCF_MODES::RESTRICTED>(keepMinimalBasis?minimalBasisController:basisController));
  auto& newDens=*newGuessDensity;
  unsigned int blockstart=0;
  for(auto atom : atoms){
    if (atom->isDummy()){
      for (auto& shell : atom->getBasisFunctions()){
        blockstart+=shell->getNContracted();
      }
      continue;
    }
    /*
     * If not already available: Get from file
     */
    if(_atomDensities.find(atom->getAtomType()->getElementSymbol())==_atomDensities.end()){
      Eigen::MatrixXd atomDensMat;
      if(_scf==GUESSMODES::SCF_INPLACE){
        std::string pathRes=systemController->getSettings().path+atom->getAtomType()->getElementSymbol()+"_FREE/"+atom->getAtomType()->getElementSymbol()+"_FREE.dmat.res.h5";
        std::string pathUnres=systemController->getSettings().path+atom->getAtomType()->getElementSymbol()+"_FREE/"+atom->getAtomType()->getElementSymbol()+"_FREE.dmat.unres.h5";
        struct stat buffer;
        if(stat (pathRes.c_str(), &buffer) == 0){
          HDF5::Filepath name(pathRes);
          HDF5::H5File file(name.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
          HDF5::dataset_exists(file,"densityMatrix");
          HDF5::load(file,"densityMatrix",atomDensMat);
          file.close();
        }else if(stat (pathUnres.c_str(), &buffer) == 0){
          Eigen::MatrixXd dummy;
          HDF5::Filepath name(pathUnres);
          HDF5::H5File file(name.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
          HDF5::dataset_exists(file,"densityMatrix_alpha");
          HDF5::load(file,"densityMatrix_alpha",dummy);
          atomDensMat=dummy;
          HDF5::dataset_exists(file,"densityMatrix_beta");
          HDF5::load(file,"densityMatrix_beta",dummy);
          atomDensMat+=dummy;
          file.close();
        }else{
          // Copy system settings
          Settings atomSettings(systemController->getSettings());
          atomSettings.name=atom->getAtomType()->getElementSymbol()+"_FREE";
          // Uncharged atoms
          atomSettings.charge=0;
          if(atom->getNuclearCharge()%2!=0){
            atomSettings.spin=1;
          }else{
            atomSettings.spin=0;
          }
          atomSettings.scf.initialguess=Options::INITIAL_GUESSES::ATOM_SCF;
          atomSettings.scf.energyThreshold=1e-6;
          atomSettings.scf.rmsdThreshold=1e-6;
          atomSettings.scf.diisThreshold=1e-6;
          atomSettings.scf.damping=Options::DAMPING_ALGORITHMS::STATIC;
          atomSettings.scf.staticDampingFactor=0.7;
          atomSettings.scf.useOffDiagLevelshift=false;
          atomSettings.dft.dispersion=Options::DFT_DISPERSION_CORRECTIONS::NONE;
          std::vector<std::shared_ptr<Atom>> atomVec={atom};
          auto atomGeom=std::make_shared<Geometry>(atomVec);
          auto atomSys=std::make_shared<SystemController>(atomGeom,atomSettings);
          if(atomSettings.spin==1){
            ScfTask<UNRESTRICTED> scf(atomSys);
            scf.settings.fractionalDegeneracy=true;
            scf.run();
            atomDensMat=atomSys->getElectronicStructure<UNRESTRICTED>()->getDensityMatrix().total();
          }else{
            ScfTask<RESTRICTED> scf(atomSys);
            scf.settings.fractionalDegeneracy=true;
            scf.run();
            atomDensMat=atomSys->getElectronicStructure<RESTRICTED>()->getDensityMatrix();
          }
        }
      }else{
        HDF5::Filepath name(pathToGuessDensities+atom->getAtomType()->getElementSymbol()+".dmat.res.h5");
        HDF5::H5File file(name.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
        HDF5::dataset_exists(file,"densityMatrix");
        HDF5::load(file,"densityMatrix",atomDensMat);
        file.close();
      }
      /*
       * If basis change is wanted: Project to new basis
       */
      if (!keepMinimalBasis and _scf!=GUESSMODES::SCF_INPLACE){
        std::vector<std::shared_ptr<Atom>> dummyVec={atom};
        auto dummyGeom=std::make_shared<Geometry>(dummyVec);
        /*
         * Build Projection operation
         */
        auto atomMinBas=AtomCenteredBasisControllerFactory::produce(
            dummyGeom,
            systemController->getSettings().basis.basisLibPath,
            systemController->getSettings().basis.makeSphericalBasis,
            false,
            systemController->getSettings().basis.firstECP,
            (_scf==GUESSMODES::SCF)?"DEF2-QZVP":"STO-3G");
        auto atomTargetBas=AtomCenteredBasisControllerFactory::produce(
            dummyGeom,
            systemController->getSettings().basis.basisLibPath,
            systemController->getSettings().basis.makeSphericalBasis,
            false,
            systemController->getSettings().basis.firstECP,
            systemController->getSettings().basis.label);
        unsigned int nBasFunc=atomTargetBas->getNBasisFunctions();
        Eigen::MatrixXd overlapB=targetBasisOverlap.block(blockstart,blockstart,nBasFunc,nBasFunc);
        auto overlapAB = libint.compute1eInts(libint2::Operator::overlap,atomMinBas,atomTargetBas);
        //Calculate inverse of AO overlap integrals in basis B. Use SVD,
        //i.e. S_B = U * D * V^T, since overlapB could be ill-conditioned
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(overlapB,Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(1e-6);
        //Calculate projection operator P_{BA}
        Eigen::MatrixXd projectionOperator;
        projectionOperator = svd.solve(overlapAB);
        auto atomdens = (projectionOperator*atomDensMat*projectionOperator.transpose()).eval();
        const double factor = scale?atom->getEffectiveCharge()/overlapB.cwiseProduct(atomdens).sum():1.0;
        atomDensMat = atomdens*factor;
      }
      /*
       * Save matrix for this atom
       */
      _atomDensities[atom->getAtomType()->getElementSymbol()]=atomDensMat;

    }
    /*
     * Get atom density matrix and fill systems density matrix
     */
    auto atomDensMat=_atomDensities[atom->getAtomType()->getElementSymbol()];
    unsigned int nBasFunc=atomDensMat.rows();
    newDens.block(blockstart,blockstart,nBasFunc,nBasFunc)=atomDensMat;
    blockstart+=nBasFunc;
  }


  libint.freeEngines(libint2::Operator::overlap,0,2);

  return newGuessDensity;
}



} /* namespace Serenity */
