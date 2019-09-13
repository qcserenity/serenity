/**
 * @file   CubeFileTask.cpp
 *
 * @date   Nov 24, 2015
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
#include "tasks/CubeFileTask.h"
/* Include Serenity Internal Headers */
#include "geometry/Atom.h"
#include "basis/BasisController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/grid/CoulombPotentialOnGridCalculator.h"
#include "io/CubeFileWriter.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/ElectronicStructure.h"
#include "io/FormattedOutput.h"
#include "grid/Grid.h"
#include "data/grid/GridData.h"
#include "data/grid/MOCalculator.h"
//#include "postHF/LRSCF/Analysis/NTOCalculator.h"
#include "data/OrbitalController.h"
#include "analysis/localizationFunctions/SEDD.h"
#include "analysis/localizationFunctions/ELFCalculator.h"
#include "settings/Settings.h"
#include "data/SpinPolarizedData.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <iostream>
#include <string>


namespace Serenity {
using namespace std;

template <>
CubeFileTask<RESTRICTED>::CubeFileTask(
    const std::vector<std::shared_ptr<SystemController> >& systems,
    const std::vector<std::shared_ptr<SystemController> >& environmentSystems):
    _systems(systems),
    _environmentSystems(environmentSystems),
    _fnameSuffix("_"){
}

template <>
CubeFileTask<UNRESTRICTED>::CubeFileTask(
    const std::vector<std::shared_ptr<SystemController> >& systems,
    const std::vector<std::shared_ptr<SystemController> >& environmentSystems):
    _systems(systems),
    _environmentSystems(environmentSystems),
    _fnameSuffix(){
  _fnameSuffix.alpha = "_alpha_";
  _fnameSuffix.beta = "_beta_";
}

template <Options::SCF_MODES SCFMode>
void CubeFileTask<SCFMode>::run() {
  printSubSectionTitle((string)"Printing Data to Cube Files" );

  Settings sysSettings=_systems[0]->getSettings();
  sysSettings.grid.blockAveThreshold=0.0;
  sysSettings.grid.basFuncRadialThreshold=0.0;
  CubeFileWriter cubeWriter(sysSettings);
  cubeWriter.setBorderWidth(settings.cubeBorder);
  cubeWriter.setStepSize(Eigen::Vector3d::Ones()*settings.cubeSpacing);
  std::shared_ptr<SystemController> supersystem = _systems[0];
  for (unsigned int i=1;i<_systems.size();i++) {
    supersystem = (*supersystem)+(*_systems[i]);
  }
  auto geom = supersystem->getGeometry();
  //Create supersystem geometry
  for (auto sys : _environmentSystems) {
    *geom += (*sys->getGeometry());
  }
  geom->deleteIdenticalAtoms();
  /*
   * print alpha/beta density
   */
  if(settings.density){
    printSmallCaption("Densities");
    const MatrixInBasis<SCFMode>& density = supersystem->template getElectronicStructure<SCFMode>()->getDensityMatrix();
    auto basis = density.getBasisController();
    for_spin(density,_fnameSuffix){
      print((string)"Printing electron density to file: "+supersystem->getSystemName()+_fnameSuffix_spin+"Density.cube");
      cubeWriter.writeMatrixToCube(
          supersystem->getSettings().path+(string)supersystem->getSystemName()+_fnameSuffix_spin+"Density",
          geom,
          basis,
          density_spin);
    };
    if (SCFMode == UNRESTRICTED){
      // print total density
      print((string)"Printing total electron density to file: "+supersystem->getSystemName()+"_TotalDensity.cube");
      cubeWriter.writeMatrixToCube(supersystem->getSettings().path+(string)supersystem->getSystemName()+"_TotalDensity",
          geom,
          density.total());
      print((string)"Printing spin difference density to file: "+supersystem->getSystemName()+"_SpinDensity.cube");
      cubeWriter.writeMatrixToCube(supersystem->getSettings().path+(string)supersystem->getSystemName()+"_SpinDensity",
          geom,
          density.difference());
    }
  }

  /*
   * Print Orbitals
   */
  if(settings.allOrbitals or settings.occOrbitals){
    // set range for the lower loop
    auto range = supersystem->getNOccupiedOrbitals<SCFMode>();
    if (settings.allOrbitals) {
      for_spin(range){
        range_spin =  supersystem->getBasisController()->getNBasisFunctions();
      };
    }
    // print all requested orbitals
    const auto& coeffMat = supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
    printSmallCaption("Molecular Orbitals");
    print((string)"Printing MOs to files: "+supersystem->getSystemName()+"_<a/b>MO<number>.cube");
    for_spin(coeffMat,range,_fnameSuffix){
      Eigen::MatrixXd moCoefficients(coeffMat.rows(),range_spin);
      std::vector<std::string> fileNames;
      for (unsigned int i=0;i<range_spin;i++){
        moCoefficients.col(i) = coeffMat_spin.col(i);
        fileNames.push_back(
            supersystem->getSettings().path+(string)supersystem->getSystemName()+_fnameSuffix_spin+"MO"+std::to_string(i+1));
      }// for i
      cubeWriter.writeVectorSetToCube(
          fileNames,
          geom,
          supersystem->getBasisController(),
          moCoefficients);
    };
  }

  if(settings.orbitals.size() > 0){
    printSmallCaption("Molecular Orbitals");
    print((string)"Printing MOs to files: "+supersystem->getSystemName()+"_<a/b>MO<number>.cube");
    const auto& coeffMat = supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
    auto nOrbs=supersystem->getBasisController()->getNBasisFunctions();
    unsigned int nError = 0;
    for(unsigned int orbital : settings.orbitals){
      if(orbital<(unsigned int)1 or orbital>(unsigned int)nOrbs){
        std::cout << "Orbital " << orbital << " not within range!" << std::endl;
        ++nError;
        continue;
      }
    }// for orbital

    for_spin(coeffMat,_fnameSuffix) {
      Eigen::MatrixXd moCoefficients(nOrbs,settings.orbitals.size()-nError);
      std::vector<std::string> fileNames;
      unsigned int col=0;
      for(unsigned int orbital : settings.orbitals){
        if(orbital>=(unsigned int)1 and orbital<=(unsigned int)nOrbs) {
          moCoefficients.col(col) = coeffMat_spin.col(orbital-1);
          fileNames.push_back(
              supersystem->getSettings().path+(string)supersystem->getSystemName()+_fnameSuffix_spin+"MO"+std::to_string(orbital));
          ++col;
        }
      }// for orbital
      cubeWriter.writeVectorSetToCube(
          fileNames,
          geom,
          supersystem->getBasisController(),
          moCoefficients);
    };
  }

  if(settings.electrostaticPot){
    printSmallCaption("Electrostatic Potential");
    //define lambda function for CubeFileWriter (using its GridController)
    auto calculateESPOnGrid = [&]( std::shared_ptr<GridController> gridController ){
      //some stuff needed
      auto densMat=supersystem->getElectronicStructure<SCFMode>()->getDensityMatrix();
      auto atoms=geom->getAtoms();
      GridPotential<RESTRICTED> potGrid(gridController);
      auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
          supersystem->getSettings(),
          supersystem->getBasisController(),
          gridController);
      CoulombPotentialOnGridCalculator coulOnGridCalc(basFuncOnGridController);

      //add up potential induced by the electrons...
      coulOnGridCalc.calculateElectronElectron<SCFMode>(potGrid,densMat);
      //...and the nuclei
      coulOnGridCalc.calculateElectronNuclei(potGrid,atoms);
      /*
       * alpha and beta ESP will each contain the full ESP, but have the wrong sign, since
       * the CoulombPotentialOnGridCalculator considers electrons rather than positive point
       * charges
       */
      potGrid.array() *= -1.0;
      return potGrid;
    };
    cubeWriter.writeCube(supersystem->getSettings().path+supersystem->getSystemName() + "_ESP",geom,calculateESPOnGrid);
  }
  if (settings.elf){
    ELFCalculator<SCFMode> elf(supersystem);
    auto calculateElfOnGrid = [&](std::shared_ptr<GridController> gridController){
      GridPotential<RESTRICTED> elfGrid(gridController);
      elfGrid = elf.calculateTotalELFOnGrid(gridController);
      return elfGrid;
    };
    cubeWriter.writeCube(supersystem->getSettings().path+supersystem->getSystemName() + "_ELF",geom,
        calculateElfOnGrid);
  }
  if (settings.elfts){
    ELFCalculator<SCFMode> elf(supersystem);
    auto calculateElfTSOnGrid = [&](std::shared_ptr<GridController> gridController){
      GridPotential<RESTRICTED> elfGrid(gridController);
      elfGrid = elf.calculateTotalELFTSOnGrid(gridController);
      return elfGrid;
    };
    cubeWriter.writeCube(supersystem->getSettings().path+supersystem->getSystemName() + "_ELF",geom,
        calculateElfTSOnGrid);
  }
  if (settings.sedd){
    SEDD<SCFMode> sedd;
    auto lambda = sedd.getSEDDLambda(supersystem);
    cubeWriter.writeCube(supersystem->getSettings().path+supersystem->getSystemName() + "_SEDD",geom,
        lambda);
  }
  if (settings.dori){
    SEDD<SCFMode> dori;
    auto lambda = dori.getDORILambda(supersystem);
    cubeWriter.writeCube(supersystem->getSettings().path+supersystem->getSystemName() + "_DORI",geom,
        lambda);
  }
  if (settings.signedDensity){
    SEDD<SCFMode> dori;
    auto lambda = dori.getSignedDensityLambda(supersystem);
    cubeWriter.writeCube(supersystem->getSettings().path+supersystem->getSystemName() + "_signedDensity",geom,
        lambda);
  }

}

template class CubeFileTask<Options::SCF_MODES::RESTRICTED>;
template class CubeFileTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
