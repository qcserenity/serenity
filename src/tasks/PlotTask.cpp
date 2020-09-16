/**
 * @file   PlotTask.cpp
 *
 * @date   Nov 24, 2015
 * @author Jan Unsleber, Anja Massolle
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

/* Include Class Header*/
#include "tasks/PlotTask.h"
/* Include Serenity Internal Headers */
#include "analysis/localizationFunctions/ELFCalculator.h"     //ELF
#include "analysis/localizationFunctions/SEDD.h"              //SEDD
#include "data/ElectronicStructure.h"                         //ElectronicStructure
#include "data/OrbitalController.h"                           //OrbitalController
#include "data/grid/BasisFunctionOnGridControllerFactory.h"   //BasisFunctionOnGridControllerFactory
#include "data/grid/CoulombPotentialOnGridCalculator.h"       //CoulombPotentialOnGridCalculator
#include "data/grid/ElectrostaticPotentialOnGridController.h" //Plot electrostatic potential.
#include "geometry/Geometry.h"                                //Geometry construction.
#include "geometry/MolecularSurfaceController.h"              //Cavity grid getter.
#include "io/CubeFileWriter.h"                                //CubeFileWriter
#include "io/FormattedOutputStream.h"                         //OutputControl
#include "io/GeneralGridFileWriter.h"                         //Plot cavity.
#include "io/PlaneFileWriter.h"                               //PlaneFileWriter
#include "misc/WarningTracker.h"                              //WarningTracker
#include "settings/Settings.h"                                //Settings.
#include "system/SystemController.h"                          //System controller definition.
#include "tasks/SystemAdditionTask.h"                         //SystemAdditionTask

namespace Serenity {
using namespace std;

template<>
PlotTask<RESTRICTED>::PlotTask(const std::vector<std::shared_ptr<SystemController>>& systems,
                               const std::vector<std::shared_ptr<SystemController>>& environmentSystems, std::string filename)
  : _systems(systems), _environmentSystems(environmentSystems), _filename(filename), _fnameSuffix("_") {
}

template<>
PlotTask<UNRESTRICTED>::PlotTask(const std::vector<std::shared_ptr<SystemController>>& systems,
                                 const std::vector<std::shared_ptr<SystemController>>& environmentSystems, std::string filename)
  : _systems(systems), _environmentSystems(environmentSystems), _filename(filename), _fnameSuffix() {
  _fnameSuffix.alpha = "_alpha_";
  _fnameSuffix.beta = "_beta_";
}

template<Options::SCF_MODES SCFMode>
void PlotTask<SCFMode>::run() {
  std::vector<double> defaultCubeSpacing = {0.12, 0.12, 0.12};
  std::vector<double> defaultPlaneSpacing = {0.08, 0.08, 0.0};

  if (std::find(settings.p1.begin(), settings.p1.end(), std::numeric_limits<double>::infinity()) != settings.p1.end() and
      settings.atom1 == std::numeric_limits<int>::infinity()) {
    _cubeTask = true;
  }
  else if (std::find(settings.p2.begin(), settings.p2.end(), std::numeric_limits<double>::infinity()) != settings.p2.end() and
           settings.atom2 == std::numeric_limits<int>::infinity()) {
    throw SerenityError((string) "The plot for one point is not implemented yet!");
    _pointTask = true;
  }
  else if (std::find(settings.p3.begin(), settings.p3.end(), std::numeric_limits<double>::infinity()) != settings.p3.end() and
           settings.atom3 == std::numeric_limits<int>::infinity()) {
    throw SerenityError((string) "The plot for a linear grid is not implemented yet!");
    _lineTask = true;
  }
  else if (std::find(settings.p4.begin(), settings.p4.end(), std::numeric_limits<double>::infinity()) != settings.p4.end() and
           settings.atom4 == std::numeric_limits<int>::infinity()) {
    _planeTask = true;
  }
  else {
    throw SerenityError((string) "The plot for a cubic grid defined with four points is not implemented yet!");
    _cubeTask = true;
  }

  std::shared_ptr<SystemController> supersystem = nullptr;
  if (_systems.size() > 1) {
    Settings superSettings = _systems[0]->getSettings();
    superSettings.name = (_systems[0]->getSettings()).name + "+" + (_systems[1]->getSettings()).name;
    superSettings.path = "./";
    superSettings.charge = 0;
    superSettings.spin = 0;

    supersystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), superSettings);
    SystemAdditionTask<SCFMode> additionTask(supersystem, _systems);
    additionTask.settings.addOccupiedOrbitals = true;
    additionTask.run();
  }
  else {
    supersystem = _systems[0];
  }
  auto geom = supersystem->getGeometry();
  // Create supersystem geometry
  for (auto sys : _environmentSystems) {
    *geom += (*sys->getGeometry());
  }
  geom->deleteIdenticalAtoms();

  if (_filename.empty()) {
    _filename = supersystem->getSystemName();
  }

  Settings sysSettings = _systems[0]->getSettings();
  sysSettings.grid.blockAveThreshold = 0.0;
  sysSettings.grid.basFuncRadialThreshold = 0.0;
  // Null pointer to the base class DataOnGridWriter
  std::shared_ptr<DataOnGridWriter> plotWriter = nullptr;
  std::string fileExtension;
  if (settings.cavity) {
    printSubSectionTitle((string) "Printing Data on the Molecular Cavity");
    auto cavityGrid = supersystem->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE)->getGridController();
    plotWriter = std::make_shared<GeneralGridFileWriter>(sysSettings, cavityGrid);
  }
  else if (_cubeTask) {
    printSubSectionTitle((string) "Printing Data to Cube Files");
    if (isinf(settings.gridSpacing[0]) or isinf(settings.gridSpacing[1]) or isinf(settings.gridSpacing[2])) {
      settings.gridSpacing = defaultCubeSpacing;
    }
    plotWriter = std::make_shared<CubeFileWriter>(sysSettings, settings);
    fileExtension = ".cube";
  }
  else if (_planeTask) {
    printSubSectionTitle((string) "Printing Data to Dat Files");
    if (isinf(settings.gridSpacing[0]) or isinf(settings.gridSpacing[1]) or isinf(settings.gridSpacing[2])) {
      settings.gridSpacing = defaultPlaneSpacing;
    }
    plotWriter = std::make_shared<PlaneFileWriter>(sysSettings, settings);
    fileExtension = ".dat";
  }

  /*
   * print alpha/beta density
   */
  if (settings.density) {
    printSmallCaption("Densities");
    const MatrixInBasis<SCFMode>& density = supersystem->template getElectronicStructure<SCFMode>()->getDensityMatrix();
    auto basis = density.getBasisController();
    for_spin(density, _fnameSuffix) {
      print((string) "Printing electron density to file: " + _filename + _fnameSuffix_spin + "Density" + fileExtension);
      plotWriter->writeMatrixToGrid(supersystem->getSettings().path + _filename + _fnameSuffix_spin + "Density", geom,
                                    basis, density_spin);
    };
    if (SCFMode == UNRESTRICTED) {
      // print total density
      print((string) "Printing total electron density to file: " + _filename + "_TotalDensity" + fileExtension);
      plotWriter->writeMatrixToGrid(supersystem->getSettings().path + _filename + "_TotalDensity", geom, density.total());
      print((string) "Printing spin difference density to file: " + _filename + "_SpinDensity" + fileExtension);
      plotWriter->writeMatrixToGrid(supersystem->getSettings().path + _filename + "_SpinDensity", geom, density.difference());
    }
  }

  /*
   * Print Orbitals
   */
  if (settings.allOrbitals or settings.occOrbitals) {
    // set range for the lower loop
    auto range = supersystem->getNOccupiedOrbitals<SCFMode>();
    if (settings.allOrbitals) {
      for_spin(range) {
        range_spin = supersystem->getBasisController()->getNBasisFunctions();
      };
    }
    // print all requested orbitals
    const auto& coeffMat = supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
    printSmallCaption("Molecular Orbitals");
    print((string) "Printing MOs to files: " + _filename + "_<a/b>MO<number>" + fileExtension);
    if (_systems.size() > 1) {
      OutputControl::nOut << "The orbitals of the active systems are numbered as follows:" << std::endl;
      OutputControl::nOut << "  * all occupied orbitals of the first active system" << std::endl;
      OutputControl::nOut << "  * all occupied orbitals of the second active system etc." << std::endl;
      OutputControl::nOut << "If a SUBsystem basis is used the orbitals are further numbered as:" << std::endl;
      OutputControl::nOut << "  * all virtual orbitals of the first active system" << std::endl;
      OutputControl::nOut << "  * all virtual orbitals of the second active system etc." << std::endl;
      OutputControl::nOut << "If a SUPERsystem basis is used:" << std::endl;
      OutputControl::nOut << "  * all virtual orbitals are set to ZERO, if you want to print the virtual orbitals,"
                          << std::endl;
      OutputControl::nOut
          << "    print the orbitals for each active system seperately (and set the other active systems" << std::endl;
      OutputControl::nOut << "    as environment systems." << std::endl;
    }
    for_spin(coeffMat, range, _fnameSuffix) {
      Eigen::MatrixXd moCoefficients(coeffMat.rows(), range_spin);
      std::vector<std::string> fileNames;
      for (unsigned int i = 0; i < range_spin; i++) {
        moCoefficients.col(i) = coeffMat_spin.col(i);
        fileNames.push_back(supersystem->getSettings().path + _filename + _fnameSuffix_spin + "MO" + std::to_string(i + 1));
      } // for i
      plotWriter->writeVectorSetToGrid(fileNames, geom, supersystem->getBasisController(), moCoefficients);
    };
  }

  if (settings.orbitals.size() > 0) {
    printSmallCaption("Molecular Orbitals");
    print((string) "Printing MOs to files: " + _filename + "_<a/b>MO<number>" + fileExtension);
    if (_systems.size() > 1) {
      WarningTracker::printWarning(
          "\n\n The orbitals of the active systems are numbered as follows: \n   * all occupied orbitals of the first "
          "active system \n   * all occupied orbitals of the second active system etc. \n If a SUBsystem basis is used "
          "the orbitals are further numbered as:\n   * all virtual orbitals of the first active system \n   * all "
          "virtual orbitals of the second active system etc. \n If a SUPERsystem basis is used:\n   * all virtual "
          "orbitals are set to ZERO, if you want to print the virtual orbitals,\n     print the orbitals for each "
          "active system seperately (and set the other active systems\n     as environment systems).",
          true);
    }
    const auto& coeffMat = supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
    auto nOrbs = supersystem->getBasisController()->getNBasisFunctions();
    unsigned int nError = 0;
    for (unsigned int orbital : settings.orbitals) {
      if (orbital < (unsigned int)1 or orbital > (unsigned int)nOrbs) {
        std::cout << "Orbital " << orbital << " not within range!" << std::endl;
        ++nError;
        continue;
      }
    } // for orbital

    for_spin(coeffMat, _fnameSuffix) {
      Eigen::MatrixXd moCoefficients(nOrbs, settings.orbitals.size() - nError);
      std::vector<std::string> fileNames;
      unsigned int col = 0;
      for (unsigned int orbital : settings.orbitals) {
        if (orbital >= (unsigned int)1 and orbital <= (unsigned int)nOrbs) {
          moCoefficients.col(col) = coeffMat_spin.col(orbital - 1);
          fileNames.push_back(supersystem->getSettings().path + _filename + _fnameSuffix_spin + "MO" + std::to_string(orbital));
          ++col;
        }
      } // for orbital
      plotWriter->writeVectorSetToGrid(fileNames, geom, supersystem->getBasisController(), moCoefficients);
    };
  }

  if (settings.electrostaticPot) {
    printSmallCaption("Electrostatic Potential");
    // define lambda function for PlotWriter (using its
    // GridController)
    auto calculateESPOnGrid = [&](std::shared_ptr<GridController> gridController) {
      // some stuff needed
      auto densMat = supersystem->getElectronicStructure<SCFMode>()->getDensityMatrix();
      auto atoms = geom->getAtoms();
      GridPotential<RESTRICTED> potGrid(gridController);

      // add up potential induced by the electrons...
      CoulombPotentialOnGridCalculator::calculateElectronElectron<SCFMode>(potGrid, densMat);
      //...and the nuclei
      CoulombPotentialOnGridCalculator::calculateElectronNuclei(potGrid, atoms);
      /*
       * alpha and beta ESP will each contain the full ESP, but have
       * the wrong sign, since the CoulombPotentialOnGridCalculator
       * considers electrons rather than positive point charges
       */
      potGrid.array() *= -1.0;
      return potGrid;
    };
    plotWriter->writeFile(supersystem->getSettings().path + _filename + "_ESP", geom, calculateESPOnGrid);
  }
  if (settings.elf) {
    ELFCalculator<SCFMode> elf(supersystem);
    auto calculateElfOnGrid = [&](std::shared_ptr<GridController> gridController) {
      GridPotential<RESTRICTED> elfGrid(gridController);
      elfGrid = elf.calculateTotalELFOnGrid(gridController);
      return elfGrid;
    };
    plotWriter->writeFile(supersystem->getSettings().path + _filename + "_ELF", geom, calculateElfOnGrid);
  }
  if (settings.elfts) {
    ELFCalculator<SCFMode> elf(supersystem);
    auto calculateElfTSOnGrid = [&](std::shared_ptr<GridController> gridController) {
      GridPotential<RESTRICTED> elfGrid(gridController);
      elfGrid = elf.calculateTotalELFTSOnGrid(gridController);
      return elfGrid;
    };
    plotWriter->writeFile(supersystem->getSettings().path + _filename + "_ELF", geom, calculateElfTSOnGrid);
  }
  if (settings.sedd) {
    SEDD<SCFMode> sedd;
    auto lambda = sedd.getSEDDLambda(supersystem);
    plotWriter->writeFile(supersystem->getSettings().path + _filename + "_SEDD", geom, lambda);
  }
  if (settings.dori) {
    SEDD<SCFMode> dori;
    auto lambda = dori.getDORILambda(supersystem);
    plotWriter->writeFile(supersystem->getSettings().path + _filename + "_DORI", geom, lambda);
  }
  if (settings.signedDensity) {
    SEDD<SCFMode> dori;
    auto lambda = dori.getSignedDensityLambda(supersystem);
    plotWriter->writeFile(supersystem->getSettings().path + _filename + "_signedDensity", geom, lambda);
  }
  if (settings.gridCoordinates) {
    auto lambda = [&](std::shared_ptr<GridController> gridController) { return gridController->getWeights(); };
    plotWriter->writeFile(supersystem->getSettings().path + _filename + "_coords", geom, lambda);
  }
} /* run */

template class PlotTask<Options::SCF_MODES::RESTRICTED>;
template class PlotTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
