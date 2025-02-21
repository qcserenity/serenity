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
#include "analysis/localizationFunctions/ELFCalculator.h" //ELF
#include "analysis/localizationFunctions/SEDD.h"          //SEDD
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"                         //ElectronicStructure
#include "data/OrbitalController.h"                           //OrbitalController
#include "data/grid/BasisFunctionOnGridControllerFactory.h"   //BasisFunctionOnGridControllerFactory
#include "data/grid/CoulombPotentialOnGridCalculator.h"       //CoulombPotentialOnGridCalculator
#include "data/grid/ElectrostaticPotentialOnGridController.h" //Plot electrostatic potential.
#include "geometry/Geometry.h"                                //Geometry construction.
#include "geometry/MolecularSurfaceController.h"              //Cavity grid getter.
#include "integrals/OneElectronIntegralController.h"
#include "io/CubeFileWriter.h"                 //CubeFileWriter
#include "io/FormattedOutputStream.h"          //OutputControl
#include "io/GeneralGridFileWriter.h"          //Plot cavity.
#include "io/MolecularSurfaceGridFileWriter.h" //Plot cavity
#include "io/PlaneFileWriter.h"                //PlaneFileWriter
#include "misc/WarningTracker.h"               //WarningTracker
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
#include "postHF/LRSCF/Analysis/NROCalculator.h"
#include "postHF/LRSCF/Analysis/NTOCalculator.h" //NTO
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"       //Settings.
#include "system/SystemController.h" //System controller definition.
#include "tasks/LRSCFTask.h"
#include "tasks/SystemAdditionTask.h" //SystemAdditionTask

namespace Serenity {

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
  this->avoidMixedSCFModes(SCFMode, _systems);
  std::vector<double> defaultCubeSpacing = {0.12, 0.12, 0.12};
  std::vector<double> defaultPlaneSpacing = {0.08, 0.08, 0.0};

  if (std::find(settings.p1.begin(), settings.p1.end(), std::numeric_limits<double>::infinity()) != settings.p1.end() and
      settings.atom1 == std::numeric_limits<int>::infinity()) {
    _cubeTask = true;
  }
  else if (std::find(settings.p2.begin(), settings.p2.end(), std::numeric_limits<double>::infinity()) != settings.p2.end() and
           settings.atom2 == std::numeric_limits<int>::infinity()) {
    throw SerenityError((std::string) "The plot for one point is not implemented yet!");
    _pointTask = true;
  }
  else if (std::find(settings.p3.begin(), settings.p3.end(), std::numeric_limits<double>::infinity()) != settings.p3.end() and
           settings.atom3 == std::numeric_limits<int>::infinity()) {
    throw SerenityError((std::string) "The plot for a linear grid is not implemented yet!");
    _lineTask = true;
  }
  else if (std::find(settings.p4.begin(), settings.p4.end(), std::numeric_limits<double>::infinity()) != settings.p4.end() and
           settings.atom4 == std::numeric_limits<int>::infinity()) {
    _planeTask = true;
  }
  else {
    throw SerenityError((std::string) "The plot for a cubic grid defined with four points is not implemented yet!");
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
  std::shared_ptr<Geometry> geom = std::make_shared<Geometry>(*(supersystem->getGeometry()));
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
    printSubSectionTitle((std::string) "Printing Data on the Molecular Cavity");
    auto cavityGrid = supersystem->getMolecularSurface(MOLECULAR_SURFACE_TYPES::ACTIVE);
    plotWriter = std::make_shared<MolecularSurfaceGridFileWriter>(sysSettings, cavityGrid);
  }
  else if (_cubeTask) {
    printSubSectionTitle((std::string) "Printing Data to Cube Files");
    if (std::isinf(settings.gridSpacing[0]) or std::isinf(settings.gridSpacing[1]) or std::isinf(settings.gridSpacing[2])) {
      settings.gridSpacing = defaultCubeSpacing;
    }
    plotWriter = std::make_shared<CubeFileWriter>(sysSettings, settings);
    fileExtension = ".cube";
  }
  else if (_planeTask) {
    printSubSectionTitle((std::string) "Printing Data to Dat Files");
    if (std::isinf(settings.gridSpacing[0]) or std::isinf(settings.gridSpacing[1]) or std::isinf(settings.gridSpacing[2])) {
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
      print((std::string) "Printing electron density to file: " + _filename + _fnameSuffix_spin + "Density" + fileExtension);
      plotWriter->writeMatrixToGrid(supersystem->getSystemPath() + _filename + _fnameSuffix_spin + "Density", geom,
                                    basis, density_spin);
    };
    if (SCFMode == UNRESTRICTED) {
      // print total density
      print((std::string) "Printing total electron density to file: " + _filename + "_TotalDensity" + fileExtension);
      plotWriter->writeMatrixToGrid(supersystem->getSystemPath() + _filename + "_TotalDensity", geom, density.total());
      print((std::string) "Printing spin difference density to file: " + _filename + "_SpinDensity" + fileExtension);
      plotWriter->writeMatrixToGrid(supersystem->getSystemPath() + _filename + "_SpinDensity", geom, density.difference());
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
    print((std::string) "Printing MOs to files: " + _filename + "_<a/b>MO<number>" + fileExtension);
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
        fileNames.push_back(supersystem->getSystemPath() + _filename + _fnameSuffix_spin + "MO" + std::to_string(i + 1));
      } // for i
      plotWriter->writeVectorSetToGrid(fileNames, geom, supersystem->getBasisController(), moCoefficients);
    };
  }

  if (settings.orbitals.size() > 0) {
    printSmallCaption("Molecular Orbitals");
    print((std::string) "Printing MOs to files: " + _filename + "_<a/b>MO<number>" + fileExtension);
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
          fileNames.push_back(supersystem->getSystemPath() + _filename + _fnameSuffix_spin + "MO" + std::to_string(orbital));
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
    plotWriter->writeFile(supersystem->getSystemPath() + _filename + "_ESP", geom, calculateESPOnGrid);
  }
  if (settings.elf) {
    ELFCalculator<SCFMode> elf(supersystem);
    auto calculateElfOnGrid = [&](std::shared_ptr<GridController> gridController) {
      GridPotential<RESTRICTED> elfGrid(gridController);
      elfGrid = elf.calculateTotalELFOnGrid(gridController);
      return elfGrid;
    };
    plotWriter->writeFile(supersystem->getSystemPath() + _filename + "_ELF", geom, calculateElfOnGrid);
  }
  if (settings.elfts) {
    ELFCalculator<SCFMode> elf(supersystem);
    auto calculateElfTSOnGrid = [&](std::shared_ptr<GridController> gridController) {
      GridPotential<RESTRICTED> elfGrid(gridController);
      elfGrid = elf.calculateTotalELFTSOnGrid(gridController);
      return elfGrid;
    };
    plotWriter->writeFile(supersystem->getSystemPath() + _filename + "_ELF", geom, calculateElfTSOnGrid);
  }
  if (settings.sedd) {
    SEDD<SCFMode> sedd;
    auto lambda = sedd.getSEDDLambda(supersystem);
    plotWriter->writeFile(supersystem->getSystemPath() + _filename + "_SEDD", geom, lambda);
  }
  if (settings.dori) {
    SEDD<SCFMode> dori;
    auto lambda = dori.getDORILambda(supersystem);
    plotWriter->writeFile(supersystem->getSystemPath() + _filename + "_DORI", geom, lambda);
  }
  if (settings.signedDensity) {
    SEDD<SCFMode> dori;
    auto lambda = dori.getSignedDensityLambda(supersystem);
    plotWriter->writeFile(supersystem->getSystemPath() + _filename + "_signedDensity", geom, lambda);
  }
  if (settings.ntos) {
    NTOCalculator<SCFMode> ntoCalculator(_systems, _environmentSystems, settings.ntoPlotThreshold);
    for (unsigned i : settings.excitations) {
      if (!i)
        // this error also catches the case of the input not specifying any excitations, since the default is {0}
        throw SerenityError("You need to specify the desired excitations in the PlotTask's input (with 1 denoting the "
                            "lowest-lying excited state)");
      if (ntoCalculator.getNumberOfStates() < i)
        throw SerenityError("You need to run an LRSCFTask with neigen greater or equal to the PlotTask's highest "
                            "excitationnumber! The requested excitation " +
                            std::to_string(i) + " does not fall within the calculated " +
                            std::to_string(ntoCalculator.getNumberOfStates()));
      i--;
      // From here on, excitations are (internally) indexed starting from 0, directly referring to LRSCF
      // excitationvector columns
      for (unsigned int iSys = 0; iSys < _systems.size(); iSys++) {
        geom = _systems[iSys]->getGeometry();
        std::shared_ptr<BasisController> basisController = _systems[iSys]->getBasisController();
        // NTO data
        const std::string dirName = ntoCalculator.getDir(i, iSys);
        auto oNTOs = ntoCalculator.getOccNTOs(i, iSys);
        auto vNTOs = ntoCalculator.getVirtNTOs(i, iSys);
        auto oEigenvalues = ntoCalculator.getOccEigenvalues(i);
        auto vEigenvalues = ntoCalculator.getVirtEigenvalues(i);
        for_spin(oNTOs, vNTOs, oEigenvalues, vEigenvalues, _fnameSuffix) {
          std::vector<std::string> fileNames;
          std::vector<unsigned int> indicesToPlot;
          for (unsigned int j = 0; j < oNTOs_spin.cols(); j++) {
            if (oEigenvalues_spin(j) > settings.ntoPlotThreshold) {
              int tmpState = j + 1;
              std::string fileName = dirName + tmpState + _fnameSuffix_spin + "occ";
              fileNames.push_back(fileName);
              indicesToPlot.push_back(j);
            }
          }

          Eigen::MatrixXd ntoToPlot = Eigen::MatrixXd::Zero(oNTOs_spin.rows(), indicesToPlot.size());
          unsigned int counter = 0;
          for (auto index : indicesToPlot) {
            ntoToPlot.col(counter) = oNTOs_spin.col(index);
            counter++;
          }
          plotWriter->writeVectorSetToGrid(fileNames, geom, basisController, ntoToPlot);
          fileNames.resize(0);
          indicesToPlot.resize(0);
          for (unsigned int j = 0; j < vNTOs_spin.cols(); j++) {
            if (vEigenvalues_spin(j) > settings.ntoPlotThreshold) {
              std::string tmpState = std::to_string(vEigenvalues_spin.size() - j + oEigenvalues_spin.size());
              std::string fileName = dirName + tmpState + _fnameSuffix_spin + "virt";
              fileNames.push_back(fileName);
              indicesToPlot.push_back(j);
            }
          }
          ntoToPlot = Eigen::MatrixXd::Zero(vNTOs_spin.rows(), indicesToPlot.size());
          counter = 0;
          for (auto index : indicesToPlot) {
            ntoToPlot.col(counter) = vNTOs_spin.col(index);
            counter++;
          }
          plotWriter->writeVectorSetToGrid(fileNames, geom, basisController, ntoToPlot);
        };

        // Transition density
        const MatrixInBasis<SCFMode>& transitionDensity = ntoCalculator.getTransitionDensity(i, iSys);
        for_spin(transitionDensity, _fnameSuffix) {
          std::string filename = _systems[iSys]->getSystemPath() + _systems[iSys]->getSystemName() +
                                 ((_systems.size() > 1) ? "_FDEc" : "") + "_TransitionDensity" + _fnameSuffix_spin +
                                 std::to_string(i + 1);
          OutputControl::nOut << "\n  Printing transition density to file: " + filename + ".cube" << std::endl;
          plotWriter->writeMatrixToGrid(filename, geom, basisController, transitionDensity_spin);
        };

        // Hole- and particle-densities
        const MatrixInBasis<SCFMode>& holeDensity = ntoCalculator.getHoleDensity(i, iSys);
        const MatrixInBasis<SCFMode>& particleDensity = ntoCalculator.getParticleDensity(i, iSys);
        const MatrixInBasis<SCFMode>& holeDensityCorr = ntoCalculator.getHoleDensityCorrection(i, iSys);
        for_spin(holeDensity, particleDensity, holeDensityCorr, _fnameSuffix) {
          std::string holename = _systems[iSys]->getSystemPath() + _systems[iSys]->getSystemName() + "_HoleDensity" +
                                 ((_systems.size() > 1) ? "_FDEc" : "") + _fnameSuffix_spin + std::to_string(i + 1);
          std::string particlename = _systems[iSys]->getSystemPath() + _systems[iSys]->getSystemName() + "_ParticleDensity" +
                                     ((_systems.size() > 1) ? "_FDEc" : "") + _fnameSuffix_spin + std::to_string(i + 1);
          OutputControl::nOut << "  Printing hole density to file: " + holename + ".cube" << std::endl;
          plotWriter->writeMatrixToGrid(holename, geom, basisController, holeDensity_spin);
          OutputControl::nOut << "  Printing particle density to file: " + particlename + ".cube" << std::endl;
          plotWriter->writeMatrixToGrid(particlename, geom, basisController, particleDensity_spin);
          if (_systems.size() > 1) {
            OutputControl::nOut << "  Printing hole density correction to file: " + holename + "corr.cube" << std::endl;
            plotWriter->writeMatrixToGrid(holename + "corr", geom, basisController, holeDensityCorr_spin);
          }
        };
      }
    }
  }
  if (settings.cctrdens || settings.ccexdens) {
    printBigCaption("CC2 Densities");
    for (auto sys : _systems) {
      auto basis = sys->getBasisController();
      std::string fileName = sys->getSystemPath() + sys->getSystemName() + "_cc2_dens.";
      if (_systems.size() > 1)
        fileName += "fdec.";
      fileName += (SCFMode == RESTRICTED) ? "res." : "unres.";
      fileName += "h5";

      const auto& coeffMat = sys->getActiveOrbitalController<SCFMode>()->getCoefficients();

      HDF5::H5File file(fileName.c_str(), H5F_ACC_RDONLY);

      Eigen::MatrixXd stateDens, leftTransDens, rightTransDens;

      HDF5::dataset_exists(file, "State Densities");
      HDF5::load(file, "State Densities", stateDens);
      // coeffMat.size() is used deliberately here, since the number of rows should be nMO * nMO
      if (stateDens.rows() != coeffMat.size())
        throw SerenityError("The dimension of your loaded state densities (" + std::to_string(stateDens.rows()) +
                            ") does not match with the size of the coefficient matrix (" +
                            std::to_string(coeffMat.size()) + ").");
      fileName = fileName.substr(0, fileName.size() - 3);
      unsigned moStart = 0;
      for_spin(coeffMat, _fnameSuffix) {
        Eigen::MatrixXd stateDensToPlot =
            (Eigen::MatrixXd)coeffMat_spin *
            Eigen::Map<Eigen::MatrixXd>(stateDens.col(0).data() + moStart, coeffMat_spin.rows(), coeffMat_spin.cols()) *
            coeffMat_spin.transpose();
        OutputControl::nOut << "  Printing ground state CC2 density to file: " + fileName + _fnameSuffix_spin + "gs.cube\n"
                            << std::endl;
        plotWriter->writeMatrixToGrid(fileName + _fnameSuffix_spin + "gs", geom, basis, stateDensToPlot);
        moStart += coeffMat_spin.size();
      };
      if (settings.cctrdens) {
        HDF5::load(file, "Right Transition Densities", rightTransDens);
        if (rightTransDens.rows() != coeffMat.size())
          throw SerenityError(
              "The dimension of your loaded right transition densities (" + std::to_string(rightTransDens.rows()) +
              ") does not match with the size of the coefficient matrix (" + std::to_string(coeffMat.size()) + ").");
        HDF5::load(file, "Left Transition Densities", leftTransDens);
        if (leftTransDens.rows() != coeffMat.size())
          throw SerenityError(
              "The dimension of your loaded left transition densities (" + std::to_string(leftTransDens.rows()) +
              ") does not match with the size of the coefficient matrix (" + std::to_string(coeffMat.size()) + ").");
        file.close();
        for (unsigned iExc : settings.excitations) {
          if (!iExc)
            throw SerenityError(
                "You need to specify the desired excitations in the PlotTask's input (with 1 denoting the "
                "lowest-lying excited state)");
          iExc--;
          moStart = 0;
          for_spin(coeffMat, _fnameSuffix) {
            Eigen::MatrixXd rightTransDensToPlot = (Eigen::MatrixXd)coeffMat_spin *
                                                   Eigen::Map<Eigen::MatrixXd>(rightTransDens.col(iExc).data() + moStart,
                                                                               coeffMat_spin.rows(), coeffMat_spin.cols()) *
                                                   coeffMat_spin.transpose();
            OutputControl::nOut << "  Printing right CC2 transition density to file: " + fileName + _fnameSuffix_spin +
                                       "rightTrans_" + std::to_string(iExc + 1) + ".cube\n"
                                << std::endl;
            plotWriter->writeMatrixToGrid(fileName + _fnameSuffix_spin + "rightTrans_" + std::to_string(iExc + 1), geom,
                                          basis, rightTransDensToPlot);
            Eigen::MatrixXd leftTransDensToPlot = (Eigen::MatrixXd)coeffMat_spin *
                                                  Eigen::Map<Eigen::MatrixXd>(leftTransDens.col(iExc).data() + moStart,
                                                                              coeffMat_spin.rows(), coeffMat_spin.cols()) *
                                                  coeffMat_spin.transpose();
            OutputControl::nOut << "  Printing left CC2 transition density to file: " + fileName + _fnameSuffix_spin +
                                       "leftTrans_" + std::to_string(iExc + 1) + ".cube\n"
                                << std::endl;
            plotWriter->writeMatrixToGrid(fileName + _fnameSuffix_spin + "leftTrans_" + std::to_string(iExc + 1), geom,
                                          basis, leftTransDensToPlot);
            moStart += coeffMat_spin.size();
          };
        } /* for excitations */
      }
      if (settings.ccexdens) {
        for (unsigned iExc : settings.excitations) {
          if (!iExc)
            throw SerenityError(
                "You need to specify the desired excitations in the PlotTask's input (with 1 denoting the "
                "lowest-lying excited state)");
          iExc--;
          moStart = 0;
          for_spin(coeffMat, _fnameSuffix) {
            Eigen::MatrixXd stateDensToPlot = (Eigen::MatrixXd)coeffMat_spin *
                                              Eigen::Map<Eigen::MatrixXd>(stateDens.col(iExc + 1).data() + moStart,
                                                                          coeffMat_spin.rows(), coeffMat_spin.cols()) *
                                              coeffMat_spin.transpose();
            OutputControl::nOut << "  Printing excited state CC2 density to file: " + fileName + _fnameSuffix_spin +
                                       std::to_string(iExc + 1) + ".cube\n"
                                << std::endl;
            plotWriter->writeMatrixToGrid(fileName + _fnameSuffix_spin + std::to_string(iExc + 1), geom, basis, stateDensToPlot);
            moStart += coeffMat_spin.size();
          };
        }
      }
    } /* for systems */
  }
  if (settings.gridCoordinates) {
    auto lambda = [&](std::shared_ptr<GridController> gridController) { return gridController->getWeights(); };
    plotWriter->writeFile(supersystem->getSystemPath() + _filename + "_coords", geom, lambda);
  }

  if (settings.nros) {
    for (auto system : _systems) {
      LRSCFTaskSettings lrscfSettings;
      lrscfSettings.loadType = Options::LRSCF_TYPE::ISOLATED;
      auto lrscfcontroller = std::make_shared<LRSCFController<SCFMode>>(system, lrscfSettings);
      std::vector<Eigen::MatrixXd> XY(2);
      Eigen::VectorXd freque;
      std::string fileName = system->getSystemPath() + system->getSystemName() + "_lrscf_resp.";
      if (lrscfSettings.loadType == Options::LRSCF_TYPE::ISOLATED) {
        fileName += "iso.";
      }
      else if (lrscfSettings.loadType == Options::LRSCF_TYPE::UNCOUPLED) {
        fileName += "fdeu.";
      }
      fileName += (SCFMode == RESTRICTED) ? "res." : "unres.";
      fileName += "h5";
      HDF5::H5File afile(fileName, H5F_ACC_RDONLY);
      HDF5::dataset_exists(afile, "X+Y");
      HDF5::dataset_exists(afile, "X-Y");
      HDF5::dataset_exists(afile, "frequencies");
      HDF5::load(afile, "X+Y", XY[0]);
      HDF5::load(afile, "X-Y", XY[1]);
      HDF5::load(afile, "frequencies", freque);
      afile.close();
      auto nocc = lrscfcontroller->getNOccupied();
      auto nvirt = lrscfcontroller->getNVirtual();
      auto geom = system->getGeometry();
      auto basis = lrscfcontroller->getBasisController();
      DipoleIntegrals<SCFMode> dip({lrscfcontroller}, Point(0.0, 0.0, 0.0));
      Eigen::MatrixXd dipolelengths = *(dip.getLengths());
      NROCalculator<SCFMode> nro(XY, lrscfcontroller);
      for (unsigned iFreq = 0; iFreq < XY[0].cols() / 3; iFreq++) {
        auto NROs = nro.getNROs(iFreq);
        Eigen::MatrixXd xpy = XY[0].middleCols(iFreq * 3, 3);
        Eigen::MatrixXd xmy = XY[1].middleCols(iFreq * 3, 3);
        SpinPolarizedData<SCFMode, Eigen::MatrixXd> singularvalues;
        unsigned ia_start = 0;
        for_spin(nocc, nvirt, singularvalues) {
          singularvalues_spin = Eigen::MatrixXd::Zero(nocc_spin, 3);
          for (unsigned j = 0; j < 3; j++) {
            Eigen::MatrixXd xpymatrix =
                (Eigen::Map<Eigen::MatrixXd>(xpy.col(j).data() + ia_start, nvirt_spin, nocc_spin)).transpose();
            Eigen::JacobiSVD<Eigen::MatrixXd> svdp(xpymatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
            singularvalues_spin.col(j) = svdp.singularValues() / (svdp.singularValues().sum());
          }
        };

        // write NRO-pairs to disk
        PlotTaskSettings plottasksettings;
        plottasksettings.gridSpacing = {0.12, 0.12, 0.12};
        CubeFileWriter writer = CubeFileWriter(system->getSettings(), plottasksettings);
        std::vector<std::string> partfilenames;
        std::vector<std::string> holefilenames;
        std::vector<std::string> SpatDirection = {"x", "y", "z"};

        std::string spin = (SCFMode == RESTRICTED) ? "" : " alpha";
        for_spin(singularvalues, NROs) {
          for (unsigned row = 0; row < 3; row++) {
            double accSingularValues = 0;
            unsigned plotcounter = 0;
            for (unsigned i = 0; i < singularvalues_spin.rows(); i++) {
              partfilenames = {};
              holefilenames = {};
              if (accSingularValues < settings.nrominimum) {
                partfilenames.push_back(system->getSystemPath() + "freq_" + std::to_string(iFreq + 1) +
                                        SpatDirection[row] + "particleNRO_" + std::to_string(i + 1) + spin);
                holefilenames.push_back(system->getSystemPath() + "freq_" + std::to_string(iFreq + 1) +
                                        SpatDirection[row] + "holeNRO_" + std::to_string(i + 1) + spin);
                OutputControl::nOut << "Printing to files " << partfilenames[0] << " and " << holefilenames[0] << std::endl;
                Eigen::MatrixXd parts = NROs_spin[row * 2 + 1].col(i);
                Eigen::MatrixXd holes = NROs_spin[row * 2].col(i);

                writer.writeVectorSetToGrid(partfilenames, geom, basis, parts);
                writer.writeVectorSetToGrid(holefilenames, geom, basis, holes);
                plotcounter++;
                accSingularValues += singularvalues_spin(i, row);
              }
            }
            if (generalSettings.printLevel >= Options::GLOBAL_PRINT_LEVELS::NORMAL) {
              printf("For the frequency %f in direction %s, %i NRO pairs were printed to cube files accounting for "
                     "approximately %.2f%% of the%s response.\n",
                     freque[iFreq], SpatDirection[row].c_str(), plotcounter, accSingularValues * 100, spin.c_str());
            }
          }
          spin = (SCFMode == RESTRICTED) ? "" : " beta";
        };
      } /* frequency loop */
    }   /* system loop */
  }     /* settings.nrotest */

} /* run */

template class PlotTask<Options::SCF_MODES::RESTRICTED>;
template class PlotTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
