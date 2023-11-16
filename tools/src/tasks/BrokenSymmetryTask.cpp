/**
 * @file BrokenSymmetryTask.cpp
 *
 * @date Feb 24, 2020
 * @author Anja Massolle
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
#include "tasks/BrokenSymmetryTask.h"
/* Include Serenity Internal Headers */
#include "analysis/brokenSymmetry/BrokenSymmetry.h"
#include "geometry/Geometry.h"
#include "io/FormattedOutputStream.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/SystemAdditionTask.h"
/* Include Std and External Headers */

namespace Serenity {

BrokenSymmetryTask::BrokenSymmetryTask(std::vector<std::shared_ptr<SystemController>> hsSystemController,
                                       std::vector<std::shared_ptr<SystemController>> bsSystemController)
  : _hsSystemController(hsSystemController), _bsSystemController(bsSystemController) {
}

void BrokenSymmetryTask::run() {
  printSectionTitle("Broken-Symmetry calculation");
  std::shared_ptr<BrokenSymmetry> bsCalculation = nullptr;
  if (_bsSystemController.size() != 0) {
    if (_bsSystemController.size() != _hsSystemController.size()) {
      throw SerenityError("Loading of a different number of HS and BS states is currently not implemented, please use "
                          "the same amount of systems for the HS and BS states!");
    }
    if (_bsSystemController.size() == 1) {
      if (_hsSystemController[0]->hasElectronicStructure<UNRESTRICTED>()) {
        OutputControl::nOut << "load " << (_hsSystemController[0]->getSettings()).name << " as HS system " << std::endl;
      }
      else {
        throw SerenityError("Please load only systems which possess an unrestricted electronic structure!");
      }
      if (_bsSystemController[0]->hasElectronicStructure<UNRESTRICTED>()) {
        OutputControl::nOut << "load " << (_bsSystemController[0]->getSettings()).name << " as BS system " << std::endl;
      }
      else {
        throw SerenityError("Please load only systems which possess an unrestricted electronic structure!");
      }
    }
    else if (_bsSystemController.size() == 2) {
      if (_hsSystemController[0]->hasElectronicStructure<UNRESTRICTED>() and
          _hsSystemController[1]->hasElectronicStructure<UNRESTRICTED>()) {
        OutputControl::nOut << "load " << (_hsSystemController[0]->getSettings()).name << " and "
                            << (_hsSystemController[1]->getSettings()).name << " as HS system " << std::endl;
      }
      else {
        throw SerenityError("Please load only systems which possess an unrestricted electronic structure!");
      }
      if (_bsSystemController[0]->hasElectronicStructure<UNRESTRICTED>() and
          _bsSystemController[1]->hasElectronicStructure<UNRESTRICTED>()) {
        OutputControl::nOut << "load " << (_bsSystemController[0]->getSettings()).name << " and "
                            << (_bsSystemController[1]->getSettings()).name << " as BS system " << std::endl;
      }
      else {
        throw SerenityError("Please load only systems which possess an unrestricted electronic structure!");
      }
    }
    else {
      throw SerenityError("Loading of more than two BS spin sites not implemented!");
    }
    bsCalculation =
        std::make_shared<BrokenSymmetry>(_hsSystemController, settings.embedding, settings.evalTsOrtho,
                                         settings.evalAllOrtho, settings.orthogonalizationScheme, _bsSystemController);
  }
  else {
    if (_hsSystemController.size() == 1) {
      bsCalculation = std::make_shared<BrokenSymmetry>(_hsSystemController);
      bsCalculation->bsDFT(settings.nA, settings.nB);
    } /*_hsSystemController.size() == 1*/
    else if (_hsSystemController.size() == 2) {
      switch (settings.embeddingScheme) {
        case Options::EMBEDDING_SCHEME::NONE: {
          /*
           * Setup supersystem
           */
          auto superSystemSettings = _hsSystemController[0]->getSettings();
          superSystemSettings.name = "HighSpinSystem";
          superSystemSettings.charge = 0;
          superSystemSettings.spin = 0;
          auto superSystem = std::make_shared<SystemController>(std::make_shared<Geometry>(), superSystemSettings);

          SystemAdditionTask<UNRESTRICTED> hsSysAdd(superSystem, _hsSystemController);
          hsSysAdd.settings.addOccupiedOrbitals = false;
          hsSysAdd.run();

          settings.nA = _hsSystemController[0]->getSpin();
          settings.nB = _hsSystemController[1]->getSpin();

          _hsSystemController = {superSystem};
          bsCalculation = std::make_shared<BrokenSymmetry>(_hsSystemController);
          bsCalculation->bsDFT(settings.nA, settings.nB);
          break;
        } /*standard BS-DFT*/

        case Options::EMBEDDING_SCHEME::ISOLATED: {
          bsCalculation = std::make_shared<BrokenSymmetry>(_hsSystemController, settings.embedding, settings.evalTsOrtho,
                                                           settings.evalAllOrtho, settings.orthogonalizationScheme);
          bsCalculation->bsIsolated();
          break;
        }

        case Options::EMBEDDING_SCHEME::FDE: {
          bsCalculation = std::make_shared<BrokenSymmetry>(_hsSystemController, settings.embedding, settings.evalTsOrtho,
                                                           settings.evalAllOrtho, settings.orthogonalizationScheme);
          bsCalculation->bsFDE();
          break;
        }

        case Options::EMBEDDING_SCHEME::FAT: {
          FreezeAndThawTaskSettings fatSettings;
          fatSettings.maxCycles = settings.maxCycles;
          fatSettings.convThresh = settings.convThresh;
          bsCalculation = std::make_shared<BrokenSymmetry>(_hsSystemController, settings.embedding, settings.evalTsOrtho,
                                                           settings.evalAllOrtho, settings.orthogonalizationScheme);
          bsCalculation->bsFAT(fatSettings);
          break;
        }

        default: {
          SerenityError("No implemented embedding scheme was selected.");
        }
      } /* switch(settings.embeddingScheme) */
    }   /*_hsSystemController.size() == 2*/
  }
  auto s2HS = bsCalculation->getS2HS();
  printSmallCaption("<S*S>");
  printf("\n      <S*S> Highspin System              = %1.3f \n", s2HS);
  printf("      <S*S> Broken-Symmetry System       = %1.3f \n\n", bsCalculation->getS2BS());
  printSmallCaption("J-coupling");
  OutputControl::nOut << "J(1) = (E_BS - E_HS)/(S_max * S_max)" << std::endl;
  OutputControl::nOut << "J(2) = (E_BS - E_HS)/(S_max * (S_max + 1))" << std::endl;
  OutputControl::nOut << "J(3) = (E_BS - E_HS)/(<S*S>_HS - <S*S>_BS)" << std::endl;
  if (settings.evalAllOrtho or _hsSystemController.size() == 1) {
    OutputControl::nOut << "J(4) = (E_BS - E_HS)/(1 + (S*S)_AB)" << std::endl;
  }
  printf("\n      J(1)                               = %18.2f cm^(-1)\n", bsCalculation->getJ1());
  printf("      J(2)                               = %18.2f cm^(-1)\n", bsCalculation->getJ2());
  printf("      J(3) (<S*S> from electron density) = %18.2f cm^(-1)\n\n", bsCalculation->getJ3());
  if (settings.evalAllOrtho or _hsSystemController.size() == 1) {
    printf("\n      <S*S> UHF Highspin System              = %1.3f \n", bsCalculation->getS2uhfHS());
    printf("      <S*S> UHF Broken-Symmetry System       = %1.3f \n\n", bsCalculation->getS2uhfBS());
    printf("      J(3)   (<S*S> from UKS orbitals)   = %18.2f cm^(-1)\n", bsCalculation->getJ3UHF());
    printf("      Sab                                = %2.3f\n", bsCalculation->getSab());
    printf("      J(4)                               = %18.2f cm^(-1)\n", bsCalculation->getJ4());
  }
} /*run*/

} /*namespace Serenity*/
