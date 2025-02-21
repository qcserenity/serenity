/**
 * @file DipoleApproximationToPairEnergies_test.cpp
 *
 * @date Feb 4, 2019
 * @author Moritz Bensberg
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

/* Include Serenity Internal Headers */
#include "postHF/MPn/DipoleApproximationToPairEnergies.h"
#include "analysis/PAOSelection/BoughtonPulayAlgorithm.h"             //Construct PAO domains.
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h" //Mulliken populations for BP-algorithm.
#include "basis/AtomCenteredBasisController.h"                        //Mulliken populations for BP-algorithm.
#include "data/ElectronicStructure.h"                                 //Orbitals.
#include "data/OrbitalController.h"                                   //Orbital coefficients.
#include "data/OrbitalPair.h"                                         //Results will be stored in the pairs.
#include "data/PAOController.h"                       //Needed for construction of DipoleApproximationToPairEnergies
#include "integrals/OneElectronIntegralController.h"  //Get overlap integrals.
#include "postHF/MPn/LocalMP2.h"                      //Compare to LMP2 pair energies.
#include "potentials/bundles/HFPotentials.h"          //Fock matrix.
#include "system/SystemController.h"                  //Test systems.
#include "tasks/LocalizationTask.h"                   //Orbital localization.
#include "tasks/ScfTask.h"                            //Clean electronic structures.
#include "testsupply/SystemController__TEST_SUPPLY.h" //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //Testing framework.

namespace Serenity {

class DipoleApproximationToPairEnergiesTest : public ::testing::Test {
 protected:
  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Minimal test example. Only one orbital pair.
 */
TEST_F(DipoleApproximationToPairEnergiesTest, WaterDimer) {
  const auto scfMode = Options::SCF_MODES::RESTRICTED;
  auto act1 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonOne_Def2_SVP);
  auto act2 = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::WaterMonTwo_Def2_SVP);
  auto waterDimer = *act1 + *act2;
  double e = waterDimer->getElectronicStructure<scfMode>()->getEnergy();
  (void)e;
  ScfTask<scfMode> scfTask(waterDimer);
  scfTask.run();
  LocalizationTask locTask(waterDimer);
  locTask.settings.locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::IBO;
  locTask.run();
  auto s_AO =
      std::make_shared<MatrixInBasis<scfMode>>(waterDimer->getOneElectronIntegralController()->getOverlapIntegrals());
  unsigned int nOcc = waterDimer->getNOccupiedOrbitals<scfMode>();
  auto occCoef = std::make_shared<Eigen::MatrixXd>(
      waterDimer->getActiveOrbitalController<scfMode>()->getCoefficients().leftCols(nOcc).eval());
  auto paoController = std::make_shared<PAOController>(
      std::make_shared<DensityMatrix<scfMode>>(waterDimer->getElectronicStructure<scfMode>()->getDensityMatrix()), s_AO, 1e-6);
  auto sysPot = waterDimer->getPotentials<scfMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>();
  auto fockMatrix = std::make_shared<FockMatrix<scfMode>>(sysPot->getFockMatrix(
      waterDimer->getElectronicStructure<scfMode>()->getDensityMatrix(), std::make_shared<EnergyComponentController>()));

  auto mullikenPops =
      std::make_shared<Eigen::MatrixXd>(MullikenPopulationCalculator<scfMode>::calculateAtomwiseOrbitalPopulations(
                                            waterDimer->getActiveOrbitalController<scfMode>()->getCoefficients(), *s_AO,
                                            waterDimer->getAtomCenteredBasisController()->getBasisIndices())
                                            .leftCols(nOcc)
                                            .eval());
  BoughtonPulayAlgorithm bpAlgorithm(waterDimer->getOneElectronIntegralController(),
                                     waterDimer->getAtomCenteredBasisController(), mullikenPops, occCoef, 0.002);

  auto paoToOcc = bpAlgorithm.selectPAOs();
  DipoleApproximationToPairEnergies approxCalculator(waterDimer->getBasisController(), s_AO, occCoef, paoController,
                                                     fockMatrix, paoToOcc, 1e-6);
  // Full local MP2 treatment
  LocalCorrelationSettings lCSettings;
  lCSettings.completenessThreshold = 0.002;
  lCSettings.orbitalToShellThreshold = 1e-5;
  lCSettings.mullikenThreshold = 1e-3;
  lCSettings.ccsdPairThreshold = 0.0;
  auto lCController = std::make_shared<LocalCorrelationController>(waterDimer, lCSettings);
  LocalMP2 localMP2(lCController);
  double localMP2Energy = localMP2.calculateEnergyCorrection().sum();
  auto localMP2OrbitalPairs = lCController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  (void)localMP2Energy;
  // Select orbital pairs between water molecules.
  std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs;
  for (const auto& lMP2Pair : localMP2OrbitalPairs) {
    unsigned int i = lMP2Pair->i;
    unsigned int j = lMP2Pair->j;
    int test = paoToOcc->col(i).cwiseProduct(paoToOcc->col(j)).sum();
    if (test == 0) {
      orbitalPairs.push_back(lMP2Pair);
    }
  }
  approxCalculator.calculateDipoleApproximation(orbitalPairs);

  double dipolePairEnergies = 0.0;
  for (const auto& pair : orbitalPairs) {
    printf("Pair: %d%s%d %s %18f %s\n", pair->i, "/", pair->j, "Pair energy app.", pair->dipolePairEnergy * 2625.50, " kJ/mol");
    printf("Pair: %d%s%d %s %18f %s\n", pair->i, "/", pair->j, "Pair energy LMP2", pair->lMP2PairEnergy * 2625.50, " kJ/mol");
    // Approximation and actual MP2 pair energies should be in the same order of magnitude.
    EXPECT_NEAR(pair->dipolePairEnergy, pair->lMP2PairEnergy, 5e-5);
    dipolePairEnergies += pair->dipolePairEnergy;
  }
  EXPECT_NEAR(dipolePairEnergies, -7.474047162993051e-05, 3e-7);
  std::cout << "Total interaction energy app. [kJ/mol]: " << dipolePairEnergies * 2625.50 << std::endl;
  // clean up
  auto supersystemName = waterDimer->getSystemName();
  SystemController__TEST_SUPPLY::cleanUpSystemDirectory(waterDimer);
  SystemController__TEST_SUPPLY::cleanUp();
}

} /* namespace Serenity */
