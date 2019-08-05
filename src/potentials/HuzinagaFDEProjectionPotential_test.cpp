/**
 * @file HuzinagaFDEProjectionPotential_test.cpp
 *
 * @date Feb 23, 2018
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
/* Include Serenity Internal Headers */
#include "potentials/HuzinagaFDEProjectionPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
#include "grid/GridControllerFactory.h"
#include "settings/Settings.h"
#include "geometry/Geometry.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>
namespace Serenity {

/**
 * @class HuzinagaFDEProjectionPotentialTest
 * @brief Sets everything up for the tests of HuzinagaFDEProjectionPotential.h/.cpp .
 */
class HuzinagaFDEProjectionPotentialTest : public ::testing::Test {
 protected:
	HuzinagaFDEProjectionPotentialTest() {
  }

  virtual ~HuzinagaFDEProjectionPotentialTest() = default;

  static void TearDownTestCase() {
   SystemController__TEST_SUPPLY::cleanUp();
  }
};

/**
 * @test
 * @brief Tests HuzinagaFDEProjectionPotential.h/.cpp: Tests restricted fock matrix (Huzinaga).
 */
TEST_F(HuzinagaFDEProjectionPotentialTest, restricted_huzinaga) {
	auto act = SystemController__TEST_SUPPLY::getSystemController(
			TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
	auto env = SystemController__TEST_SUPPLY::getSystemController(
			TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

	HuzinagaFDEProjectionPotential<Options::SCF_MODES::RESTRICTED> huzFDEPotential(
			act,{env},FunctionalClassResolver::resolveFunctional(act->getSettings().dft.functional));
	auto f = huzFDEPotential.getMatrix();
	EXPECT_NEAR(0.052226824418551698,f(0,0),1e-7);
	EXPECT_NEAR(0.11198231040394366,f(0,1),1e-7);
	EXPECT_NEAR(0.11198231040394366,f(1,0),1e-7);
	EXPECT_NEAR(0.22587526848483189,f(1,1),1e-7);

	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
	SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests HuzinagaFDEProjectionPotential.h/.cpp: Tests unrestricted fock matrix (Huzinaga).
 */
TEST_F(HuzinagaFDEProjectionPotentialTest, unrestricted_huzinaga) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  HuzinagaFDEProjectionPotential<Options::SCF_MODES::UNRESTRICTED> huzFDEPotential(
      act,{env},
      FunctionalClassResolver::resolveFunctional(act->getSettings().dft.functional),
      false,0.0,
      false);
  auto f = huzFDEPotential.getMatrix();
  EXPECT_NEAR(0.05166910663948323,f.alpha(0,0),1e-7);
  EXPECT_NEAR(0.10986760211800925,f.alpha(0,1),1e-7);
  EXPECT_NEAR(0.10986760211800925,f.alpha(1,0),1e-7);
  EXPECT_NEAR(0.21856343934609146,f.alpha(1,1),1e-7);

  EXPECT_NEAR(0.05166910663948323,f.beta(0,0),1e-7);
  EXPECT_NEAR(0.10986760211800925,f.beta(0,1),1e-7);
  EXPECT_NEAR(0.10986760211800925,f.beta(1,0),1e-7);
  EXPECT_NEAR(0.21856343934609146,f.beta(1,1),1e-7);
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests HuzinagaFDEProjectionPotential.h/.cpp: Tests restricted fock matrix (Hoffmann) with unrestricted environment.
 */
TEST_F(HuzinagaFDEProjectionPotentialTest, restricted_unrestricted_hoffmann) {
	auto act = SystemController__TEST_SUPPLY::getSystemController(
			TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
	auto env = SystemController__TEST_SUPPLY::getSystemController(
			TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
	env->setSCFMode(Options::SCF_MODES::UNRESTRICTED);
	HuzinagaFDEProjectionPotential<Options::SCF_MODES::RESTRICTED> huzFDEPotential(
			act,{env},
			FunctionalClassResolver::resolveFunctional(act->getSettings().dft.functional),
			false,0.0,
			true);
	auto f = huzFDEPotential.getMatrix();
	EXPECT_NEAR(0.041308143680642548,f(0,0),1e-7);
	EXPECT_NEAR(0.081439053107599857,f(0,1),1e-7);
	EXPECT_NEAR(0.081439053107599857,f(1,0),1e-7);
	EXPECT_NEAR(0.14191046489412251,f(1,1),1e-7);
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.settings").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.xyz").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.energies.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.orbs.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE.dmat.res.h5").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE/out").c_str());
	std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+OH_MINBAS_PBE").c_str());
	SystemController__TEST_SUPPLY::cleanUp();
}


/**
 * @test
 * @brief Tests HuzinagaFDEProjectionPotential.h/.cpp: Tests unrestricted fock matrix (Hoffmann).
 */
TEST_F(HuzinagaFDEProjectionPotentialTest, unrestricted_hoffmann) {
  auto act = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

  env->getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>();
  HuzinagaFDEProjectionPotential<Options::SCF_MODES::UNRESTRICTED> huzFDEPotential(
      act,{env},
      FunctionalClassResolver::resolveFunctional(act->getSettings().dft.functional),
      false,0.0,
      true);
  auto f = huzFDEPotential.getMatrix();
  EXPECT_NEAR(0.039352521735825291,f.alpha(0,0),1e-7);
  EXPECT_NEAR(0.077029478641347057,f.alpha(0,1),1e-7);
  EXPECT_NEAR(0.077029478641347057,f.alpha(1,0),1e-7);
  EXPECT_NEAR(0.13101138153920278,f.alpha(1,1),1e-7);

  EXPECT_NEAR(0.039352521735825291,f.beta(0,0),1e-7);
  EXPECT_NEAR(0.077029478641347057,f.beta(0,1),1e-7);
  EXPECT_NEAR(0.077029478641347057,f.beta(1,0),1e-7);
  EXPECT_NEAR(0.13101138153920278,f.beta(1,1),1e-7);
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.settings").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.xyz").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.energies.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.orbs.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE.dmat.unres.h5").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE/out").c_str());
  std::remove((env->getSettings().path+"TestSystem_H2_6_31Gs_ACTIVE_FDE+TestSystem_H2_6_31Gs_ENVIRONMENT_FDE").c_str());
  SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests HuzinagaFDEProjectionPotential.h/.cpp: Tests restricted fock matrix; with completely truncated projector (Huzinaga).
 */
TEST_F(HuzinagaFDEProjectionPotentialTest, restricted_huzinaga_truncProjection) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
	auto act = SystemController__TEST_SUPPLY::getSystemController(
			TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
	auto env = SystemController__TEST_SUPPLY::getSystemController(
			TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);

	HuzinagaFDEProjectionPotential<SPIN> huzFDEPotential(
			act,{env},
			FunctionalClassResolver::resolveFunctional(act->getSettings().dft.functional),
			true,100,
			false);

	auto f = huzFDEPotential.getMatrix();
	EXPECT_NEAR(0.0,f(0,0),1e-7);
	EXPECT_NEAR(0.0,f(0,1),1e-7);
	EXPECT_NEAR(0.0,f(1,0),1e-7);
	EXPECT_NEAR(0.0,f(1,1),1e-7);

	SystemController__TEST_SUPPLY::cleanUp();
}

/**
 * @test
 * @brief Tests HuzinagaFDEProjectionPotential.h/.cpp: Tests restricted fock matrix; with truncated projector and kin. energy functional (Huzinaga).
 *
 * Identical to the use of a kinetic energy functional for the Pauli repulsion due to the complete truncation of the projection operator.
 */
TEST_F(HuzinagaFDEProjectionPotentialTest, restricted_huzinaga_truncProjection_distantKinFunc) {
  const auto SPIN = Options::SCF_MODES::RESTRICTED;
  auto act = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ACTIVE_FDE);
  auto env = SystemController__TEST_SUPPLY::getSystemController(
      TEST_SYSTEM_CONTROLLERS::H2_6_31Gs_ENVIRONMENT_FDE);
  std::vector<std::shared_ptr<Atom> > superSystemAtoms;
  superSystemAtoms.push_back(act->getGeometry()->getAtoms()[0]);
  superSystemAtoms.push_back(act->getGeometry()->getAtoms()[1]);
  superSystemAtoms.push_back(env->getGeometry()->getAtoms()[0]);
  superSystemAtoms.push_back(env->getGeometry()->getAtoms()[1]);
  auto superGeom = std::make_shared<Geometry>(superSystemAtoms);
  auto supersystemGrid = GridControllerFactory::produce(superGeom,act->getSettings(), Options::GRID_PURPOSES::DEFAULT);

  HuzinagaFDEProjectionPotential<Options::SCF_MODES::RESTRICTED> huzFDEPotential(
      act,{env},
      FunctionalClassResolver::resolveFunctional(act->getSettings().dft.functional),
      true,100,
      false,
      false,
      true,
      2.0,
      2.0,
      supersystemGrid,
      Options::KINFUNCTIONALS::TF,
      0.0,
      {std::make_shared<EnergyComponentController>()});

  std::vector<std::shared_ptr<DensityMatrixController<SPIN> > > envDens
  = {env->getElectronicStructure<SPIN>()->getDensityMatrixController()};
  NAddFuncPotential<SPIN> naddFuncPot(
      act,
      act->getElectronicStructure<SPIN>()->getDensityMatrixController(),
      envDens,
      supersystemGrid,
      FunctionalClassResolver::resolveFunctional(Options::KINFUNCTIONALS::TF));

  auto fHuz = huzFDEPotential.getMatrix();
  auto fNadd=naddFuncPot.getMatrix();

  EXPECT_NEAR((fHuz-fNadd).sum(),0.0,1e-7);

  SystemController__TEST_SUPPLY::cleanUp();
}

} /*namespace Serenity*/
