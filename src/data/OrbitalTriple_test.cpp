/**
 * @file OrbitalTriple_test.cpp
 *
 * @author Moritz Bensberg
 * @date Nov 6, 2019
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
#include "data/OrbitalTriple.h"                                     //To be tested.
#include "analysis/PAOSelection/TNOConstructor.h"                   //TNO construction.
#include "basis/BasisController.h"                                  //BasisController.
#include "data/ElectronicStructure.h"                               //Access to orbital coefficients.
#include "data/OrbitalController.h"                                 //Access to orbital coefficients.
#include "data/PAOController.h"                                     //Access to PAOs.
#include "energies/EnergyContributions.h"                           //Get the canonical triples correction.
#include "integrals/MO3CenterIntegralController.h"                  //Integrals over MOs and PAOs.
#include "integrals/looper/TwoElecFourCenterIntLooper.h"            //References.
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h" //Calculate two center aux.-metric.
#include "integrals/transformer/Ao2MoTransformer.h"                 //References.
#include "io/FormattedOutputStream.h"                               //Print level.
#include "math/RegularRankFourTensor.h"                             //References.
#include "postHF/CC/DLPNO_CCSD.h"                                   //Initial doubles amplitudes.
#include "postHF/LocalCorrelation/LocalCorrelationController.h"     //Access to PAOs.
#include "system/SystemController.h"                                //Test systems.
#include "tasks/CoupledClusterTask.h"                               //Cross-check for triples correction.
#include "tasks/ScfTask.h"                                          //Clean electronic structure.
#include "testsupply/SystemController__TEST_SUPPLY.h"               //Test systems.
/* Include Std and External Headers */
#include <gtest/gtest.h> //gtest
#include <memory>        //smrt_ptr

namespace Serenity {
class OrbitalTripleTest : public ::testing::Test {
 protected:
  OrbitalTripleTest() {
  }

  virtual ~OrbitalTripleTest() = default;
};
TEST_F(OrbitalTripleTest, calculateIntegralsAndEnergy) {
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::DEBUGGING;
  auto system = SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2Dimer_Def2_TZVP, true);
  auto basisController = system->getBasisController();
  auto nBasisFunc = basisController->getNBasisFunctions();
  ScfTask<Options::SCF_MODES::RESTRICTED> scf(system);
  scf.run();
  LocalCorrelationSettings lCSettings;
  lCSettings.pnoThreshold = 1e-14;
  lCSettings.doiPAOThreshold = 1e-14;
  lCSettings.mullikenThreshold = 0.0;
  lCSettings.orbitalToShellThreshold = 0.0;
  lCSettings.dumpIntegrals = false;
  auto localCorrelationController = std::make_shared<LocalCorrelationController>(system, lCSettings);
  DLPNO_CCSD localCCSD(localCorrelationController, 1e-5, 100);
  localCCSD.calculateElectronicEnergyCorrections().sum();

  localCorrelationController->getOrbitalTripleSets();
  std::vector<std::shared_ptr<OrbitalTriple>> triples = localCorrelationController->getOrbitalTriples();

  auto activeSystem = localCorrelationController->getActiveSystemController();
  auto auxBasisController = activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
  auto tnoConstructor = localCorrelationController->produceTNOConstructor();
  auto mo3CenterIntegralController = localCorrelationController->getMO3CenterIntegralController();
  // Make sure that the integrals are available in order to avoid any parallel construction of
  // the integral sets, since these functions are not thread save.
  Eigen::SparseVector<int> auxSuperDomain = Eigen::VectorXi::Constant(auxBasisController->getBasis().size(), 1).sparseView();
  const auto& iaK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain);
  const auto& abK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ab_K, auxSuperDomain);
  const auto& klK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::kl_K, auxSuperDomain);

  MatrixInBasis<RESTRICTED> metric(auxBasisController);
  Ao2MoExchangeIntegralTransformer::calculateTwoCenterIntegrals(metric);

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, basisController->getPrescreeningThreshold());
  RegularRankFourTensor<double> eris(nBasisFunc);
  auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                 const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId; // no warning, please
    eris(b, a, i, j) = integral(0);
    eris(b, a, j, i) = integral(0);
    eris(a, b, j, i) = integral(0);
    eris(a, b, i, j) = integral(0);
    eris(i, j, b, a) = integral(0);
    eris(i, j, a, b) = integral(0);
    eris(j, i, b, a) = integral(0);
    eris(j, i, a, b) = integral(0);
  };
  looper.loop(storeERIS);
  double triplesCorrection = 0.0;
  for (const auto& triple : triples) {
    tnoConstructor->transformToTNOBasis(triple);
    triple->calculateIntegrals(auxBasisController, mo3CenterIntegralController, iaK, abK, klK, metric);

    // Reference:
    Eigen::MatrixXd R_00 = localCorrelationController->getPAOController()->getAllPAOs();

    CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeffs =
        system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
    const Eigen::MatrixXd tnoCoefficients = R_00 * triple->getTNOCoefficients();
    OutputControl::dOut << tnoCoefficients.rows() << "x" << tnoCoefficients.cols() << " NBasis " << nBasisFunc << std::endl;
    coeffs.block(0, 2, nBasisFunc, nBasisFunc - 2) = tnoCoefficients;
    // Calculate fully transformed integrals
    Ao2MoTransformer aoToMo(basisController);
    RegularRankFourTensor<double> eris_copy = eris;
    RegularRankFourTensor<double> result(nBasisFunc, 0.0);
    aoToMo.transformTwoElectronIntegrals(eris_copy, result, coeffs, nBasisFunc);

    unsigned int i = triple->getI();
    unsigned int j = triple->getJ();
    unsigned int k = triple->getK();
    OutputControl::dOut << "Triplet: " << i << " " << j << " " << k << std::endl;
    OutputControl::dOut << "IL: ";
    for (const auto l : triple->getILIndices())
      OutputControl::dOut << l << " ";
    OutputControl::dOut << std::endl << "JL: ";
    for (const auto l : triple->getJLIndices())
      OutputControl::dOut << l << " ";
    OutputControl::dOut << std::endl << "KL: ";
    for (const auto l : triple->getKLIndices())
      OutputControl::dOut << l << " ";
    OutputControl::dOut << std::endl;
    unsigned int nTNOs = tnoCoefficients.cols();
    const auto& ia_bc = triple->get_ia_bc();
    const auto& ja_bc = triple->get_ja_bc();
    const auto& ka_bc = triple->get_ka_bc();

    const auto& ja_il = triple->get_ja_il();
    const auto& ia_jl = triple->get_ia_jl();
    const auto& ka_il = triple->get_ka_il();
    const auto& ia_kl = triple->get_ia_kl();
    const auto& ka_jl = triple->get_ka_jl();
    const auto& ja_kl = triple->get_ja_kl();

    const auto& ia_jb = triple->get_ia_jb();
    const auto& ia_kb = triple->get_ia_kb();
    const auto& ja_kb = triple->get_ja_kb();

    // The more virtual coefficients, the worse the RI approximation.
    for (unsigned int c = 0; c < nTNOs; ++c) {
      for (unsigned int b = 0; b < nTNOs; ++b) {
        for (unsigned int a = 0; a < nTNOs; ++a) {
          double diff = ia_bc[c](a, b) - result(i, a + 2, b + 2, c + 2);
          EXPECT_NEAR(0.0, diff, 5e-3);
          diff = ja_bc[c](a, b) - result(j, a + 2, b + 2, c + 2);
          EXPECT_NEAR(0.0, diff, 5e-3);
          diff = ka_bc[c](a, b) - result(k, a + 2, b + 2, c + 2);
          EXPECT_NEAR(0.0, diff, 5e-3);
        } // for a
      }   // for b
    }     // for c

    for (unsigned int l = 0; l < 2; ++l) {
      for (unsigned int a = 0; a < nTNOs; ++a) {
        double diff = ja_il(a, l) - result(j, a + 2, i, l);
        EXPECT_NEAR(0.0, diff, 1e-4);
        diff = ia_jl(a, l) - result(i, a + 2, j, l);
        EXPECT_NEAR(0.0, diff, 1e-4);
        diff = ka_il(a, l) - result(k, a + 2, i, l);
        EXPECT_NEAR(0.0, diff, 1e-4);
        diff = ia_kl(a, l) - result(i, a + 2, k, l);
        EXPECT_NEAR(0.0, diff, 1e-4);
        diff = ka_jl(a, l) - result(k, a + 2, j, l);
        EXPECT_NEAR(0.0, diff, 1e-4);
        diff = ja_kl(a, l) - result(j, a + 2, k, l);
        EXPECT_NEAR(0.0, diff, 1e-4);
      } // for a
    }   // for l

    for (unsigned int b = 0; b < nTNOs; ++b) {
      for (unsigned int a = 0; a < nTNOs; ++a) {
        double diff = ia_jb(a, b) - result(i, a + 2, j, b + 2);
        EXPECT_NEAR(0.0, diff, 1.1e-3);
        diff = ia_kb(a, b) - result(i, a + 2, k, b + 2);
        EXPECT_NEAR(0.0, diff, 7e-4);
        diff = ja_kb(a, b) - result(j, a + 2, k, b + 2);
        EXPECT_NEAR(0.0, diff, 7e-4);
      } // for a
    }   // for b

    // Check triplet energy increment.
    double tripletEnergy = triple->getTripleEnergy();
    triplesCorrection += tripletEnergy;
    OutputControl::dOut << tripletEnergy << std::endl;
  } // for triple

  CoupledClusterTask task(system);
  task.settings.level = Options::CC_LEVEL::CCSD_T;
  task.run();
  double canonicalTriplesCorrection =
      system->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getEnergy(ENERGY_CONTRIBUTIONS::TRIPLES_CORRECTION);
  EXPECT_NEAR(canonicalTriplesCorrection, triplesCorrection, 5e-5);

  SystemController__TEST_SUPPLY::cleanUp();
  GLOBAL_PRINT_LEVEL = Options::GLOBAL_PRINT_LEVELS::NORMAL;
}

} /* namespace Serenity */
