/**
 * @file NROCalculator_test.cpp
 *
 * @date Oct 13, 2021
 * @author Anton Rikus
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
#include "postHF/LRSCF/Analysis/NROCalculator.h"
#include "io/HDF5.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
#include "tasks/PlotTask.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

/**
 * @class NROCalculatorTest
 * @brief Sets everything up for the tests of NROCalculator.h/.cpp .
 */
class NROCalculatorTest : public ::testing::Test {
 protected:
  NROCalculatorTest() {
  }

  virtual ~NROCalculatorTest() = default;

  static void TearDownTestCase() {
    SystemController__TEST_SUPPLY::cleanUp();
  }
};
TEST_F(NROCalculatorTest, NROs) {
  auto systemController =
      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::Formaldehyde_HF_AUG_CC_PVDZ, true);
  std::vector<std::shared_ptr<SystemController>> active = {systemController};
  LRSCFTask<Options::SCF_MODES::RESTRICTED> lrscf(active);
  auto task = PlotTask<RESTRICTED>({systemController}, {});
  lrscf.settings.densFitK = Options::DENS_FITS::RI;
  lrscf.settings.frequencies = {1.0};
  lrscf.run();
  task.settings.nros = true;
  task.run();
  Eigen::MatrixXd expecSingularValues = Eigen::MatrixXd::Zero(8, 3);
  expecSingularValues.row(0) << 3.4037098170495778e-01, 3.3876373553641048e-01, 4.0339674414029036e-01;
  expecSingularValues.row(1) << 2.7101060750640849e-01, 2.1576636067542329e-01, 1.9878081795264674e-01;
  expecSingularValues.row(2) << 1.8296561013728360e-01, 1.9309857097168870e-01, 1.5528983405325902e-01;
  expecSingularValues.row(3) << 1.0254532000170059e-01, 1.2031978576498660e-01, 1.4190689056027833e-01;
  expecSingularValues.row(4) << 8.4466609493293140e-02, 1.0327941177658408e-01, 7.1094377364036690e-02;
  expecSingularValues.row(5) << 1.6856270126986050e-02, 2.7964558976879869e-02, 2.7325756981525628e-02;
  expecSingularValues.row(6) << 1.3929429454001467e-03, 5.0897558009538118e-04, 1.8429345057537230e-03;
  expecSingularValues.row(7) << 3.9165808397019855e-04, 2.9860071793164826e-04, 3.6264444220940870e-04;

  LRSCFTaskSettings lrscfSettings;
  lrscfSettings.loadType = Options::LRSCF_TYPE::ISOLATED;
  auto lrscfcontroller = std::make_shared<LRSCFController<RESTRICTED>>(systemController, lrscfSettings);
  std::vector<Eigen::MatrixXd> XY(2);
  std::string filename = systemController->getSystemPath() + systemController->getSystemName() + "_lrscf_resp.iso.res.h5";
  HDF5::H5File afile(filename, H5F_ACC_RDONLY);
  HDF5::dataset_exists(afile, "X+Y");
  HDF5::dataset_exists(afile, "X-Y");
  HDF5::dataset_exists(afile, "frequencies");

  HDF5::load(afile, "X+Y", XY[0]);
  HDF5::load(afile, "X-Y", XY[1]);
  afile.close();
  NROCalculator<RESTRICTED> nro(XY, lrscfcontroller);
  nro.getNROs(0);

  for (unsigned iocc = 0; iocc < 8; iocc++) {
    for (unsigned j = 0; j < 3; j++) {
      EXPECT_NEAR(expecSingularValues(iocc, j), nro.getSingularValues()[0](iocc, j), 1.0e-6);
    }
  }
  std::vector<std::string> SpatDirection = {"x", "y", "z"};
  for (unsigned iFreq = 0; iFreq < XY[0].cols() / 3; iFreq++) {
    for (unsigned row = 0; row < 3; row++) {
      double accSingularValues = 0.0;
      for (unsigned i = 0; i < expecSingularValues.rows(); i++) {
        if (accSingularValues < task.settings.nrominimum) {
          std::string partfilename = (systemController->getSystemPath() + "freq_" + std::to_string(iFreq + 1) +
                                      SpatDirection[row] + "particleNRO_" + std::to_string(i + 1) + ".cube");
          std::string holefilename = (systemController->getSystemPath() + "freq_" + std::to_string(iFreq + 1) +
                                      SpatDirection[row] + "holeNRO_" + std::to_string(i + 1) + ".cube");
          EXPECT_EQ(0, std::remove(partfilename.c_str()));
          EXPECT_EQ(0, std::remove(holefilename.c_str()));
          accSingularValues += nro.getSingularValues()[iFreq](i, row);
        }
      }
    }
  }
}

} /* namespace Serenity */
