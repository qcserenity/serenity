/**
 * @file OrbitalPair_test.cpp
 *
 * @date May 20, 2019
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
#include "data/OrbitalPair.h"
#include "data/SingleSubstitution.h"
#include "io/HDF5.h"
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h" //k-sets.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"       //kl-sets.
/* Include Std and External Headers */
#include <gtest/gtest.h>
#include <memory> //smrt_ptr

namespace Serenity {
class OrbitalPairTest : public ::testing::Test {
 protected:
  OrbitalPairTest() {
    _testPair = std::make_shared<OrbitalPair>(0, 0, 3.33e-7, 1e-4, 1e-6);
    Eigen::MatrixXd dummyMatrix(2, 2);
    dummyMatrix << 0.6, 0.8, 1.7, 0.003;
    _testPair->toPAODomain = 0.06 * dummyMatrix;
    _testPair->k_ij = dummyMatrix;
    _testPair->t_ij = (-dummyMatrix.array() / 1.2).matrix();
    _testPair->ac_bd.reset(new Matrix<Eigen::MatrixXd>(2, 2, dummyMatrix));
    (*_testPair->ac_bd)(0, 0) = -dummyMatrix * dummyMatrix.transpose(); // Is expected to be symmetric!
    (*_testPair->ac_bd)(1, 0) = -1.2 * dummyMatrix;
    _testPair->jc_ab.push_back(dummyMatrix * dummyMatrix);
    _testPair->jc_ab.push_back(dummyMatrix);
    _testPair->ic_ab.push_back(dummyMatrix * dummyMatrix);
    _testPair->ic_ab.push_back(dummyMatrix);
    auto ikPair = std::make_shared<OrbitalPair>(1, 0, 3.33e-7, 1e-4, 1e-6);
    auto kkPair = std::make_shared<OrbitalPair>(1, 1, 3.33e-7, 1e-4, 1e-6);
    auto zzPair = std::make_shared<OrbitalPair>(0, 0, 3.33e-7, 1e-4, 1e-6);
    auto singles_k = std::make_shared<SingleSubstitution>(kkPair, 0.01);
    auto singles_0 = std::make_shared<SingleSubstitution>(zzPair, 0.01);
    ikPair->singles_i = singles_0;
    ikPair->singles_j = singles_k;
    _testPair->singles_i = singles_0;
    _testPair->singles_j = singles_0;
    auto kSet = std::make_shared<CouplingOrbitalSet>(_testPair, ikPair, ikPair, 1);
    _testPair->coupledPairs.push_back(kSet);

    Eigen::MatrixXd dummyMatrix2 = dummyMatrix.transpose() * dummyMatrix;
    auto klPairSet = std::make_shared<KLOrbitalSet>(kkPair, kkPair);
    _testPair->klPairSets.push_back(klPairSet);
    Eigen::MatrixXd dummyMatrixK = 16.0 * dummyMatrix;

    _testPair->coupledPairs[0]->ka_bc.push_back(dummyMatrix);
    _testPair->coupledPairs[0]->ka_bc.push_back(dummyMatrix2);
    _testPair->coupledPairs[0]->ab_kcX2_M_ak_bc.push_back(0.1 * dummyMatrix);
    _testPair->coupledPairs[0]->ab_kcX2_M_ak_bc.push_back(0.8 * dummyMatrix2);
    _testPair->ij_ab = 0.1 * dummyMatrix;
    _testPair->coupledPairs[0]->ia_kc = 9.6 * dummyMatrix;
    _testPair->coupledPairs[0]->ja_kc = 0.587 * dummyMatrix;
    _testPair->coupledPairs[0]->ik_ca = dummyMatrix;
    _testPair->coupledPairs[0]->jk_ca = 0.7 * dummyMatrix;
    _testPair->ia_bc.push_back(dummyMatrix2);
    _testPair->ia_bc.push_back(dummyMatrix);
    _testPair->ja_bc.push_back(dummyMatrix2);
    _testPair->ja_bc.push_back(dummyMatrix);
    _testPair->iaS_jbSX2_M_ij_aSbS = dummyMatrix;
  }

  virtual ~OrbitalPairTest() = default;

  std::shared_ptr<OrbitalPair> _testPair;
};

/**
 * @test Tests the writing and loading of pair integrals from disk.
 */
TEST_F(OrbitalPairTest, writeAndRead_minimum) {
  // Save the initial integral values.
  Eigen::MatrixXd ac_bd_00 = (*_testPair->ac_bd)(0, 0);
  Eigen::MatrixXd ac_bd_10 = (*_testPair->ac_bd)(0, 1);
  Eigen::VectorXd ja_ik_0 = _testPair->coupledPairs[0]->ja_ik;
  Eigen::MatrixXd ia_jk_0 = _testPair->coupledPairs[0]->ia_jk;
  Eigen::MatrixXd ij_ak_0 = _testPair->coupledPairs[0]->ij_ak;
  Eigen::MatrixXd ia_kc_0 = _testPair->coupledPairs[0]->ia_kc;
  Eigen::MatrixXd ik_ca_0 = _testPair->coupledPairs[0]->ik_ca;
  Eigen::MatrixXd ja_kc_0 = _testPair->coupledPairs[0]->ja_kc;
  Eigen::MatrixXd jk_ca_0 = _testPair->coupledPairs[0]->jk_ca;

  // Write the integrals to a file and flush everything in memory.
  std::string fileName = "PairIntegrals.h5";
  HDF5::H5File file(fileName.c_str(), H5F_ACC_TRUNC);
  _testPair->writeIntegralsToFile(file);
  file.close();
  _testPair->flushIntegrals();
  EXPECT_TRUE(!_testPair->ac_bd);
  EXPECT_EQ((unsigned int)_testPair->coupledPairs[0]->ia_kc.size(), 0);
  EXPECT_EQ((unsigned int)_testPair->coupledPairs[0]->ik_ca.size(), 0);
  EXPECT_EQ((unsigned int)_testPair->coupledPairs[0]->ja_kc.size(), 0);
  EXPECT_EQ((unsigned int)_testPair->coupledPairs[0]->jk_ca.size(), 0);
  // Load the integrals from disk and compare them to the initial integrals.
  HDF5::Filepath name(fileName);
  HDF5::H5File readFile(name.c_str(), H5F_ACC_RDONLY);
  _testPair->loadIntegralsFromFile(readFile);
  readFile.close();
  Eigen::MatrixXd diff = (*_testPair->ac_bd)(0, 0) - ac_bd_00;
  EXPECT_NEAR((double)diff.array().abs().sum(), 0.0, 1e-12);
  diff = (*_testPair->ac_bd)(0, 1) - ac_bd_10;
  EXPECT_NEAR((double)diff.array().abs().sum(), 0.0, 1e-12);
  diff = ia_kc_0 - _testPair->coupledPairs[0]->ia_kc;
  EXPECT_NEAR((double)diff.array().abs().sum(), 0.0, 1e-12);
  diff = ik_ca_0 - _testPair->coupledPairs[0]->ik_ca;
  EXPECT_NEAR((double)diff.array().abs().sum(), 0.0, 1e-12);
  diff = ja_kc_0 - _testPair->coupledPairs[0]->ja_kc;
  EXPECT_NEAR((double)diff.array().abs().sum(), 0.0, 1e-12);
  diff = jk_ca_0 - _testPair->coupledPairs[0]->jk_ca;
  EXPECT_NEAR((double)diff.array().abs().sum(), 0.0, 1e-12);
  std::remove(fileName.c_str());
}

} /* namespace Serenity */
