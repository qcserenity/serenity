/**
 * @file Kernel_test.cpp
 *
 * @date Nov 09, 2017
 * @author Michael Boeckers
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
#include "testsupply/GridController__TEST_SUPPLY.h"
#include "postHF/LRSCF/Kernel.h"
#include "testsupply/SystemController__TEST_SUPPLY.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>

namespace Serenity {

///**
// * @class KernelTest
// * @brief Sets everything up for the tests of Kernel.h/.cpp .
// */
//class KernelTest : public ::testing::Test {
// protected:
//  KernelTest() {
//  }
//
//  virtual ~KernelTest() = default;
//
//  /// system
//  static void TearDownTestCase() {
//    SystemController__TEST_SUPPLY::cleanUp();
//  }
//};
//
//TEST(KernelTest,RKERNEL) {
//  //Setup test systems
//  auto activeSystem =
//      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ACTIVE,false);
//  auto environmentSystem =
//      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT,false);
//  //Change grids to tiny grid
//  auto testGrid =
//      GridController__TEST_SUPPLY::getGridController(TEST_GRID_CONTROLLERS::TINY);
//  activeSystem->setGridController(testGrid);
//  environmentSystem->setGridController(testGrid);
//  //Build kernel
//  Kernel<Options::SCF_MODES::RESTRICTED> kernel(
//      activeSystem,{environmentSystem},false,false,
//      Options::FUNCTIONALS::PW91,Options::FUNCTIONALS::PW91K,Options::FUNCTIONALS::PW91);
//
//  //Get derivatives
//  auto d2FdRho2 = kernel.getD2F_dRho2();
//  auto dFdSigma = kernel.getDF_dSigma();
//  auto d2FdSigma2 = kernel.getD2F_dSigma2();
//  auto d2FdRhodSigma = kernel.getD2F_dRhodSigma();
//  auto totalNaddD2FdSigma2 = kernel.getTotalNaddD2F_dSigma2();
//  auto totalNaddD2FdRhodSigma = kernel.getTotalNaddD2F_dRhodSigma();
//  auto activeDensityGradient = kernel.getActiveDensityGradient();
//  auto totalDensityGradient = kernel.getTotalDensityGradient();
//
//  //Results (calculated with correct version)
//  Eigen::VectorXd cD2FdRho2(4);
//  Eigen::VectorXd cDFdSigma(4);
//  Eigen::VectorXd cD2FdSigma2(4);
//  Eigen::VectorXd cD2FdRhodSigma(4);
//  Eigen::VectorXd cTotaldNaddD2FdSigma2(4);
//  Eigen::VectorXd cTotaldNaddD2FdRhodSigma(4);
//  Eigen::VectorXd cActiveGradientX(4);
//  Eigen::VectorXd cActiveGradientY(4);
//  Eigen::VectorXd cActiveGradientZ(4);
//  Eigen::VectorXd cTotalGradientX(4);
//  Eigen::VectorXd cTotalGradientY(4);
//  Eigen::VectorXd cTotalGradientZ(4);
//
//  cD2FdRho2                <<   -0.88273492178765967,     -2.6983379747064409,      -3.8954725605859934,     282.67226752447915       ;
//  cDFdSigma                <<   0.012551679175773375,     -0.057231304646680192,    -1.3548039405401666,     5.6369904722738262       ;
//  cD2FdSigma2              <<   -3.842618711205402,       2.0132399903674436,       14662.020696712108,      2513721.0373972347       ;
//  cD2FdRhodSigma           <<   0.20952523975316723,      0.87827121885866521,      -154.04780993524059,     -7185.698235906937       ;
//  cTotaldNaddD2FdRhodSigma <<   -0.26281763799186497,     0.073373513561957093,     -45.01633649773629,      -43561.562363833968      ;
//  cTotaldNaddD2FdSigma2    <<   -0.11652722207780819,     -1.8983739755630518,      369.18353005626614,      2256074.8956254171       ;
//  cActiveGradientX         <<   0,                        -0.14441701679769445,     0.0063116060517008901,   -0.0010981989064538485   ;
//  cActiveGradientY         <<   0,                        0,                        -0.003155803025850445,   -0.0001830331510756414   ;
//  cActiveGradientZ         <<   2.7755575615628914e-16,   4.163336342344337e-17,    -0.002212731536127728,   0.00028521695070329802   ;
//  cTotalGradientX          <<   0,                        -0.15057756675336337,     0.011954844341751782,    -0.0011550166089871474   ;
//  cTotalGradientY          <<   0,                        0,                        -0.0059774221708758912,  -0.00019250276816452455  ;
//  cTotalGradientZ          <<   0.024345754211872385,     0.011523563674304498,     0.00048025708807551519,  0.00034042313355003265   ;
//
//  for (unsigned int i = 0; i < 4; ++i) {
//    EXPECT_NEAR((*d2FdRho2)(i)/cD2FdRho2(i),   1.0             ,1.0e-5);
//    EXPECT_NEAR((*dFdSigma)(i)/cDFdSigma(i),    1.0           ,1.0e-5);
//    EXPECT_NEAR((*d2FdSigma2)(i)/cD2FdSigma2(i),  1.0            ,1.0e-5);
//    EXPECT_NEAR((*d2FdRhodSigma)(i)/cD2FdRhodSigma(i),  1.0         ,1.0e-5);
//    EXPECT_NEAR((*totalNaddD2FdRhodSigma)(i)/cTotaldNaddD2FdRhodSigma(i),1.0 ,1.0e-5);
//    EXPECT_NEAR((*totalNaddD2FdSigma2)(i)/cTotaldNaddD2FdSigma2(i), 1.0   ,1.0e-5);
//    EXPECT_NEAR((*activeDensityGradient).x(i),   cActiveGradientX(i)         ,1.0e-5);
//    EXPECT_NEAR((*activeDensityGradient).y(i),   cActiveGradientY(i)         ,1.0e-5);
//    EXPECT_NEAR((*activeDensityGradient).z(i),   cActiveGradientZ(i)         ,1.0e-5);
//    EXPECT_NEAR((*totalDensityGradient).x(i),    cTotalGradientX(i)          ,1.0e-5);
//    EXPECT_NEAR((*totalDensityGradient).y(i),    cTotalGradientY(i)          ,1.0e-5);
//    EXPECT_NEAR((*totalDensityGradient).z(i),    cTotalGradientZ(i)          ,1.0e-5);
//  }
//}
//
//TEST(KernelTest,UKERNEL) {
//  //Setup test systems
//  auto activeSystem =
//      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::OH_MINBAS_PBE,false);
//  auto environmentSystem =
//      SystemController__TEST_SUPPLY::getSystemController(TEST_SYSTEM_CONTROLLERS::H2_MINBAS_ENVIRONMENT,false);
//  //Change grids to tiny grid
//  auto testGrid =
//      GridController__TEST_SUPPLY::getGridController(TEST_GRID_CONTROLLERS::TINY);
//  activeSystem->setGridController(testGrid);
//  environmentSystem->setGridController(testGrid);
//  //Build kernel
//  Kernel<Options::SCF_MODES::UNRESTRICTED> kernel(
//      activeSystem,{environmentSystem},false,false,
//      Options::FUNCTIONALS::PBE,Options::FUNCTIONALS::PW91K,Options::FUNCTIONALS::PW91);
//
//  //Get derivatives
//  auto d2FdRho2 = kernel.getD2F_dRho2();
//  auto dFdSigma = kernel.getDF_dSigma();
//  auto d2FdSigma2 = kernel.getD2F_dSigma2();
//  auto d2FdRhodSigma = kernel.getD2F_dRhodSigma();
//  auto totalNaddD2FdSigma2 = kernel.getTotalNaddD2F_dSigma2();
//  auto totalNaddD2FdRhodSigma = kernel.getTotalNaddD2F_dRhodSigma();
//  auto activeDensityGradient = kernel.getActiveDensityGradient();
//  auto totalDensityGradient = kernel.getTotalDensityGradient();
//
//  Eigen::VectorXd cD2FdRho2_aa(4);
//  Eigen::VectorXd cD2FdRho2_ab(4);
//  Eigen::VectorXd cD2FdRho2_bb(4);
//  Eigen::VectorXd cDFdSigma_aa(4);
//  Eigen::VectorXd cDFdSigma_ab(4);
//  Eigen::VectorXd cDFdSigma_bb(4);
//  Eigen::VectorXd cD2FdSigma2_aaaa(4);
//  Eigen::VectorXd cD2FdSigma2_aaab(4);
//  Eigen::VectorXd cD2FdSigma2_aabb(4);
//  Eigen::VectorXd cD2FdSigma2_abab(4);
//  Eigen::VectorXd cD2FdSigma2_abbb(4);
//  Eigen::VectorXd cD2FdSigma2_bbbb(4);
//  Eigen::VectorXd cD2FdRhodSigma_aaa(4);
//  Eigen::VectorXd cD2FdRhodSigma_aab(4);
//  Eigen::VectorXd cD2FdRhodSigma_abb(4);
//  Eigen::VectorXd cD2FdRhodSigma_baa(4);
//  Eigen::VectorXd cD2FdRhodSigma_bab(4);
//  Eigen::VectorXd cD2FdRhodSigma_bbb(4);
//  Eigen::VectorXd cTotalD2FdSigma2_aaaa(4);
//  Eigen::VectorXd cTotalD2FdSigma2_aaab(4);
//  Eigen::VectorXd cTotalD2FdSigma2_aabb(4);
//  Eigen::VectorXd cTotalD2FdSigma2_abab(4);
//  Eigen::VectorXd cTotalD2FdSigma2_abbb(4);
//  Eigen::VectorXd cTotalD2FdSigma2_bbbb(4);
//  Eigen::VectorXd cTotalD2FdRhodSigma_aaa(4);
//  Eigen::VectorXd cTotalD2FdRhodSigma_aab(4);
//  Eigen::VectorXd cTotalD2FdRhodSigma_abb(4);
//  Eigen::VectorXd cTotalD2FdRhodSigma_baa(4);
//  Eigen::VectorXd cTotalD2FdRhodSigma_bab(4);
//  Eigen::VectorXd cTotalD2FdRhodSigma_bbb(4);
//  Eigen::VectorXd cActiveGradient_ax(4);
//  Eigen::VectorXd cActiveGradient_ay(4);
//  Eigen::VectorXd cActiveGradient_az(4);
//  Eigen::VectorXd cActiveGradient_bx(4);
//  Eigen::VectorXd cActiveGradient_by(4);
//  Eigen::VectorXd cActiveGradient_bz(4);
//  Eigen::VectorXd cTotalGradient_ax(4);
//  Eigen::VectorXd cTotalGradient_ay(4);
//  Eigen::VectorXd cTotalGradient_az(4);
//  Eigen::VectorXd cTotalGradient_bx(4);
//  Eigen::VectorXd cTotalGradient_by(4);
//  Eigen::VectorXd cTotalGradient_bz(4);
//  cD2FdRho2_aa            <<    -0.01929460646546426,          -1.1523637841217613,           -26.774562822844089,        70.499624733445899       ;
//  cD2FdRho2_ab            <<    -0.0004336600003688839,        -0.078130675590302748,         -6.0812420737700377,        -1.0445304244669629      ;
//  cD2FdRho2_bb            <<    -0.019303892653148819,         -1.1564310745590394,           -33.137758436059165,        76.636325837494496       ;
//  cDFdSigma_aa            <<    -3.7835142037557544e-06,       -0.012162823945211341,         -1.1132167743244357,        1.5133144882380325       ;
//  cDFdSigma_ab            <<    7.5592554574350796e-06,        0.0094875246950914607,         0.12550077643073285,        0.0067987968553448203    ;
//  cDFdSigma_bb            <<    -3.7910744279978578e-06,       -0.01228974492059842,          -1.2499274732388912,        1.483499479060125        ;
//  cD2FdSigma2_aaaa        <<    -7.8281000975933434e-10,       0.015954396744621055,          2646.2719928464139,         66156076.199792035       ;
//  cD2FdSigma2_aaab        <<    6.8290066126191062e-12,        -0.00013476707972395535,       67.045566812666948,         1034.5605453037952       ;
//  cD2FdSigma2_aabb        <<    3.4145033063095531e-12,        -6.7383539861977674e-05,       33.522783406333474,         517.28027265189758       ;
//  cD2FdSigma2_abab        <<    1.3658013225238212e-11,        -0.0002695341594479107,        134.0911336253339,          2069.1210906075903       ;
//  cD2FdSigma2_abbb        <<    6.8290066126191062e-12,        -0.00013476707972395535,       67.045566812666948,         1034.5605453037952       ;
//  cD2FdSigma2_bbbb        <<    -7.8497073745926813e-10,       0.016283347371892007,          2657.5802249181506,         67766374.783372074       ;
//  cD2FdRhodSigma_aaa      <<    1.5300687614504435e-06,        0.13099388731073916,           325.58554337783426,         -324133.3479296586       ;
//  cD2FdRhodSigma_aab      <<    9.3214410553396133e-09,        -0.0010871286955383835,        -7.8215419497936267,        -12.02239520774976       ;
//  cD2FdRhodSigma_abb      <<    4.6607205276698067e-09,        -0.00054356434776919177,       -3.9107709748968134,        -6.0111976038748836      ;
//  cD2FdRhodSigma_baa      <<    4.6601409733624834e-09,        -0.00054206887753044039,       -4.030864467364669,         -6.0102514238333598      ;
//  cD2FdRhodSigma_bab      <<    9.3202819467249668e-09,        -0.0010841377550608808,        -8.0617289347293379,        -12.02050284766672       ;
//  cD2FdRhodSigma_bbb      <<    1.5323874538567207e-06,        0.13252430773125407,           465.75706652688655,         -326688.64594800753      ;
//  cTotalD2FdSigma2_aaaa   <<    7.8063716546710092e-10,        -0.012833557550978904,         2658.6669699034192,         -57196757.811327487      ;
//  cTotalD2FdSigma2_aaab   <<    -1.1573001338881965e-11,       -0.0028579356726201794,        -547.44359268114738,        -13459.942732968926      ;
//  cTotalD2FdSigma2_aabb   <<    -5.7865006694409824e-12,       -0.0014289678363100897,        -273.72179634057369,        -6729.971366484463       ;
//  cTotalD2FdSigma2_abab   <<    -2.314600267776393e-11,        -0.0057158713452403588,        -1094.8871853622948,        -26919.885465937852      ;
//  cTotalD2FdSigma2_abbb   <<    -1.1573001338881965e-11,       -0.0028579356726201794,        -547.44359268114738,        -13459.942732968926      ;
//  cTotalD2FdSigma2_bbbb   <<    7.8279856392013313e-10,        -0.013058505362955403,         3561.0329937724928,         -58333173.738859832      ;
//  cTotalD2FdRhodSigma_aaa <<    -1.4519361139267468e-06,       -0.068188092353588969,         -254.68931368800622,        266900.74622622598       ;
//  cTotalD2FdRhodSigma_aab <<    -6.1368655416495541e-08,       -0.00044999763025732088,       69.552154676785122,         151.55851755537947       ;
//  cTotalD2FdRhodSigma_abb <<    -3.068432770824777e-08,        -0.00022499881512866044,       34.776077338392561,         75.779258777689734       ;
//  cTotalD2FdRhodSigma_baa <<    -3.0680498023347187e-08,       -0.00020517974084644907,       35.478416857486032,         75.76847890746663        ;
//  cTotalD2FdRhodSigma_bab <<    -6.1360996046694373e-08,       -0.00041035948169289814,       70.956833714972063,         151.53695781493326       ;
//  cTotalD2FdRhodSigma_bbb <<    -1.4540693549677533e-06,       -0.068851278462072424,         -308.16839503912013,        267831.37274707097       ;
//  cActiveGradient_ax      <<    1.6836653925556997e-10,        -0.64522691568957358,          0.0079188515760304449,      -0.00055618858694658378  ;
//  cActiveGradient_ay      <<    5.3803954980560482e-09,        -1.3528170332466428e-10,       -0.0039594257865946522,     -9.2698097725689546e-05  ;
//  cActiveGradient_az      <<    1.440848798429184,             -0.022145961923094495,         -0.0033955986768306238,     0.00020084476645736732   ;
//  cActiveGradient_bx      <<    4.4225538571729158e-10,        -0.64009344131406576,          0.0069200257993836632,      -0.00055067657245409167  ;
//  cActiveGradient_by      <<    6.0547975628419298e-08,        -1.3074655059971187e-05,       -0.0041062351291749032,     -0.00010326286607512522  ;
//  cActiveGradient_bz      <<    1.6558956584296241,            -0.025104740311780987,         -0.0028400466308103824,     0.00020018006269221044   ;
//  cTotalGradient_ax       <<    1.6836653925556997e-10,        -0.64830719066740805,          0.010740470721055891,       -0.00058459743821323324  ;
//  cTotalGradient_ay       <<    5.3803954980560482e-09,        -1.3528170332466428e-10,       -0.005370235359107375,      -9.7432906270131123e-05  ;
//  cTotalGradient_az       <<    1.45302167553512,              -0.016384180085942265,         -0.0020491043647290022,     0.00022844785788073463   ;
//  cTotalGradient_bx       <<    4.4225538571729158e-10,        -0.64317371629190023,          0.0097416449444091089,      -0.00057908542372074113  ;
//  cTotalGradient_by       <<    6.0547975628419298e-08,        -1.3074655059971187e-05,       -0.005517044701687626,      -0.0001079976746195668   ;
//  cTotalGradient_bz       <<    1.66806853553556,              -0.01934295847462876,          -0.0014935523187087608,     0.00022778315411557776   ;
//
//
//  //Results (calculated with correct version)
//  for (unsigned int i = 0; i < 4; ++i) {
//    EXPECT_NEAR((*d2FdRho2).aa(i)                  ,  cD2FdRho2_aa(i)             ,1.0e-4);
//    EXPECT_NEAR((*d2FdRho2).ab(i)                  ,  cD2FdRho2_ab(i)             ,1.0e-4);
//    EXPECT_NEAR((*d2FdRho2).bb(i)                  ,  cD2FdRho2_bb(i)             ,1.0e-4);
//    EXPECT_NEAR((*dFdSigma).aa(i)                  ,  cDFdSigma_aa(i)             ,1.0e-4);
//    EXPECT_NEAR((*dFdSigma).ab(i)                  ,  cDFdSigma_ab(i)             ,1.0e-4);
//    EXPECT_NEAR((*dFdSigma).bb(i)                  ,  cDFdSigma_bb(i)             ,1.0e-4);
//    EXPECT_NEAR((*d2FdSigma2).aaaa(i)              ,  cD2FdSigma2_aaaa(i)         ,1.0e-3);
//    EXPECT_NEAR((*d2FdSigma2).aaab(i)              ,  cD2FdSigma2_aaab(i)         ,1.0e-4);
//    EXPECT_NEAR((*d2FdSigma2).aabb(i)              ,  cD2FdSigma2_aabb(i)         ,1.0e-4);
//    EXPECT_NEAR((*d2FdSigma2).abab(i)              ,  cD2FdSigma2_abab(i)         ,1.0e-4);
//    EXPECT_NEAR((*d2FdSigma2).abbb(i)              ,  cD2FdSigma2_abbb(i)         ,1.0e-4);
//    EXPECT_NEAR((*d2FdSigma2).bbbb(i)              ,  cD2FdSigma2_bbbb(i)         ,1.0e-3);
//    EXPECT_NEAR((*d2FdRhodSigma).aaa(i)            ,  cD2FdRhodSigma_aaa(i)       ,1.0e-4);
//    EXPECT_NEAR((*d2FdRhodSigma).aab(i)            ,  cD2FdRhodSigma_aab(i)       ,1.0e-4);
//    EXPECT_NEAR((*d2FdRhodSigma).abb(i)            ,  cD2FdRhodSigma_abb(i)       ,1.0e-4);
//    EXPECT_NEAR((*d2FdRhodSigma).baa(i)            ,  cD2FdRhodSigma_baa(i)       ,1.0e-4);
//    EXPECT_NEAR((*d2FdRhodSigma).bab(i)            ,  cD2FdRhodSigma_bab(i)       ,1.0e-4);
//    EXPECT_NEAR((*d2FdRhodSigma).bbb(i)            ,  cD2FdRhodSigma_bbb(i)       ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdSigma2).aaaa(i)     ,  cTotalD2FdSigma2_aaaa(i)   ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdSigma2).aaab(i)     ,  cTotalD2FdSigma2_aaab(i)    ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdSigma2).aabb(i)     ,  cTotalD2FdSigma2_aabb(i)    ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdSigma2).abab(i)     ,  cTotalD2FdSigma2_abab(i)    ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdSigma2).abbb(i)     ,  cTotalD2FdSigma2_abbb(i)    ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdSigma2).bbbb(i)     ,  cTotalD2FdSigma2_bbbb(i)    ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdRhodSigma).aaa(i)   ,  cTotalD2FdRhodSigma_aaa(i)  ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdRhodSigma).aab(i)   ,  cTotalD2FdRhodSigma_aab(i)  ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdRhodSigma).abb(i)   ,  cTotalD2FdRhodSigma_abb(i)  ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdRhodSigma).baa(i)   ,  cTotalD2FdRhodSigma_baa(i)  ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdRhodSigma).bab(i)   ,  cTotalD2FdRhodSigma_bab(i)  ,1.0e-4);
//    EXPECT_NEAR((*totalNaddD2FdRhodSigma).bbb(i)   ,  cTotalD2FdRhodSigma_bbb(i)  ,1.0e-4);
//    EXPECT_NEAR((*activeDensityGradient).x.alpha(i),  cActiveGradient_ax(i)       ,1.0e-4);
//    EXPECT_NEAR((*activeDensityGradient).y.alpha(i),  cActiveGradient_ay(i)       ,1.0e-4);
//    EXPECT_NEAR((*activeDensityGradient).z.alpha(i),  cActiveGradient_az(i)       ,1.0e-4);
//    EXPECT_NEAR((*activeDensityGradient).x.beta(i) ,  cActiveGradient_bx(i)       ,1.0e-4);
//    EXPECT_NEAR((*activeDensityGradient).y.beta(i) ,  cActiveGradient_by(i)       ,1.0e-4);
//    EXPECT_NEAR((*activeDensityGradient).z.beta(i) ,  cActiveGradient_bz(i)       ,1.0e-4);
//    EXPECT_NEAR((*totalDensityGradient).x.alpha(i) ,  cTotalGradient_ax(i)        ,1.0e-4);
//    EXPECT_NEAR((*totalDensityGradient).y.alpha(i) ,  cTotalGradient_ay(i)        ,1.0e-4);
//    EXPECT_NEAR((*totalDensityGradient).z.alpha(i) ,  cTotalGradient_az(i)        ,1.0e-4);
//    EXPECT_NEAR((*totalDensityGradient).x.beta(i)  ,  cTotalGradient_bx(i)        ,1.0e-4);
//    EXPECT_NEAR((*totalDensityGradient).y.beta(i)  ,  cTotalGradient_by(i)        ,1.0e-4);
//    EXPECT_NEAR((*totalDensityGradient).z.beta(i)  ,  cTotalGradient_bz(i)        ,1.0e-4);
//  }
//}
}
