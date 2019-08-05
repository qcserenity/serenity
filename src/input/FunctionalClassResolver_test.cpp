/**
 * @file   FunctionalClassResolver_test.cpp
 * @author David Schnieders
 *
 * @date   Aug 2, 2018
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
#include "input/FunctionalClassResolver.h"
#include "dft/Functional.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <gtest/gtest.h>


namespace Serenity {


class FunctionalClassResolverTest : public ::testing::Test {
protected:
  FunctionalClassResolverTest(){

  }
  virtual ~FunctionalClassResolverTest() = default;

};
/**
 * @test
 * @brief Tests XC LDA.
 */
TEST_F(FunctionalClassResolverTest, LDA) {


  for(unsigned int i=1; i<=5; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no errors.
   */

}

/**
 * @test
 * @brief Tests XC GGA.
 */
TEST_F(FunctionalClassResolverTest, GGA) {


  for(unsigned int i=101; i<=111; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no errors.
   */

}

/**
 * @test
 * @brief Tests XC Hybrid.
 */
TEST_F(FunctionalClassResolverTest, Hybrid) {


  for(unsigned int i=301; i<=306; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no errors.
   */

}

/**
 * @test
 * @brief Tests XC Range Separated Hybrid.
 */
TEST_F(FunctionalClassResolverTest, RSHybrid) {


  for(unsigned int i=501; i<=502; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no errors.
   */

}

/**
 * @test
 * @brief Tests XC Double Hybrid.
 */
TEST_F(FunctionalClassResolverTest, DoubleHybrid) {


  for(unsigned int i=601; i<=606; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  for(unsigned int i=609; i<=614; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no errors.
   */

}

/**
 * @test
 * @brief Tests XC Model Potentials.
 */
TEST_F(FunctionalClassResolverTest, Model) {


  for(unsigned int i=801; i<=801; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no errors.
   */

}

/**
 * @test
 * @brief Tests KIN LDA.
 */
TEST_F(FunctionalClassResolverTest, KinLDA) {


  for(unsigned int i=1001; i<=1001; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no errors.
   */

}

/**
 * @test
 * @brief Tests KIN GGA.
 */
TEST_F(FunctionalClassResolverTest, KinGGA) {


  for(unsigned int i=1101; i<=1109; i++){
    FunctionalClassResolver::resolveFunctional(static_cast<Options::FUNCTIONALS>(i));
  }

  /*
   * Nothing much to expect here, rather than no runtime errors.
   */

}

}
