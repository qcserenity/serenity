/**
 * @file   Normalization.h
 *
 * @date   21. Juli 2013, 22:17
 * @author Thomas Dresselhaus
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
#ifndef NORMALIZATION_H
#define NORMALIZATION_H
/* Include Serenity Internal Headers */
#include "math/IntegerMaths.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <cmath>
#include <vector>

namespace Serenity {
/**
 * @class Normalization Normalization.h
 * @brief All functions needed for normalization of basis functions are in here. Purely static.
 */
class Normalization {
  Normalization() = delete;

 public:
  /*
   * This function is switched off because normalization is always done in two
   * steps. (Computational efficiency reasons)
   */
  //    static double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
  //                    double alpha1, const double* A) {
  //        const static unsigned int MAXFAC = 100;
  //        static vector<double> df(2*MAXFAC);
  //        df[0] = 1ul;
  //        df[1] = 1ul;
  //        df[2] = 1ul;
  //        for(unsigned long i=3; i<MAXFAC*2; i++) {
  //          df[i] = (i-1)*df[i-2];
  //        }
  //        return pow(2*alpha1/M_PI,0.75)*pow(4*alpha1,0.5*(l1+m1+n1))/sqrt(df[2*l1]*df[2*m1]*df[2*n1]);
  //    }
  /**
   * @brief This is the first step of the normalization.
   *
   * The radial normalization is done here (needs the gaussian exponent) and the part of the
   * angular normalization which is common to a whole shell of functions (e.g. all d-functions).
   * Must be used for each primitive inside a shell but is the same for varying magnetic quantum
   * numbers.
   *
   * @param amTotal the angular momentum of the basis function, e.g. 1 for all p-type functions.
   * @param alpha the exponent of a primitive gaussian function
   * @returns a normalization factor, which still needs to be combined with one from the
   *          finalNormalization().
   */
  static double normalizeTotalAngMom(unsigned int amTotal, double alpha) {
    return pow(2.0 * alpha / M_PI, 0.75) * pow(4.0 * alpha, 0.5 * (amTotal));
  }
  /**
   * @brief This is the second step of the normalization
   *
   * which considers the actual orientation of a function (i.e. the magnetic quantum number, e.g.
   * d_x² and d_xy differ in this). It can be done after the integration because no information of
   * the primitives is needed here.
   *
   * @param amX, amY, amZ the x y and z parts of the angular momentum, e.g. (1,0,0) for a p_x
   *                      function or (1,0,1) for a d_xz function.
   * @returns a normalization factor which needs to be combined with the result of
   *          normalizeTotalAngMom
   */
  inline static long double finalNormalization(unsigned int amX, unsigned int amY, unsigned int amZ) {
    long double tmp = std::pow((long double)(double_factorial(2 * amX - 1)) * (long double)(double_factorial(2 * amY - 1)) *
                                   (long double)(double_factorial(2 * amZ - 1)),
                               -0.5);
    return tmp;
  }

 public:
  /**
   * This does the finalNormalization() for a whole
   * shell set of 1-el-integrals, which are computed together anyway. The
   * ordering of the integrals is properly considered.
   *
   * @param integrals these are multiplied by the resulting normalization factors.
   * @param amA, amB  the angular momenta of the basis functions (e.g 0 for s and 2 for d).
   */
  static void normalizeShell(Eigen::Ref<Eigen::VectorXd> integrals, unsigned int amA, unsigned int amB) {
    const long double undo = std::pow(double_factorial(2 * amA - 1) * double_factorial(2 * amB - 1), 0.5);
    for (int aX = amA, index = 0; aX >= 0; --aX) {
      for (int aY = amA - aX; aY >= 0; --aY) {
        int aZ = amA - aY - aX;
        for (int bX = amB; bX >= 0; --bX) {
          for (int bY = amB - bX; bY >= 0; --bY) {
            int bZ = amB - bY - bX;
            integrals[index] *= (long double)(finalNormalization(aX, aY, aZ) * finalNormalization(bX, bY, bZ) * undo);
            ++index;
          }
        }
      }
    }
  }

  /**
   * This does the finalNormalization() (see above) for a whole
   * shell set of 3-center-2-el-integrals
   * TODO: Not sure if it can be done this way
   *
   * See normalizeShell(vector<double>*, unsigned int, unsigned int).
   *
   * @param integrals
   * @param amA, amB, amC
   */
  static void normalizeShell(Eigen::Ref<Eigen::VectorXd> integrals, unsigned int amA, unsigned int amB, unsigned int amC) {
    const long double undo =
        std::pow(double_factorial(2 * amA - 1) * double_factorial(2 * amB - 1) * double_factorial(2 * amC - 1), 0.5);
    for (int aX = amA, index = 0; aX >= 0; --aX) {
      for (int aY = amA - aX; aY >= 0; --aY) {
        int aZ = amA - aY - aX;
        const long double normA = finalNormalization(aX, aY, aZ);
        for (int bX = amB; bX >= 0; --bX) {
          for (int bY = amB - bX; bY >= 0; --bY) {
            int bZ = amB - bY - bX;
            const long double normB = normA * finalNormalization(bX, bY, bZ);
            for (int cX = amC; cX >= 0; --cX) {
              for (int cY = amC - cX; cY >= 0; --cY) {
                int cZ = amC - cY - cX;
                integrals[index] *= (long double)(normB * finalNormalization(cX, cY, cZ) * undo);
                ++index;
              }
            }
          }
        }
      }
    }
  };
  /**
   * This does the finalNormalization() (see above) for a whole
   * shell set of 3-center-2-el-integrals. Here the normalization
   * of the auxiliary shell is omitted.
   *
   * See normalizeShell(vector<double>*, unsigned int, unsigned int).
   *
   * @param integrals
   * @param amA, amB, amC
   */
  static void normalizeShellNoAux(Eigen::Ref<Eigen::VectorXd> integrals, unsigned int amA, unsigned int amB, unsigned int amC) {
    const long double undo = std::pow(1.0 * double_factorial(2 * amB - 1) * double_factorial(2 * amC - 1), 0.5);
    for (int aX = amA, index = 0; aX >= 0; --aX) {
      for (int aY = amA - aX; aY >= 0; --aY) {
        const long double normA = 1.0;
        for (int bX = amB; bX >= 0; --bX) {
          for (int bY = amB - bX; bY >= 0; --bY) {
            int bZ = amB - bY - bX;
            const long double normB = normA * finalNormalization(bX, bY, bZ);
            for (int cX = amC; cX >= 0; --cX) {
              for (int cY = amC - cX; cY >= 0; --cY) {
                int cZ = amC - cY - cX;
                integrals[index] *= (long double)(normB * finalNormalization(cX, cY, cZ) * undo);
                ++index;
              }
            }
          }
        }
      }
    }
  };
  /**
   * This does the finalNormalization() (see above) for a whole
   * shell set of 2-el-integrals, which are computed together anyway. The
   * ordering of the integrals is properly considered.
   *
   * See normalizeShell(vector<double>*, unsigned int, unsigned int).
   *
   * @param integrals
   * @param amA, amB, amC, amD
   */
  static void normalizeShell(Eigen::Ref<Eigen::VectorXd> integrals, unsigned int amA, unsigned int amB,
                             unsigned int amC, unsigned int amD) {
    const long double undo = std::pow(double_factorial(2 * amA - 1) * double_factorial(2 * amB - 1) *
                                          double_factorial(2 * amC - 1) * double_factorial(2 * amD - 1),
                                      0.5);
    for (int aX = amA, index = 0; aX >= 0; --aX) {
      for (int aY = amA - aX; aY >= 0; --aY) {
        int aZ = amA - aY - aX;
        const long double normA = finalNormalization(aX, aY, aZ);
        for (int bX = amB; bX >= 0; --bX) {
          for (int bY = amB - bX; bY >= 0; --bY) {
            int bZ = amB - bY - bX;
            const long double normB = normA * finalNormalization(bX, bY, bZ);
            for (int cX = amC; cX >= 0; --cX) {
              for (int cY = amC - cX; cY >= 0; --cY) {
                int cZ = amC - cY - cX;
                const long double normC = normB * finalNormalization(cX, cY, cZ);
                for (int dX = amD; dX >= 0; --dX) {
                  for (int dY = amD - dX; dY >= 0; --dY) {
                    int dZ = amD - dY - dX;
                    integrals[index] *= (long double)(normC * finalNormalization(dX, dY, dZ) * undo);
                    ++index;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};

} /* namespace Serenity */
#endif /* NORMALIZATION_H */
