/**
 * @file CCSD_T.cpp
 *
 * @date Apr 14, 2016
 * @author Jan Unsleber
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
#include "postHF/CC/CCSD_T.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "math/Matrix.h"
#include "system/SystemController.h"

namespace Serenity {

double CCSD_T::calculateTripplesCorrection() {
  double correction = 0.0;
  assert(this->_converged);

  auto& t1 = *_t1;
  auto& t2 = *_t2;
  auto& eris = *_eris;
  const unsigned int nocc = t1.rows();
  const unsigned int nvirt = t1.cols();

  Matrix<Matrix<Matrix<double>>> r(nocc, nocc, Matrix<Matrix<double>>(nocc, nvirt, Eigen::MatrixXd::Zero(nvirt, nvirt)));
  Matrix<Matrix<Matrix<double>>> w(nocc, nocc, Matrix<Matrix<double>>(nocc, nvirt, Eigen::MatrixXd::Zero(nvirt, nvirt)));

#pragma omp parallel shared(r, w, t1, t2, eris)
  {
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                r(i, j)(k, a)(b, c) = 0.0;
                for (unsigned int l = 0; l != nocc; ++l) {
                  r(i, j)(k, a)(b, c) -= eris(i, a + nocc, j, l) * t2(l, b)(k, c);
                }
                for (unsigned int d = 0; d != nvirt; ++d) {
                  r(i, j)(k, a)(b, c) += eris(i, a + nocc, b + nocc, d + nocc) * t2(k, c)(j, d);
                }
              }
            }
          }
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                w(i, j)(k, a)(b, c) = r(i, j)(k, a)(b, c);
                w(i, j)(k, a)(b, c) += 0.5 * eris(i, a + nocc, j, b + nocc) * t1(k, c);
              }
            }
          }
        }
      }
    }
  } /* End of parallel region */

  p6(r);
  r6(r);
  p6(w);

  const unsigned int nBasisFunc = _systemController->getBasisController()->getNBasisFunctions();
  const auto& orbitals = _systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>();
  const auto& orbitalEnergies = orbitals->getEigenvalues();

  // Orbital energy differences
  Matrix<double> eia(nocc, nvirt);
  for (unsigned int i = 0; i != nocc; ++i) {
    for (unsigned int a = nocc; a != nBasisFunc; ++a) {
      eia(i, a - nocc) = orbitalEnergies[i] - orbitalEnergies[a];
    }
  }

#pragma omp parallel shared(w, r, correction)
  {
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                w(i, j)(k, a)(b, c) /= (eia(i, a) + eia(j, b) + eia(k, c));
              }
            }
          }
        }
      }
    }
#pragma omp for schedule(static) reduction(+ : correction)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                correction += 1.0 / 3.0 * r(i, j)(k, a)(b, c) * w(i, j)(k, a)(b, c);
              }
            }
          }
        }
      }
    }
  } /* End of parallel region */

  return correction;
}

void CCSD_T::p6(Matrix<Matrix<Matrix<double>>>& mat) {
  const unsigned int nocc = _t1->rows();
  const unsigned int nvirt = _t1->cols();
  Matrix<Matrix<Matrix<double>>> tmp(nocc, nocc, Matrix<Matrix<double>>(nocc, nvirt, Eigen::MatrixXd::Zero(nvirt, nvirt)));
#pragma omp parallel shared(tmp, mat)
  {
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                tmp(i, j)(k, a)(b, c) = mat(i, j)(k, a)(b, c) + mat(i, k)(j, a)(c, b);
              }
            }
          }
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                mat(i, j)(k, a)(b, c) = tmp(i, j)(k, a)(b, c) + tmp(j, i)(k, b)(a, c) + tmp(k, i)(j, c)(a, b);
              }
            }
          }
        }
      }
    }
  } /* End of parallel region */
}

void CCSD_T::r6(Matrix<Matrix<Matrix<double>>>& mat) {
  const unsigned int nocc = _t1->rows();
  const unsigned int nvirt = _t1->cols();
  Matrix<Matrix<Matrix<double>>> tmp(nocc, nocc, Matrix<Matrix<double>>(nocc, nvirt, Eigen::MatrixXd::Zero(nvirt, nvirt)));
#pragma omp parallel shared(tmp, mat)
  {
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                tmp(i, j)(k, a)(b, c) = mat(i, j)(k, a)(b, c);
              }
            }
          }
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                mat(i, j)(k, a)(b, c) = 4.0 * tmp(i, j)(k, a)(b, c) + tmp(i, j)(k, c)(a, b) + tmp(i, j)(k, b)(c, a) -
                                        2.0 * tmp(i, j)(k, c)(b, a) - 2.0 * tmp(i, j)(k, a)(c, b) -
                                        2.0 * tmp(i, j)(k, b)(a, c);
              }
            }
          }
        }
      }
    }
  } /* End of parallel region */
}

} /* namespace Serenity */
