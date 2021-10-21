/**
 * @file CCSD.cpp
 *
 * @date Mar 31, 2016
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
#include "postHF/CC/CCSD.h"
/* Include Serenity Internal Headers */
#include "basis/Basis.h"
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/transformer/Ao2MoTransformer.h"
#include "io/FormattedOutput.h"
#include "math/Matrix.h"
#include "math/RegularRankFourTensor.h"
#include "system/SystemController.h"

namespace Serenity {

CCSD::CCSD(std::shared_ptr<SystemController> systemController, double normThreshold, unsigned int maxCycles)
  : _eris(nullptr),
    _t1(nullptr),
    _t2(nullptr),
    _systemController(systemController),
    _converged(false),
    _MP2EnergyCorrection(0.0),
    _foo(nullptr),
    _fov(nullptr),
    _fvv(nullptr),
    _normThreshold(normThreshold),
    _maxCycles(maxCycles) {
  assert(_systemController);
}

std::pair<double, double> CCSD::calculateElectronicEnergyCorrections() {
  if (_converged)
    return std::make_pair(_MP2EnergyCorrection, calculateCCSDEnergyCorrection());
  prepERIS();
  initializeAmplitudes();
  double oldEnergy = calculateCCSDEnergyCorrection();
  _MP2EnergyCorrection = oldEnergy;
  unsigned int cycle = 0;
  std::cout << "Maximum number of cycles: " << _maxCycles << std::endl;
  printSmallCaption("Amplitude Optimization");
  printf("%4s %2s %10s      %12s  %10s\n", "", "#", "E_ccsd", "Delta E_ccsd", "Amp. Norm");
  while (!_converged) {
    if (cycle == _maxCycles)
      break;
    updateF();
    double norm = updateAmplitudes();
    double energy = calculateCCSDEnergyCorrection();
    printf("%4s %2d %-13.10f  %-13.10f  %-13.10f\n", "", cycle + 1, energy, energy - oldEnergy, norm);
    if (fabs(norm) < _normThreshold) {
      std::cout << std::endl << "    CCSD equations converged, exiting." << std::endl;
      _converged = true;
      return std::make_pair(_MP2EnergyCorrection, energy);
    }
    ++cycle;
    oldEnergy = energy;
  }
  // TODO handle failed convergence;
  throw SerenityError("Failed to converge CCSD amplitudes");
}

double CCSD::calculateCCSDEnergyCorrection() {
  auto& eris = *_eris;
  auto& t1 = *_t1;
  auto& t2 = *_t2;

  auto nocc = t1.rows();
  auto nvirt = t1.cols();
  double CCSDEnergyCorrection = 0.0;
  for (unsigned int i = 0; i != nocc; ++i) {
    for (unsigned int a = 0; a != nvirt; ++a) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int b = 0; b != nvirt; ++b) {
          double theta = 2.0 * (t2(i, a)(j, b) + t1(i, a) * t1(j, b));
          theta -= (t2(j, a)(i, b) + t1(j, a) * t1(i, b));
          CCSDEnergyCorrection += theta * eris(i, a + nocc, j, b + nocc);
        }
      }
    }
  }
  return CCSDEnergyCorrection;
}

double CCSD::updateAmplitudes() {
  auto& eris = *_eris;
  auto& t1 = *_t1;
  auto& t2 = *_t2;
  auto nocc = t1.rows();
  auto nvirt = t1.cols();
  auto& foo = *_foo;
  auto& fov = *_fov;
  auto& fvv = *_fvv;

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

  /*
   * t1 update
   */

  auto newt1_ptr = std::unique_ptr<Matrix<double>>(new Matrix<double>(nocc, nvirt));
  auto& newt1 = *newt1_ptr;
  newt1.setZero();
  Matrix<Matrix<double>> tmp1(nocc, nvirt, Eigen::MatrixXd::Zero(nocc, nvirt));

#pragma omp parallel shared(newt1, tmp1)
  {
#pragma omp for schedule(dynamic)
    // make sure tmp1 is 0
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        tmp1(i, a).setZero();
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            tmp1(i, a)(j, b) = t2(i, a)(j, b) * 2.0 - t2(j, a)(i, b);
          }
        }
      }
    }

#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          newt1(i, a) -= t1(j, a) * foo(j, i);
          for (unsigned int b = 0; b != nvirt; ++b) {
            newt1(i, a) += t1(j, b) * (2.0 * eris(i, a + nocc, j, b + nocc) - eris(j, i, a + nocc, b + nocc));
            newt1(i, a) += fov(j, b) * tmp1(i, a)(j, b);
            for (unsigned int c = 0; c != nvirt; ++c) {
              newt1(i, a) += tmp1(i, b)(j, c) * eris(j, c + nocc, b + nocc, a + nocc);
            }
            for (unsigned int k = 0; k != nocc; ++k) {
              newt1(i, a) -= tmp1(j, b)(k, a) * eris(j, b + nocc, k, i);
            }
          }
        }
        for (unsigned int b = 0; b != nvirt; ++b) {
          newt1(i, a) += t1(i, b) * fvv(a, b);
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        newt1(i, a) /= eia(i, a);
      }
    }
  } /* end of parallel region */

  /*
   * t2 update
   */
  auto newt2_ptr =
      std::unique_ptr<Matrix<Matrix<double>>>(new Matrix<Matrix<double>>(nocc, nvirt, Eigen::MatrixXd::Zero(nocc, nvirt)));
  auto& newt2 = *newt2_ptr;
  Matrix<double> ft_ij(foo);
  Matrix<double> ft_ab(fvv);

  Matrix<Matrix<double>> tau(nocc, nvirt, Eigen::MatrixXd::Zero(nocc, nvirt));
  { /* scope of tw */
    Matrix<Matrix<double>> tw(nocc, nvirt, Eigen::MatrixXd::Zero(nocc, nvirt));
#pragma omp parallel shared(newt2, tmp1, ft_ij, ft_ab, tau, tw)
    {
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              newt2(i, a)(j, b) = eris(i, a + nocc, j, b + nocc);
            }
          }
        }
      }

#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            ft_ij(i, j) += 0.5 * t1(j, a) * fov(i, a);
          }
        }
      }

#pragma omp for schedule(dynamic)
      for (unsigned int a = 0; a < nvirt; ++a) {
        for (unsigned int b = 0; b != nvirt; ++b) {
          for (unsigned int i = 0; i != nocc; ++i) {
            ft_ab(a, b) -= 0.5 * t1(i, a) * fov(i, b);
          }
        }
      }

      // make sure tw is 0
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          tw(i, a).setZero();
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int k = 0; k != nocc; ++k) {
                tw(i, a)(j, b) -= ft_ij(k, j) * t2(k, b)(i, a);
              }
              for (unsigned int c = 0; c != nvirt; ++c) {
                tw(i, a)(j, b) += ft_ab(b, c) * t2(i, a)(j, c);
              }
            }
          }
        }
      }

      //    ft_ij = Matrix<double>(0,0);
      //    ft_ab = Matrix<double>(0,0);
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              double tmp = 0.0;
              for (unsigned int c = 0; c != nvirt; ++c) {
                tmp += t1(i, c) * eris(k, c + nocc, j, b + nocc);
              }
              for (unsigned int a = 0; a != nvirt; ++a) {
                tw(i, a)(j, b) -= t1(k, a) * tmp;
              }
            }
          }
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int a = 0; a != nvirt; ++a) {
              double tmp = 0.0;
              for (unsigned int c = 0; c != nvirt; ++c) {
                tmp += t1(i, c) * eris(k, j, a + nocc, c + nocc);
              }
              for (unsigned int b = 0; b != nvirt; ++b) {
                tw(i, a)(j, b) -= t1(k, b) * tmp;
              }
            }
          }
        }
      }

      // make sure tau is 0
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          tau(i, a).setZero();
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              tau(i, a)(j, b) = 2.0 * t2(i, a)(j, b) - t2(j, a)(i, b) - t1(i, b) * t1(j, a) * 2.0;
            }
          }
        }
      }

      // reset tmp1
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          tmp1(i, a).setZero();
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              tmp1(i, a)(j, b) += eris(i, a + nocc, j, b + nocc);
              for (unsigned int k = 0; k != nocc; ++k) {
                for (unsigned int c = 0; c != nvirt; ++c) {
                  tmp1(i, b)(j, a) += 0.5 * tau(j, a)(k, c) * eris(i, b + nocc, k, c + nocc);
                  tmp1(i, b)(j, a) -= 0.5 * t2(j, a)(k, c) * eris(i, c + nocc, k, b + nocc);
                }
              }
              for (unsigned int k = 0; k != nocc; ++k) {
                tmp1(i, b)(j, a) -= t1(k, a) * eris(i, b + nocc, k, j);
              }
              for (unsigned int c = 0; c != nvirt; ++c) {
                tmp1(i, b)(j, a) += t1(j, c) * eris(i, b + nocc, a + nocc, c + nocc);
              }
            }
          }
        }
      }

#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int k = 0; k != nocc; ++k) {
            for (unsigned int c = 0; c != nvirt; ++c) {
              double theta = t2(i, a)(k, c) * 2.0 - t2(k, a)(i, c);
              for (unsigned int j = 0; j != nocc; ++j) {
                for (unsigned int b = 0; b != nvirt; ++b) {
                  tw(i, a)(j, b) += theta * tmp1(k, c)(j, b);
                }
              }
            }
          }
        }
      }

#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              tau(i, a)(j, b) = 0.5 * t2(i, a)(j, b) + t1(i, a) * t1(j, b);
            }
          }
        }
      }

      // reset tmp1
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          tmp1(i, a).setZero();
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int a = 0; a != nvirt; ++a) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              tmp1(i, a)(j, b) -= 1.0 * eris(j, i, a + nocc, b + nocc);
              for (unsigned int k = 0; k != nocc; ++k) {
                for (unsigned int c = 0; c != nvirt; ++c) {
                  tmp1(i, b)(j, a) += tau(j, c)(k, a) * eris(i, c + nocc, k, b + nocc);
                }
              }
              for (unsigned int k = 0; k != nocc; ++k) {
                tmp1(i, b)(j, a) += t1(k, a) * eris(k, b + nocc, i, j);
              }
              for (unsigned int c = 0; c != nvirt; ++c) {
                tmp1(i, b)(j, a) -= t1(j, c) * eris(i, c + nocc, a + nocc, b + nocc);
              }
            }
          }
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int k = 0; k != nocc; ++k) {
                for (unsigned int c = 0; c != nvirt; ++c) {
                  tw(i, a)(j, b) += t2(i, a)(k, c) * tmp1(k, c)(j, b);
                  tw(i, b)(j, a) += t2(i, c)(k, a) * tmp1(k, c)(j, b);
                }
              }
            }
          }
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                tw(j, b)(i, a) += t1(j, c) * eris(i, a + nocc, c + nocc, b + nocc);
              }
              for (unsigned int k = 0; k != nocc; ++k) {
                tw(j, b)(i, a) -= t1(k, a) * eris(j, b + nocc, i, k);
              }
            }
          }
        }
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              newt2(i, a)(j, b) += tw(i, a)(j, b) + tw(j, b)(i, a);
            }
          }
        }
      }

#pragma omp for schedule(dynamic)
      for (unsigned int i = 0; i < nocc; ++i) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int j = 0; j != nocc; ++j) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              tau(i, a)(j, b) = t2(i, a)(j, b) + t1(i, a) * t1(j, b);
            }
          }
        }
      }
    } /* end of parallel region */
  }   /* scope of tw */

  // new tw
  Matrix<Matrix<double>> tw(nocc, nocc, Eigen::MatrixXd::Zero(nocc, nocc));
#pragma omp parallel shared(newt2, tmp1, ft_ij, ft_ab, tau, tw)
  {
    // make sure tw is 0
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j < nocc; ++j) {
        tw(i, j).setZero();
      }
    }

#pragma omp for schedule(dynamic)
    for (unsigned int j = 0; j < nocc; ++j) {
      for (unsigned int l = 0; l != nocc; ++l) {
        for (unsigned int i = 0; i != nocc; ++i) {
          for (unsigned int k = 0; k != nocc; ++k) {
            double tmp = 0.0;
            for (unsigned int a = 0; a != nvirt; ++a) {
              tmp += (t1(l, a) * eris(j, a + nocc, i, k));
            }
#pragma omp critical
            { tw(j, l)(i, k) += eris(j, l, i, k) + tmp; }
            for (unsigned int a = 0; a != nvirt; ++a) {
              for (unsigned int b = 0; b != nvirt; ++b) {
                tmp += eris(i, a + nocc, j, b + nocc) * tau(k, a)(l, b);
              }
            }
#pragma omp critical
            { tw(i, k)(j, l) += tmp; }
          }
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            for (unsigned int k = 0; k != nocc; ++k) {
              for (unsigned int l = 0; l != nocc; ++l) {
                newt2(i, a)(j, b) += tw(k, i)(l, j) * tau(k, a)(l, b);
              }
            }
          }
        }
      }
    }

    // reset tmp1
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        tmp1(i, a).setZero();
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            double tmp = 0.0;
            for (unsigned int d = 0; d != nvirt; ++d) {
              for (unsigned int c = 0; c != nvirt; ++c) {
                tmp += tau(i, c)(j, d) * eris(k, c + nocc, b + nocc, d + nocc);
              }
            }
            for (unsigned int a = 0; a != nvirt; ++a) {
              tmp1(i, a)(j, b) += tmp * t1(k, a);
            }
          }
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            newt2(i, a)(j, b) -= (tmp1(i, a)(j, b) + tmp1(j, b)(i, a));
          }
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            for (unsigned int c = 0; c != nvirt; ++c) {
              for (unsigned int d = 0; d != nvirt; ++d) {
                newt2(i, a)(j, b) += (t2(i, c)(j, d) + t1(i, c) * t1(j, d)) * eris(a + nocc, c + nocc, b + nocc, d + nocc);
              }
            }
          }
        }
      }
    }
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            newt2(i, a)(j, b) /= (eia(i, a) + eia(j, b));
          }
        }
      }
    }
  } /* end of parallel region */

  // measure amplitude norm
  double norm = 0.0;
#pragma omp parallel shared(norm)
  {
#pragma omp parallel for schedule(static) reduction(+ : norm)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        norm += (newt2(i, a) - t2(i, a)).norm();
      }
    }
  } /* end of parallel region */
  norm += (newt1 - t1).norm();

  _t1 = std::move(newt1_ptr);
  _t2 = std::move(newt2_ptr);

  return norm;
}

void CCSD::updateF() {
  auto& eris = *_eris;
  auto& t1 = *_t1;
  auto& t2 = *_t2;
  auto nocc = t1.rows();
  auto nvirt = t1.cols();

  if (!_foo)
    _foo = std::unique_ptr<Matrix<double>>(new Matrix<double>(nocc, nocc));
  if (!_fov)
    _fov = std::unique_ptr<Matrix<double>>(new Matrix<double>(nocc, nvirt));
  if (!_fvv)
    _fvv = std::unique_ptr<Matrix<double>>(new Matrix<double>(nvirt, nvirt));
  auto& foo = *_foo;
  auto& fov = *_fov;
  auto& fvv = *_fvv;
  foo.setZero();
  fov.setZero();
  fvv.setZero();

  // Precalc theta
  Matrix<Matrix<double>> theta(nocc, nvirt, Eigen::MatrixXd::Zero(nocc, nvirt));
#pragma omp parallel shared(theta, foo, fov, fvv, t1, t2, eris)
  {
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            theta(i, a)(j, b) = 2.0 * (t2(i, a)(j, b) + 0.5 * t1(i, a) * t1(j, b));
            theta(i, a)(j, b) -= (t2(j, a)(i, b) + 0.5 * t1(j, a) * t1(i, b));
          }
        }
      }
    }

#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int j = 0; j != nocc; ++j) {
        for (unsigned int a = 0; a != nvirt; ++a) {
          for (unsigned int k = 0; k != nocc; ++k) {
            for (unsigned int b = 0; b != nvirt; ++b) {
              foo(i, j) += eris(i, a + nocc, k, b + nocc) * theta(j, a)(k, b);
            }
            foo(i, j) += t1(k, a) * (2.0 * eris(k, a + nocc, i, j) - eris(i, a + nocc, k, j));
          }
        }
      }
    }

#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int k = 0; k != nocc; ++k) {
          for (unsigned int c = 0; c != nvirt; ++c) {
            fov(i, a) += t1(k, c) * (2.0 * eris(i, a + nocc, k, c + nocc) - eris(k, a + nocc, i, c + nocc));
          }
        }
      }
    }

#pragma omp for schedule(dynamic)
    for (unsigned int a = 0; a < nvirt; ++a) {
      for (unsigned int b = 0; b != nvirt; ++b) {
        for (unsigned int k = 0; k < nocc; ++k) {
          for (unsigned int c = 0; c != nvirt; ++c) {
            fvv(a, b) += t1(k, c) * (2.0 * eris(k, c + nocc, a + nocc, b + nocc) - eris(k, b + nocc, a + nocc, c + nocc));
            for (unsigned int j = 0; j != nocc; ++j) {
              fvv(a, b) -= theta(k, c)(j, a) * eris(k, c + nocc, j, b + nocc);
            }
          }
        }
      }
    }
  } /*end of parallel region*/
}

void CCSD::initializeAmplitudes() {
  const unsigned int nBasisFunc = _systemController->getBasisController()->getNBasisFunctions();
  const auto& orbitals = _systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>();
  const auto& orbitalEnergies = orbitals->getEigenvalues();
  const auto& nElectrons = _systemController->getNElectrons<Options::SCF_MODES::RESTRICTED>();
  const auto nocc = nElectrons / 2;
  const auto nvirt = nBasisFunc - nElectrons / 2;
  auto& eris = *_eris;

  _t1 = std::unique_ptr<Matrix<double>>(new Matrix<double>(nocc, nvirt));
  auto& t1 = *_t1;
  t1.setZero();

  _t2 = std::unique_ptr<Matrix<Matrix<double>>>(new Matrix<Matrix<double>>(nocc, nvirt, Eigen::MatrixXd::Zero(nocc, nvirt)));
  auto& t2 = *_t2;

  // Orbital energy differences
  Matrix<double> eia(nocc, nvirt);
#pragma omp parallel shared(eia, t2)
  {
#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = nocc; a != nBasisFunc; ++a) {
        eia(i, a - nocc) = orbitalEnergies[i] - orbitalEnergies[a];
      }
    }

#pragma omp for schedule(dynamic)
    for (unsigned int i = 0; i < nocc; ++i) {
      for (unsigned int a = 0; a != nvirt; ++a) {
        for (unsigned int j = 0; j != nocc; ++j) {
          for (unsigned int b = 0; b != nvirt; ++b) {
            t2(i, a)(j, b) = eris(i, a + nocc, j, b + nocc) / (eia(i, a) + eia(j, b));
          }
        }
      }
    }
  } /*end of parallel region*/
}

void CCSD::prepERIS() {
  const auto& basisController = _systemController->getBasisController();
  const unsigned int nBasisFunc = basisController->getNBasisFunctions();
  if (!_eris) {
    _eris = std::unique_ptr<RegularRankFourTensor<double>>(new RegularRankFourTensor<double>(nBasisFunc, 0.0));
  }
  auto& eris = *_eris;

  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, 1E-10);

  auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                 const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId;
    eris(b, a, i, j) = integral(0);
    eris(b, a, j, i) = integral(0);
    eris(a, b, j, i) = integral(0);
    eris(a, b, i, j) = integral(0);
    eris(i, j, b, a) = integral(0);
    eris(i, j, a, b) = integral(0);
    eris(j, i, b, a) = integral(0);
    eris(j, i, a, b) = integral(0);
  };

  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coefficients =
      _systemController->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  looper.loop(storeERIS, coefficients.lpNorm<Eigen::Infinity>());

  Ao2MoTransformer ao2mo(basisController);

  ao2mo.transformTwoElectronIntegrals(eris, eris, coefficients, nBasisFunc);
}

} /* namespace Serenity */
