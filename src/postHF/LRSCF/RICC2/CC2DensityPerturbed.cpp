/**
 * @file CC2DensityPerturbed.cpp
 *
 * @date Apr 20, 2020
 * @author Niklas Niemeyer
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

/* Include Class Header*/
#include "postHF/LRSCF/RICC2/CC2Controller.h"

/* Include Serenity Internal Headers */
#include "integrals/looper/TwoElecThreeCenterCalculator.h"
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<>
void CC2Controller<RESTRICTED>::calcPerturbedDensities(std::vector<Eigen::MatrixXd>& solutionvectors,
                                                       std::vector<double> frequencies,
                                                       std::vector<Eigen::MatrixXd>& densityMatrices,
                                                       Eigen::MatrixXd dipoleIntegrals) {
  Timings::takeTime("CC2 -      Nonit. Pert. Dens.");

  unsigned mo = _no + _nv;
  unsigned nFreq = frequencies.size();
  unsigned nCart = solutionvectors[0].cols();
  densityMatrices.resize(nFreq);
  printf("  Calculating %3i eta density matrices.\n\n", nFreq * nCart);

  // CC2 density matrix intermediate.
  Eigen::MatrixXd CC2virt = 0.5 * _gsDensity.bottomRightCorner(_nv, _nv);
  Eigen::MatrixXd CC2occ = 0.5 * _gsDensity.topLeftCorner(_no, _no);
  CC2occ.diagonal() -= Eigen::VectorXd::Constant(_no, 1.0);

  for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
    densityMatrices[iFreq] = Eigen::MatrixXd::Zero(mo * mo, nCart);
  }

  Eigen::Map<Eigen::MatrixXd> gsl(_gsLagrange.data(), _nv, _no);
  this->performTransformation(_Zia, _gsLagrange, true);

  printf("  Non-iterative part.\n\n");
  printTableHead("  root    time(min)");

  // Non-iterative contribution.
  for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
    auto itStart = std::chrono::steady_clock::now();
    _Xia = (*_Bia);
    for (unsigned i = 0, ia = 0; i < _no; ++i) {
      for (unsigned a = 0; a < _nv; ++a, ++ia) {
        _Xia.row(ia) *= std::sqrt(_weights(iRoot)) * std::exp(-(_e(_no + a) - _e(i)) * _roots(iRoot));
      }
    }

    for (size_t iCart = 0; iCart < nCart; ++iCart) {
      // Dipole Integrals.
      Eigen::Map<Eigen::MatrixXd> X(dipoleIntegrals.col(iCart).data(), mo, mo);
      Eigen::Ref<Eigen::MatrixXd> Xoo = X.topLeftCorner(_no, _no);
      Eigen::Ref<Eigen::MatrixXd> Xvv = X.bottomRightCorner(_nv, _nv);

      _Yia = Eigen::MatrixXd::Zero(_nv * _no, _nx);
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> xia(_Xia.col(Q).data(), _nv, _no);
        Eigen::Map<Eigen::MatrixXd> yia(_Yia.col(Q).data(), _nv, _no);
        yia = Xvv * xia - xia * Xoo;
      }

      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = i; j < _no; ++j) {
          Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];

            Eigen::Map<Eigen::MatrixXd> D(densityMatrices[iFreq].col(iCart).data(), mo, mo);

            Eigen::Ref<Eigen::MatrixXd> Dov = D.topRightCorner(_no, _nv);
            Eigen::Ref<Eigen::MatrixXd> Dvv = D.bottomRightCorner(_nv, _nv);

            Eigen::Map<Eigen::MatrixXd> rev(solutionvectors[iFreq].col(iCart).data(), _nv, _no);

            Eigen::MatrixXd xij = this->getLTDoubles(i, j, frequency);
            Eigen::MatrixXd xxij = _soss * xij - _sss * xij.transpose();

            Dvv.noalias() += tij * xij.transpose();
            Dov.row(i).noalias() += xxij * gsl.col(j);
            if (i != j) {
              Dvv.noalias() += tij.transpose() * xij;
              Dov.row(j).noalias() += xxij.transpose() * gsl.col(i);
            }
          }
        }
      }

      this->reorderTensor(_Fia, _nv, _no);
      this->reorderTensor(_gsLagrange, _nv, _no);
      this->reorderTensor(_Xia, _nv, _no);
      this->reorderTensor(_Yia, _nv, _no);
      this->reorderTensor(_Zia, _nv, _no);
      this->reorderTensor(*_Jia, _nv, _no);
      if (_Bia != _Jia) {
        this->reorderTensor(*_Bia, _nv, _no);
      }

      for (size_t a = 0; a < _nv; ++a) {
        for (size_t b = a; b < _nv; ++b) {
          Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];
            Eigen::Map<Eigen::MatrixXd> D(densityMatrices[iFreq].col(iCart).data(), mo, mo);
            Eigen::Ref<Eigen::MatrixXd> Doo = D.topLeftCorner(_no, _no);

            Eigen::MatrixXd xab = this->getLTDoublesV(a, b, frequency);
            Doo.noalias() -= xab * tab.transpose();
            if (a != b) {
              Doo.noalias() -= xab.transpose() * tab;
            }
          }
        }
      }

      this->reorderTensor(_Fia, _no, _nv);
      this->reorderTensor(_gsLagrange, _no, _nv);
      this->reorderTensor(_Xia, _no, _nv);
      this->reorderTensor(_Zia, _no, _nv);
      this->reorderTensor(*_Jia, _no, _nv);
      if (_Bia != _Jia) {
        this->reorderTensor(*_Bia, _no, _nv);
      }
    }

    auto itEnd = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();
    printf(" %5i %12.3f\n", iRoot + 1, duration / 60.0);
  }
  Timings::timeTaken("CC2 -      Nonit. Pert. Dens.");

  Timings::takeTime("CC2 -         It. Pert. Dens.");
  printf("\n  Iterative part.\n\n");
  printf("  root  time(min)\n");
  printf(" ------------------");
  // Iterative contribution.
  for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
    double frequency = frequencies[iFreq];
    for (unsigned iCart = 0; iCart < nCart; ++iCart) {
      auto itStart = std::chrono::steady_clock::now();

      Eigen::Map<Eigen::MatrixXd> D(densityMatrices[iFreq].col(iCart).data(), mo, mo);
      Eigen::Ref<Eigen::VectorXd> solutionvector = solutionvectors[iFreq].col(iCart);

      Eigen::Ref<Eigen::MatrixXd> Doo = D.topLeftCorner(_no, _no);
      Eigen::Ref<Eigen::MatrixXd> Dov = D.topRightCorner(_no, _nv);
      Eigen::Ref<Eigen::MatrixXd> Dvv = D.bottomRightCorner(_nv, _nv);

      Eigen::Map<Eigen::MatrixXd> rev(solutionvector.data(), _nv, _no);

      // Start with easy contributions.
      Dov.noalias() += rev.transpose();
      Doo.noalias() -= rev.transpose() * gsl;
      Dvv.noalias() += gsl * rev.transpose();
      Dov.noalias() -= rev.transpose() * CC2virt;
      Dov.noalias() += CC2occ * rev.transpose();

      this->performTransformation(_Xia, solutionvector, false);
      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = i; j < _no; ++j) {
          auto tij = this->getGLagrangeAmplitudes(i, j);
          auto rij = this->getRightAmplitudesA(i, j, frequency);
          auto rrij = _soss * rij - _sss * rij.transpose();

          Dvv.noalias() += tij * rij.transpose();
          Dov.row(i).noalias() += rrij * gsl.col(j);
          if (i != j) {
            Dvv.noalias() += tij.transpose() * rij;
            Dov.row(j).noalias() += rrij.transpose() * gsl.col(i);
          }
        }
      }

      this->reorderTensor(_Fia, _nv, _no);
      this->reorderTensor(_gsLagrange, _nv, _no);
      this->reorderTensor(*_Jia, _nv, _no);
      if (_Bia != _Jia) {
        this->reorderTensor(*_Bia, _nv, _no);
      }
      this->reorderTensor(_Xia, _nv, _no);
      this->reorderTensor(_Zia, _nv, _no);

      for (size_t a = 0; a < _nv; ++a) {
        for (size_t b = a; b < _nv; ++b) {
          auto tab = this->getGLagrangeAmplitudesV(a, b);
          auto rab = this->getRightAmplitudesAV(a, b, frequency);

          Doo.noalias() -= rab * tab.transpose();
          if (a != b) {
            Doo.noalias() -= rab.transpose() * tab;
          }
        }
      }

      // done with ab-wise amps, must order ints back
      this->reorderTensor(_Fia, _no, _nv);
      this->reorderTensor(_gsLagrange, _no, _nv);
      this->reorderTensor(_Zia, _no, _nv);
      this->reorderTensor(*_Jia, _no, _nv);
      if (_Bia != _Jia) {
        this->reorderTensor(*_Bia, _no, _nv);
      }

      auto itEnd = std::chrono::steady_clock::now();
      double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();
      auto str = ((iFreq * nCart + iCart) % 5 == 0) ? "\n %5i %9.3f" : " %5i %9.3f";
      printf(str, iFreq * nCart + iCart + 1, duration / 60.0);
    } /* iCart */
  }   /* iFreq */
  printf("\n\n");
  Timings::timeTaken("CC2 -         It. Pert. Dens.");
} /* this->calcPerturbedDensities() restricted */

template<>
void CC2Controller<UNRESTRICTED>::calcPerturbedDensities(std::vector<Eigen::MatrixXd>& solutionvectors,
                                                         std::vector<double> frequencies,
                                                         std::vector<Eigen::MatrixXd>& densityMatrices,
                                                         Eigen::MatrixXd dipoleIntegrals) {
  Timings::takeTime("CC2 -      Nonit. Pert. Dens.");

  unsigned mo_a = _no_a + _nv_a;
  unsigned mo_b = _no_b + _nv_b;
  unsigned nFreq = frequencies.size();
  unsigned nCart = solutionvectors[0].cols();
  densityMatrices.resize(nFreq);
  printf("  Calculating %3i eta density matrices.\n\n", nFreq * nCart);

  // CC2 density matrix intermediate.
  Eigen::MatrixXd CC2virt_a = _gsDensity.alpha.bottomRightCorner(_nv_a, _nv_a);
  Eigen::MatrixXd CC2virt_b = _gsDensity.beta.bottomRightCorner(_nv_b, _nv_b);
  Eigen::MatrixXd CC2occ_a = _gsDensity.alpha.topLeftCorner(_no_a, _no_a);
  Eigen::MatrixXd CC2occ_b = _gsDensity.beta.topLeftCorner(_no_b, _no_b);
  CC2occ_a.diagonal() -= Eigen::VectorXd::Constant(_no_a, 1.0);
  CC2occ_b.diagonal() -= Eigen::VectorXd::Constant(_no_b, 1.0);

  for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
    densityMatrices[iFreq] = Eigen::MatrixXd::Zero(mo_a * mo_a + mo_b * mo_b, nCart);
  }

  Eigen::Map<Eigen::MatrixXd> gsl_a(_gsLagrange.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> gsl_b(_gsLagrange.data() + _alpha, _nv_b, _no_b);
  this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
  this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);

  printf("  Non-iterative part.\n\n");
  printTableHead("  root    time(min)");

  // Non-iterative contribution.
  for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
    auto itStart = std::chrono::steady_clock::now();
    _Xia.alpha = _Bia->alpha;
    _Xia.beta = _Bia->beta;
    for (unsigned i = 0, ia = 0; i < _no_a; ++i) {
      for (unsigned a = 0; a < _nv_a; ++a, ++ia) {
        _Xia.alpha.row(ia) *= std::sqrt(_weights(iRoot)) * std::exp(-(_e.alpha(_no_a + a) - _e.alpha(i)) * _roots(iRoot));
      }
    }
    for (unsigned j = 0, jb = 0; j < _no_b; ++j) {
      for (unsigned b = 0; b < _nv_b; ++b, ++jb) {
        _Xia.beta.row(jb) *= std::sqrt(_weights(iRoot)) * std::exp(-(_e.beta(_no_b + b) - _e.beta(j)) * _roots(iRoot));
      }
    }

    for (size_t iCart = 0; iCart < nCart; ++iCart) {
      // Dipole Integrals.
      Eigen::Map<Eigen::MatrixXd> X_a(dipoleIntegrals.col(iCart).data(), mo_a, mo_a);
      Eigen::Map<Eigen::MatrixXd> X_b(dipoleIntegrals.col(iCart).data() + mo_a * mo_a, mo_b, mo_b);
      Eigen::Ref<Eigen::MatrixXd> Xoo_a = X_a.topLeftCorner(_no_a, _no_a);
      Eigen::Ref<Eigen::MatrixXd> Xoo_b = X_b.topLeftCorner(_no_b, _no_b);
      Eigen::Ref<Eigen::MatrixXd> Xvv_a = X_a.bottomRightCorner(_nv_a, _nv_a);
      Eigen::Ref<Eigen::MatrixXd> Xvv_b = X_b.bottomRightCorner(_nv_b, _nv_b);

      _Yia.alpha = Eigen::MatrixXd::Zero(_nv_a * _no_a, _nx);
      _Yia.beta = Eigen::MatrixXd::Zero(_nv_b * _no_b, _nx);
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> xia_a(_Xia.alpha.col(Q).data(), _nv_a, _no_a);
        Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.col(Q).data(), _nv_a, _no_a);
        yia_a = Xvv_a * xia_a - xia_a * Xoo_a;

        Eigen::Map<Eigen::MatrixXd> xia_b(_Xia.beta.col(Q).data(), _nv_b, _no_b);
        Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.col(Q).data(), _nv_b, _no_b);
        yia_b = Xvv_b * xia_b - xia_b * Xoo_b;
      }

      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];

            Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[iFreq].col(iCart).data(), mo_a, mo_a);

            Eigen::Ref<Eigen::MatrixXd> Dov_a = D_a.topRightCorner(_no_a, _nv_a);
            Eigen::Ref<Eigen::MatrixXd> Dvv_a = D_a.bottomRightCorner(_nv_a, _nv_a);

            Eigen::Map<Eigen::MatrixXd> rev_a(solutionvectors[iFreq].col(iCart).data(), _nv_a, _no_a);

            Eigen::MatrixXd xij = this->getLTDoubles(i, j, frequency, 1);
            Eigen::MatrixXd xxij = _sss * xij - _sss * xij.transpose();

            Dvv_a.noalias() += tij * xij.transpose();
            Dov_a.row(i).noalias() += xxij * gsl_a.col(j);
            if (i != j) {
              Dvv_a.noalias() += tij.transpose() * xij;
              Dov_a.row(j).noalias() += xxij.transpose() * gsl_a.col(i);
            }
          }
        }
      }

      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];

            Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[iFreq].col(iCart).data() + mo_a * mo_a, mo_b, mo_b);

            Eigen::Ref<Eigen::MatrixXd> Dov_b = D_b.topRightCorner(_no_b, _nv_b);
            Eigen::Ref<Eigen::MatrixXd> Dvv_b = D_b.bottomRightCorner(_nv_b, _nv_b);

            Eigen::Map<Eigen::MatrixXd> rev_b(solutionvectors[iFreq].col(iCart).data() + _alpha, _nv_b, _no_b);

            Eigen::MatrixXd xij = this->getLTDoubles(i, j, frequency, -1);
            Eigen::MatrixXd xxij = _sss * xij - _sss * xij.transpose();

            Dvv_b.noalias() += tij * xij.transpose();
            Dov_b.row(i).noalias() += xxij * gsl_b.col(j);
            if (i != j) {
              Dvv_b.noalias() += tij.transpose() * xij;
              Dov_b.row(j).noalias() += xxij.transpose() * gsl_b.col(i);
            }
          }
        }
      }

      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];

            Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[iFreq].col(iCart).data(), mo_a, mo_a);
            Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[iFreq].col(iCart).data() + mo_a * mo_a, mo_b, mo_b);

            Eigen::Ref<Eigen::MatrixXd> Dov_a = D_a.topRightCorner(_no_a, _nv_a);
            Eigen::Ref<Eigen::MatrixXd> Dvv_a = D_a.bottomRightCorner(_nv_a, _nv_a);
            Eigen::Ref<Eigen::MatrixXd> Dov_b = D_b.topRightCorner(_no_b, _nv_b);
            Eigen::Ref<Eigen::MatrixXd> Dvv_b = D_b.bottomRightCorner(_nv_b, _nv_b);

            Eigen::Map<Eigen::MatrixXd> rev_a(solutionvectors[iFreq].col(iCart).data(), _nv_a, _no_a);
            Eigen::Map<Eigen::MatrixXd> rev_b(solutionvectors[iFreq].col(iCart).data() + _alpha, _nv_b, _no_b);

            Eigen::MatrixXd xij = this->getLTDoubles(i, j, frequency, 0);

            Dvv_a.noalias() += tij * xij.transpose();
            Dov_a.row(i).noalias() += _oss * xij * gsl_b.col(j);
            Dvv_b.noalias() += tij.transpose() * xij;
            Dov_b.row(j).noalias() += _oss * xij.transpose() * gsl_a.col(i);
          }
        }
      }

      this->reorderTensor(_Fia.head(_alpha), _nv_a, _no_a);
      this->reorderTensor(_gsLagrange.head(_alpha), _nv_a, _no_a);
      this->reorderTensor(_Xia.alpha, _nv_a, _no_a);
      this->reorderTensor(_Yia.alpha, _nv_a, _no_a);
      this->reorderTensor(_Zia.alpha, _nv_a, _no_a);
      this->reorderTensor(_Jia->alpha, _nv_a, _no_a);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->alpha, _nv_a, _no_a);
      }
      this->reorderTensor(_Fia.tail(_beta), _nv_b, _no_b);
      this->reorderTensor(_gsLagrange.tail(_beta), _nv_b, _no_b);
      this->reorderTensor(_Xia.beta, _nv_b, _no_b);
      this->reorderTensor(_Yia.beta, _nv_b, _no_b);
      this->reorderTensor(_Zia.beta, _nv_b, _no_b);
      this->reorderTensor(_Jia->beta, _nv_b, _no_b);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->beta, _nv_b, _no_b);
      }

      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = a; b < _nv_a; ++b) {
          Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, 1);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];
            Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[iFreq].col(iCart).data(), mo_a, mo_a);
            Eigen::Ref<Eigen::MatrixXd> Doo_a = D_a.topLeftCorner(_no_a, _no_a);

            Eigen::MatrixXd xab = this->getLTDoublesV(a, b, frequency, 1);
            Doo_a.noalias() -= xab * tab.transpose();
            if (a != b) {
              Doo_a.noalias() -= xab.transpose() * tab;
            }
          }
        }
      }
      for (size_t a = 0; a < _nv_b; ++a) {
        for (size_t b = a; b < _nv_b; ++b) {
          Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, -1);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];
            Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[iFreq].col(iCart).data() + mo_a * mo_a, mo_b, mo_b);
            Eigen::Ref<Eigen::MatrixXd> Doo_b = D_b.topLeftCorner(_no_b, _no_b);

            Eigen::MatrixXd xab = this->getLTDoublesV(a, b, frequency, -1);
            Doo_b.noalias() -= xab * tab.transpose();
            if (a != b) {
              Doo_b.noalias() -= xab.transpose() * tab;
            }
          }
        }
      }
      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = 0; b < _nv_b; ++b) {
          Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, 0);
          for (size_t iFreq = 0; iFreq < nFreq; ++iFreq) {
            double frequency = frequencies[iFreq];
            Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[iFreq].col(iCart).data(), mo_a, mo_a);
            Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[iFreq].col(iCart).data() + mo_a * mo_a, mo_b, mo_b);
            Eigen::Ref<Eigen::MatrixXd> Doo_a = D_a.topLeftCorner(_no_a, _no_a);
            Eigen::Ref<Eigen::MatrixXd> Doo_b = D_b.topLeftCorner(_no_b, _no_b);

            Eigen::MatrixXd xab = this->getLTDoublesV(a, b, frequency, 0);
            Doo_a.noalias() -= xab * tab.transpose();
            Doo_b.noalias() -= xab.transpose() * tab;
          }
        }
      }

      this->reorderTensor(_Fia.head(_alpha), _no_a, _nv_a);
      this->reorderTensor(_gsLagrange.head(_alpha), _no_a, _nv_a);
      this->reorderTensor(_Xia.alpha, _no_a, _nv_a);
      this->reorderTensor(_Zia.alpha, _no_a, _nv_a);
      this->reorderTensor(_Jia->alpha, _no_a, _nv_a);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->alpha, _no_a, _nv_a);
      }
      this->reorderTensor(_Fia.tail(_beta), _no_b, _nv_b);
      this->reorderTensor(_gsLagrange.tail(_beta), _no_b, _nv_b);
      this->reorderTensor(_Xia.beta, _no_b, _nv_b);
      this->reorderTensor(_Zia.beta, _no_b, _nv_b);
      this->reorderTensor(_Jia->beta, _no_b, _nv_b);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->beta, _no_b, _nv_b);
      }
    }

    auto itEnd = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();
    printf(" %5i %12.3f\n", iRoot + 1, duration / 60.0);
  }
  Timings::timeTaken("CC2 -      Nonit. Pert. Dens.");

  Timings::takeTime("CC2 -         It. Pert. Dens.");
  printf("\n  Iterative part.\n\n");
  printf("  root  time(min)\n");
  printf(" ------------------");
  // Iterative contribution.
  for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
    double frequency = frequencies[iFreq];
    for (unsigned iCart = 0; iCart < nCart; ++iCart) {
      auto itStart = std::chrono::steady_clock::now();
      Eigen::Ref<Eigen::VectorXd> solutionvector = solutionvectors[iFreq].col(iCart);

      Eigen::Map<Eigen::MatrixXd> D_a(densityMatrices[iFreq].col(iCart).data(), mo_a, mo_a);
      Eigen::Map<Eigen::MatrixXd> D_b(densityMatrices[iFreq].col(iCart).data() + mo_a * mo_a, mo_b, mo_b);

      Eigen::Ref<Eigen::MatrixXd> Doo_a = D_a.topLeftCorner(_no_a, _no_a);
      Eigen::Ref<Eigen::MatrixXd> Dov_a = D_a.topRightCorner(_no_a, _nv_a);
      Eigen::Ref<Eigen::MatrixXd> Dvv_a = D_a.bottomRightCorner(_nv_a, _nv_a);

      Eigen::Ref<Eigen::MatrixXd> Doo_b = D_b.topLeftCorner(_no_b, _no_b);
      Eigen::Ref<Eigen::MatrixXd> Dov_b = D_b.topRightCorner(_no_b, _nv_b);
      Eigen::Ref<Eigen::MatrixXd> Dvv_b = D_b.bottomRightCorner(_nv_b, _nv_b);

      Eigen::Map<Eigen::MatrixXd> rev_a(solutionvector.data(), _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> rev_b(solutionvector.data() + _alpha, _nv_b, _no_b);

      // Start with easy contributions.
      Dov_a.noalias() += rev_a.transpose();
      Doo_a.noalias() -= rev_a.transpose() * gsl_a;
      Dvv_a.noalias() += gsl_a * rev_a.transpose();
      Dov_a.noalias() -= rev_a.transpose() * CC2virt_a;
      Dov_a.noalias() += CC2occ_a * rev_a.transpose();

      Dov_b.noalias() += rev_b.transpose();
      Doo_b.noalias() -= rev_b.transpose() * gsl_b;
      Dvv_b.noalias() += gsl_b * rev_b.transpose();
      Dov_b.noalias() -= rev_b.transpose() * CC2virt_b;
      Dov_b.noalias() += CC2occ_b * rev_b.transpose();

      this->performTransformation(_Xia.alpha, solutionvector.head(_alpha), false, 1);
      this->performTransformation(_Xia.beta, solutionvector.tail(_beta), false, -1);

      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          auto tij = this->getGLagrangeAmplitudes(i, j, 1);
          auto rij = this->getRightAmplitudesA(i, j, frequency, 1);
          auto rrij = _sss * rij - _sss * rij.transpose();

          Dvv_a.noalias() += tij * rij.transpose();
          Dov_a.row(i).noalias() += rrij * gsl_a.col(j);
          if (i != j) {
            Dvv_a.noalias() += tij.transpose() * rij;
            Dov_a.row(j).noalias() += rrij.transpose() * gsl_a.col(i);
          }
        }
      }
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          auto tij = this->getGLagrangeAmplitudes(i, j, -1);
          auto rij = this->getRightAmplitudesA(i, j, frequency, -1);
          auto rrij = _sss * rij - _sss * rij.transpose();

          Dvv_b.noalias() += tij * rij.transpose();
          Dov_b.row(i).noalias() += rrij * gsl_b.col(j);
          if (i != j) {
            Dvv_b.noalias() += tij.transpose() * rij;
            Dov_b.row(j).noalias() += rrij.transpose() * gsl_b.col(i);
          }
        }
      }
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          auto tij = this->getGLagrangeAmplitudes(i, j, 0);
          auto rij = this->getRightAmplitudesA(i, j, frequency, 0);

          Dvv_a.noalias() += tij * rij.transpose();
          Dov_a.row(i).noalias() += _oss * rij * gsl_b.col(j);
          Dvv_b.noalias() += tij.transpose() * rij;
          Dov_b.row(j).noalias() += _oss * rij.transpose() * gsl_a.col(i);
        }
      }

      this->reorderTensor(_Fia.head(_alpha), _nv_a, _no_a);
      this->reorderTensor(_gsLagrange.head(_alpha), _nv_a, _no_a);
      this->reorderTensor(_Jia->alpha, _nv_a, _no_a);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->alpha, _nv_a, _no_a);
      }
      this->reorderTensor(_Xia.alpha, _nv_a, _no_a);
      this->reorderTensor(_Zia.alpha, _nv_a, _no_a);

      this->reorderTensor(_Fia.tail(_beta), _nv_b, _no_b);
      this->reorderTensor(_gsLagrange.tail(_beta), _nv_b, _no_b);
      this->reorderTensor(_Jia->beta, _nv_b, _no_b);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->beta, _nv_b, _no_b);
      }
      this->reorderTensor(_Xia.beta, _nv_b, _no_b);
      this->reorderTensor(_Zia.beta, _nv_b, _no_b);

      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = a; b < _nv_a; ++b) {
          auto tab = this->getGLagrangeAmplitudesV(a, b, 1);
          auto rab = this->getRightAmplitudesAV(a, b, frequency, 1);

          Doo_a.noalias() -= rab * tab.transpose();
          if (a != b) {
            Doo_a.noalias() -= rab.transpose() * tab;
          }
        }
      }
      for (size_t a = 0; a < _nv_b; ++a) {
        for (size_t b = a; b < _nv_b; ++b) {
          auto tab = this->getGLagrangeAmplitudesV(a, b, -1);
          auto rab = this->getRightAmplitudesAV(a, b, frequency, -1);

          Doo_b.noalias() -= rab * tab.transpose();
          if (a != b) {
            Doo_b.noalias() -= rab.transpose() * tab;
          }
        }
      }
      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = 0; b < _nv_b; ++b) {
          auto tab = this->getGLagrangeAmplitudesV(a, b, 0);
          auto rab = this->getRightAmplitudesAV(a, b, frequency, 0);

          Doo_a.noalias() -= rab * tab.transpose();
          Doo_b.noalias() -= rab.transpose() * tab;
        }
      }

      // done with ab-wise amps, must order ints back
      this->reorderTensor(_Fia.head(_alpha), _no_a, _nv_a);
      this->reorderTensor(_gsLagrange.head(_alpha), _no_a, _nv_a);
      this->reorderTensor(_Zia.alpha, _no_a, _nv_a);
      this->reorderTensor(_Jia->alpha, _no_a, _nv_a);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->alpha, _no_a, _nv_a);
      }

      this->reorderTensor(_Fia.tail(_beta), _no_b, _nv_b);
      this->reorderTensor(_gsLagrange.tail(_beta), _no_b, _nv_b);
      this->reorderTensor(_Zia.beta, _no_b, _nv_b);
      this->reorderTensor(_Jia->beta, _no_b, _nv_b);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->beta, _no_b, _nv_b);
      }

      auto itEnd = std::chrono::steady_clock::now();
      double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();
      auto str = ((iFreq * nCart + iCart) % 5 == 0) ? "\n %5i %9.3f" : " %5i %9.3f";
      printf(str, iFreq * nCart + iCart + 1, duration / 60.0);
    } /* iCart */
  }   /* iFreq */
  printf("\n\n");
  Timings::timeTaken("CC2 -         It. Pert. Dens.");
} /* this->calcPerturbedDensities() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity