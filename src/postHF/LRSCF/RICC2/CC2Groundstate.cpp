/**
 * @file CC2Groundstate.cpp
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
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<>
void CC2Controller<RESTRICTED>::calculateGroundstate() {
  // matrices to store vectors for diis
  Eigen::MatrixXd sglStorage, resStorage;

  // timings
  std::chrono::steady_clock::time_point itStart = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point itEnd = std::chrono::steady_clock::now();

  printf("    Convergence threshold             : %-5.1e\n", _settings.conv);
  printf("    Maximum number of iterations      : %-5i\n\n", _settings.maxCycles);
  printf("    The zeroth iteration corresponds to MP2.\n\n");

  _iter = 0;
  printf("   it.   time (min)      gs energy      corr energy      res norm      t1 norm  \n");
  printf("  ------------------------------------------------------------------------------\n");
  do {
    if (_xwfModel == Options::LR_METHOD::CC2) {
      // update singles with residuals
      for (size_t i = 0; i < _no; ++i) {
        for (size_t a = 0; a < _nv; ++a) {
          _singles(i * _nv + a) += _residual(i * _nv + a) / (_e(i) - _e(_no + a));
        }
      }
      if (_iter > 0) {
        if (_iter <= _settings.diisStore) {
          resStorage.conservativeResize(_no * _nv, resStorage.cols() + 1);
          sglStorage.conservativeResize(_no * _nv, sglStorage.cols() + 1);
        }
        resStorage.col((_iter - 1) % _settings.diisStore) = _residual;
        sglStorage.col((_iter - 1) % _settings.diisStore) = _singles;

        Eigen::VectorXd rhs(Eigen::VectorXd::Zero(_iter + 1));
        Eigen::MatrixXd lhs(Eigen::MatrixXd::Zero(_iter + 1, _iter + 1));

        for (size_t iStored = 0; iStored < _iter; ++iStored) {
          lhs(iStored, _iter) = -1.0;
          lhs(_iter, iStored) = -1.0;
        }
        lhs.topLeftCorner(_iter, _iter) = resStorage.transpose() * resStorage;
        double norm = lhs.topLeftCorner(_iter, _iter).cwiseAbs().maxCoeff();
        lhs.topLeftCorner(_iter, _iter) *= 1.0 / norm;
        rhs(_iter) = -1.0;
        _singles = sglStorage * lhs.householderQr().solve(rhs).head(_iter);
      }
      this->transformIntegrals();
      this->calculateResidual();
    }
    else {
      // Yia only needs to be calculated once
      _Yia.setZero();
      if (_settings.ltconv == 0) {
        for (size_t i = 0; i < _no; ++i) {
          for (size_t j = i; j < _no; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j);
            _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Jia->middleRows(j * _nv, _nv);
            if (i != j) {
              _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Jia->middleRows(i * _nv, _nv);
            }
          }
        }
      }
      else {
        for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
          Eigen::VectorXd param(_nv * _no);
          for (size_t i = 0, ia = 0; i < _no; ++i) {
            for (size_t a = _no; a < _no + _nv; ++a, ++ia) {
              param(ia) = std::exp(-(_e(a) - _e(i)) * _roots(iRoot));
            }
          }
          Eigen::MatrixXd QP = (*_Jia).transpose() * param.asDiagonal() * (*_Jia);
          _Yia.noalias() -= _oss * _weights(iRoot) * param.asDiagonal() * (*_Jia) * QP;
        }
      }
      _P = _C;
      _H = _C;
    }

    // calculate correlation energy
    _corrEnergy = _Yia.cwiseProduct(*_Jia).sum() + _Fia.dot(_singles);
    if (_iter == 0) {
      _mp2Energy = _corrEnergy;
    }

    // iteration time
    itEnd = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();

    // print info
    printf("%5i %11.3f %18.10f %15.10f %14.10f %10.6f\n", _iter++, duration / 60.0, _hfEnergy + _corrEnergy,
           _corrEnergy, _residual.norm(), _singles.norm());

    itStart = std::chrono::steady_clock::now();

    if (_iter > _settings.maxCycles) {
      printf("\n  Could not converge t1 amplitudes.\n");
      break;
    }
  } while (_residual.norm() > _settings.conv && _xwfModel == Options::LR_METHOD::CC2);
  printf("\n    RHF energy (a.u.): %24.10f\n", _hfEnergy);
  printf("    Final MP2 energy (a.u.): %18.10f\n", _mp2Energy);
  if (_xwfModel == Options::LR_METHOD::CC2) {
    printf("    Final CC2 energy (a.u.): %18.10f\n", _corrEnergy);
  }
  printf("\n");
} /* this->calculateGroundstate() restricted */

template<>
void CC2Controller<UNRESTRICTED>::calculateGroundstate() {
  // matrices to store vectors for diis
  Eigen::MatrixXd sglStorage, resStorage;

  // timings
  std::chrono::steady_clock::time_point itStart = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point itEnd = std::chrono::steady_clock::now();

  printf("    Convergence threshold             : %-5.1e\n", _settings.conv);
  printf("    Maximum number of iterations      : %-5i\n\n", _settings.maxCycles);
  printf("    The zeroth iteration corresponds to MP2.\n\n");

  _iter = 0;
  printf("   it.   time (min)      gs energy      corr energy      res norm      t1 norm  \n");
  printf("  ------------------------------------------------------------------------------\n");
  do {
    if (_xwfModel == Options::LR_METHOD::CC2) {
      // update singles with residuals
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t a = 0; a < _nv_a; ++a) {
          _singles(i * _nv_a + a) += _residual(i * _nv_a + a) / (_e.alpha(i) - _e.alpha(_no_a + a));
        }
      }
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t a = 0; a < _nv_b; ++a) {
          _singles(_alpha + i * _nv_b + a) += _residual(_alpha + i * _nv_b + a) / (_e.beta(i) - _e.beta(_no_b + a));
        }
      }
      if (_iter > 0) {
        if (_iter <= _settings.diisStore) {
          resStorage.conservativeResize(_nDim, resStorage.cols() + 1);
          sglStorage.conservativeResize(_nDim, sglStorage.cols() + 1);
        }
        resStorage.col((_iter - 1) % _settings.diisStore) = _residual;
        sglStorage.col((_iter - 1) % _settings.diisStore) = _singles;

        Eigen::VectorXd rhs(Eigen::VectorXd::Zero(_iter + 1));
        Eigen::MatrixXd lhs(Eigen::MatrixXd::Zero(_iter + 1, _iter + 1));

        for (size_t iStored = 0; iStored < _iter; ++iStored) {
          lhs(iStored, _iter) = -1.0;
          lhs(_iter, iStored) = -1.0;
        }
        lhs.topLeftCorner(_iter, _iter) = resStorage.transpose() * resStorage;
        double norm = lhs.topLeftCorner(_iter, _iter).cwiseAbs().maxCoeff();
        lhs.topLeftCorner(_iter, _iter) *= 1.0 / norm;
        rhs(_iter) = -1.0;
        _singles = sglStorage * lhs.householderQr().solve(rhs).head(_iter);
      }
      this->transformIntegrals();
      this->calculateResidual();
    }
    else {
      if (_settings.ltconv == 0) {
        // alpha
        for (size_t i = 0; i < _no_a; ++i) {
          for (size_t j = i; j < _no_a; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
            _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
            if (i != j) {
              _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
            }
          }
        }

        // beta
        for (size_t i = 0; i < _no_b; ++i) {
          for (size_t j = i; j < _no_b; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
            _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
            if (i != j) {
              _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
            }
          }
        }

        // mixed
        for (size_t i = 0; i < _no_a; ++i) {
          for (size_t j = 0; j < _no_b; ++j) {
            Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
            _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
            _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
          }
        }
      }
      else {
        for (int iRoot = 0; iRoot < _roots.size(); ++iRoot) {
          Eigen::VectorXd param_a(_nv_a * _no_a);
          for (size_t i = 0, ia = 0; i < _no_a; ++i) {
            for (size_t a = _no_a; a < _no_a + _nv_a; ++a, ++ia) {
              param_a(ia) = std::exp(-(_e.alpha(a) - _e.alpha(i)) * _roots(iRoot));
            }
          }
          Eigen::VectorXd param_b(_nv_b * _no_b);
          for (size_t j = 0, jb = 0; j < _no_b; ++j) {
            for (size_t b = _no_b; b < _no_b + _nv_b; ++b, ++jb) {
              param_b(jb) = std::exp(-(_e.beta(b) - _e.beta(j)) * _roots(iRoot));
            }
          }
          Eigen::MatrixXd QP_a = _Jia->alpha.transpose() * param_a.asDiagonal() * _Jia->alpha;
          Eigen::MatrixXd QP_b = _Jia->beta.transpose() * param_b.asDiagonal() * _Jia->beta;

          _Yia.alpha.noalias() -= _oss * _weights(iRoot) * param_a.asDiagonal() * _Jia->alpha * QP_b;
          _Yia.beta.noalias() -= _oss * _weights(iRoot) * param_b.asDiagonal() * _Jia->beta * QP_a;
        }
      }
      for_spin(_C, _P, _H) {
        _P_spin = _C_spin;
        _H_spin = _C_spin;
      };
    }

    // calculate correlation energy
    _corrEnergy = 0.5 * _Fia.dot(_singles);
    auto& Jia = *_Jia;
    for_spin(_Yia, Jia) {
      _corrEnergy += 0.5 * _Yia_spin.cwiseProduct(Jia_spin).sum();
    };
    if (_iter == 0) {
      _mp2Energy = _corrEnergy;
    }

    // iteration time
    itEnd = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();

    // print info
    printf("%5i %11.3f %18.10f %15.10f %14.10f %10.6f\n", _iter, duration / 60.0, _hfEnergy + _corrEnergy, _corrEnergy,
           _residual.norm(), _singles.norm());

    itStart = std::chrono::steady_clock::now();

    if (_iter > _settings.maxCycles) {
      printf("\n  Could not converge t1 amplitudes.\n");
      break;
    }
    ++_iter;
  } while (_residual.norm() > _settings.conv && _xwfModel == Options::LR_METHOD::CC2);
  printf("\n    UHF energy (a.u.): %24.10f\n", _hfEnergy);
  printf("    Final MP2 energy (a.u.): %18.10f\n", _mp2Energy);
  if (_xwfModel == Options::LR_METHOD::CC2) {
    printf("    Final CC2 energy (a.u.): %18.10f\n", _corrEnergy);
  }
  printf("\n");
} /* this->calculateGroundstate() unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity