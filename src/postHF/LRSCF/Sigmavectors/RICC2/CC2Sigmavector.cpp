/**
 * @file CC2Sigmavector.cpp
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
#include "postHF/LRSCF/Sigmavectors/RICC2/CC2Sigmavector.h"
#include "integrals/looper/TwoElecThreeCenterIntLooper.h"
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CC2Sigmavector<SCFMode>::CC2Sigmavector(std::shared_ptr<LRSCFController<SCFMode>> lrscf)
  : XWFController<SCFMode>(lrscf) {
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd CC2Sigmavector<SCFMode>::transformedABFock(Eigen::Ref<Eigen::VectorXd> trafovector,
                                                           Eigen::Ref<Eigen::VectorXd> XQ, Eigen::MatrixXd& Yia, int iSpin) {
  // declarations
  unsigned _nv = (iSpin > 0) ? this->_nv_a : this->_nv_b;
  unsigned _no = (iSpin > 0) ? this->_no_a : this->_no_b;
  double* coeffAPtr;
  double* coeffBPtr;
  double* JiaPtr;

  auto& P = this->_P;
  auto& H = this->_H;
  auto& Jia = *this->_Jia;

  int fSpin = 1;
  for_spin(P, H, Jia) {
    if (iSpin * fSpin > 0) {
      coeffAPtr = P_spin.data();
      coeffBPtr = H_spin.data();
      JiaPtr = Jia_spin.data();
    }
    fSpin *= -1;
  };

  Eigen::Map<Eigen::MatrixXd> coeffA(coeffAPtr + this->_nb * _no, this->_nb, _nv);
  Eigen::Map<Eigen::MatrixXd> coeffB(coeffBPtr + this->_nb * _no, this->_nb, _nv);

  Eigen::MatrixXd Fab = Eigen::MatrixXd::Zero(_nv, _nv);
  Eigen::Map<Eigen::MatrixXd> trafo(trafovector.data(), _nv, _no);

  /*
   * calculates following contributions to Fab block in transformed Fock
   * _Eab += 2 * bab * XQ(Q);
   * _Eab -= bab * sgl * jia.transpose();
   */
  Eigen::MatrixXd eta = coeffB * trafo;
  Yia = Eigen::MatrixXd::Zero(this->_nb * _no, this->_nxb);

  Eigen::VectorXd XP = this->_riTrafo * XQ;
  unsigned nThreads = omp_get_max_threads();
  Eigen::ArrayXd pseudoFock = Eigen::ArrayXd::Zero(nThreads * this->_nb * this->_nb);

  double* YiaPtr = Yia.data();
  double* etaptr = eta.data();
  double* fockptr = pseudoFock.data();

  // loop over all cached (mn|P)
  double* jmnptr = this->_Jmn->data();
#pragma omp parallel
  {
    size_t indexP;
    size_t offset1 = omp_get_thread_num() * this->_nb * this->_nb;
#pragma omp for schedule(dynamic)
    for (size_t P = 0; P < this->_nxs; ++P) {
      size_t offset2 = P * this->_nb * _no;
      indexP = P * this->_nb * (this->_nb + 1) / 2;
      for (size_t nu = 0; nu < this->_nb; ++nu) {
        for (size_t mu = nu; mu < this->_nb; ++mu, ++indexP) {
          double& integral = jmnptr[indexP];
          if (std::abs(integral) > this->_prescreeningThreshold) {
            fockptr[offset1 + nu * this->_nb + mu] += integral * XP(P);
            if (mu != nu) {
              fockptr[offset1 + mu * this->_nb + nu] += integral * XP(P);
            }
            for (size_t i = 0; i < _no; ++i) {
              YiaPtr[offset2 + i * this->_nb + mu] += etaptr[i * this->_nb + nu] * integral;
              if (mu != nu) {
                YiaPtr[offset2 + i * this->_nb + nu] += etaptr[i * this->_nb + mu] * integral;
              }
            } /* i */
          }   /* prescreen */
        }     /* mu */
      }       /* nu */
    }         /* P */
  }           /* parallel region */

  // loop over the rest
  if (this->_nxs < this->_nxb) {
    auto loopFunc = [&](unsigned mu, unsigned nu, size_t P, double integral, size_t threadId) {
      size_t offset1 = threadId * this->_nb * this->_nb;
      fockptr[offset1 + nu * this->_nb + mu] += integral * XP(P);
      if (mu != nu) {
        fockptr[offset1 + mu * this->_nb + nu] += integral * XP(P);
      }

      size_t offset2 = P * this->_nb * _no;
      for (unsigned i = 0; i < _no; ++i) {
        YiaPtr[offset2 + i * this->_nb + mu] += etaptr[i * this->_nb + nu] * integral;
        if (mu != nu) {
          YiaPtr[offset2 + i * this->_nb + nu] += etaptr[i * this->_nb + mu] * integral;
        }
      }
    };
    this->_looper->loopNoDerivative(loopFunc);
  }

  this->performRITransformation(Yia, _no);

  unsigned spinFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;

  // sum over all threads
  for (size_t iThread = 0; iThread < nThreads; ++iThread) {
    Eigen::Map<Eigen::MatrixXd> psf(fockptr + iThread * this->_nb * this->_nb, this->_nb, this->_nb);
    Fab += spinFactor * coeffA.transpose() * psf * coeffB;
  }

  for (size_t Q = 0; Q < this->_nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(JiaPtr + Q * _nv * _no, _nv, _no);
    Eigen::Map<Eigen::MatrixXd> yia(YiaPtr + Q * this->_nb * _no, this->_nb, _no);
    Fab -= coeffA.transpose() * yia * jia.transpose();
  }

  return Fab;
} /* this->transformedABFock() restricted / unrestricted */

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd CC2Sigmavector<SCFMode>::lagrangianSubspaceSolve(Eigen::Ref<Eigen::MatrixXd> rightHandSides,
                                                                 Eigen::Ref<Eigen::VectorXd> eigenvalues) {
  unsigned nEigen = eigenvalues.size();

  printf("\n    Number of roots to be determined  : %-5i\n", nEigen);
  printf("    Convergence threshold             : %-5.1e\n", this->_settings.conv);
  printf("    Maximum number of iterations      : %-5i\n\n", this->_settings.maxCycles);
  printTableHead("root     it.   dimension    converged     max norm");

  // here you cannot use a common subspace for all states
  // since the reference matrix depends on the excitation energy
  // that's why I solve each root one after another

  Eigen::MatrixXd guessVectors(this->_nDim, 1);
  Eigen::MatrixXd sigmaVectors(this->_nDim, 1);
  Eigen::MatrixXd ritzVectors(this->_nDim, nEigen);
  Eigen::VectorXd residual(this->_nDim);

  Eigen::MatrixXd lhs, rhs, sol;

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    // Set references
    Eigen::Ref<Eigen::VectorXd> rightHandSide = rightHandSides.col(iEigen);
    double eigenvalue = eigenvalues(iEigen);

    // Reset subspace
    guessVectors.resize(this->_nDim, 1);
    sigmaVectors.resize(this->_nDim, 1);

    // Setup initial subspace
    auto& e = this->_e;
    auto& nv_ref = this->_nv;
    auto& no_ref = this->_no;
    unsigned iStart = 0;
    for_spin(e, nv_ref, no_ref) {
      for (size_t i = 0; i < no_ref_spin; ++i) {
        for (size_t a = 0; a < nv_ref_spin; ++a) {
          guessVectors(iStart + i * nv_ref_spin + a, 0) =
              rightHandSide(iStart + i * nv_ref_spin + a) / (e_spin(no_ref_spin + a) - e_spin(i) + eigenvalue);
        }
      }
      iStart += nv_ref_spin * no_ref_spin;
    };
    guessVectors.col(0).normalize();
    sigmaVectors.col(0) = this->getLeftXWFSigma(guessVectors.col(0), -eigenvalue) + eigenvalue * guessVectors.col(0);

    unsigned _iter = 1;
    while (true) {
      // setup subspace problem
      lhs = guessVectors.transpose() * sigmaVectors;
      rhs = guessVectors.transpose() * rightHandSide;
      // solve subspace problem
      sol = lhs.householderQr().solve(rhs);
      // ritz vector
      ritzVectors.col(iEigen) = guessVectors * sol;
      // residual vectors
      residual = sigmaVectors * sol - rightHandSide;
      // print info
      printf("%5i %7i %9i %11s %16.2e\n", iEigen + 1, this->_iter, (int)lhs.rows(),
             (residual.norm() < this->_settings.conv) ? "T" : "", (double)residual.norm());
      if (residual.norm() < this->_settings.conv || _iter >= this->_settings.maxCycles) {
        break;
      }

      // precondition using zeroth order Jacobian
      iStart = 0;
      for_spin(e, nv_ref, no_ref) {
        for (size_t i = 0; i < no_ref_spin; ++i) {
          for (size_t a = 0; a < nv_ref_spin; ++a) {
            residual(iStart + i * nv_ref_spin + a) =
                residual(iStart + i * nv_ref_spin + a) / (e_spin(no_ref_spin + a) - e_spin(i) + eigenvalue);
          }
        }
        iStart += nv_ref_spin * no_ref_spin;
      };

      // orthogonalize against guess space with Gram-Schmidt
      residual.normalize();
      for (unsigned g = 0; g < guessVectors.cols(); ++g) {
        double ij = residual.dot(guessVectors.col(g));
        double jj = guessVectors.col(g).dot(guessVectors.col(g));
        residual -= ij / jj * guessVectors.col(g);
      }
      residual.normalize();

      // subspace expansion
      guessVectors.conservativeResize(this->_nDim, guessVectors.cols() + 1);
      sigmaVectors.conservativeResize(this->_nDim, sigmaVectors.cols() + 1);
      guessVectors.rightCols(1) = residual;
      sigmaVectors.rightCols(1) = this->getLeftXWFSigma(residual, -eigenvalue) + eigenvalue * residual;
      ++_iter;
    }
  }

  return ritzVectors;
} /* this->lagrangianSubspaceSolve() restricted / unrestricted */

template<>
void CC2Sigmavector<RESTRICTED>::calculateE() {
  _Eij = Eigen::MatrixXd::Zero(_no, _no);
  _Eab = Eigen::MatrixXd::Zero(_nv, _nv);

  _Eij.diagonal() = _e.segment(0, _no);
  _Eab.diagonal() = _e.segment(_no, _nv);

  Eigen::VectorXd XQ = _Jia->transpose() * _singles;
  Eigen::Map<Eigen::MatrixXd> sgl(_singles.data(), _nv, _no);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
    Eigen::Map<Eigen::MatrixXd> yia(_Yia.data() + Q * _nv * _no, _nv, _no);
    // ^Fij and ^Fab
    _Eij += 2 * bij * XQ(Q);

    _Eij -= jia.transpose() * sgl * bij;

    _Eij += jia.transpose() * yia;
    _Eab -= yia * jia.transpose();
  }

  _Eab += this->transformedABFock(_singles, XQ, _Yia);

  // don't need Yia anymore now
  _Yia.resize(0, 0);

  // gonna need these though
  _Xia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);
  _Zia = SpinPolarizedData<RESTRICTED, Eigen::MatrixXd>(_no * _nv, _nx);
} /* this->calculateE() restricted */

template<>
void CC2Sigmavector<UNRESTRICTED>::calculateE() {
  Eigen::Map<Eigen::MatrixXd> sgl_a(_singles.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sgl_b(_singles.data() + _alpha, _nv_b, _no_b);

  _Eij.alpha = Eigen::MatrixXd::Zero(_no_a, _no_a);
  _Eij.beta = Eigen::MatrixXd::Zero(_no_b, _no_b);
  _Eab.alpha = Eigen::MatrixXd::Zero(_nv_a, _nv_a);
  _Eab.beta = Eigen::MatrixXd::Zero(_nv_b, _nv_b);

  _Eij.alpha.diagonal() = _e.alpha.segment(0, _no_a);
  _Eij.beta.diagonal() = _e.beta.segment(0, _no_b);
  _Eab.alpha.diagonal() = _e.alpha.segment(_no_a, _nv_a);
  _Eab.beta.diagonal() = _e.beta.segment(_no_b, _nv_b);

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * _singles.head(_alpha) + _Jia->beta.transpose() * _singles.tail(_beta);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    // ^Fij and ^Fab
    _Eij.alpha += bij_a * XQ(Q);

    _Eij.alpha -= jia_a.transpose() * sgl_a * bij_a;

    _Eij.alpha += jia_a.transpose() * yia_a;
    _Eab.alpha -= yia_a * jia_a.transpose();

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    // ^Fij and ^Fab
    _Eij.beta += bij_b * XQ(Q);

    _Eij.beta -= jia_b.transpose() * sgl_b * bij_b;

    _Eij.beta += jia_b.transpose() * yia_b;
    _Eab.beta -= yia_b * jia_b.transpose();
  }

  _Eab.alpha += this->transformedABFock(_singles.head(_alpha), XQ, _Yia.alpha, 1);
  _Eab.beta += this->transformedABFock(_singles.tail(_beta), XQ, _Yia.beta, -1);

  _Xia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  _Zia = SpinPolarizedData<UNRESTRICTED, Eigen::MatrixXd>();
  for_spin(_Xia, _Yia, _Zia, _nv, _no) {
    _Yia_spin.resize(0, 0);
    _Xia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
    _Zia_spin = Eigen::MatrixXd::Zero(_nv_spin * _no_spin, _nx);
  };
} /* this->calculateE() unrestricted */

template<>
void CC2Sigmavector<RESTRICTED>::calculateResidual() {
  _Fia = Eigen::VectorXd::Zero(_nDim);
  _residual = Eigen::VectorXd::Zero(_nDim);

  Eigen::Map<Eigen::MatrixXd> sgl(_singles.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> res(_residual.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> fia(_Fia.data(), _nv, _no);

  Eigen::VectorXd XQ = _Jia->transpose() * _singles;
  _Fia = 2 * _Jia->operator*(XQ);
  _residual = 2 * _Bia->operator*(XQ);

  res.noalias() += _e.segment(_no, _nv).asDiagonal() * sgl;
  res.noalias() -= sgl * _e.segment(0, _no).asDiagonal();

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _no * _nv, _nv, _no);
    fia.noalias() -= jia * sgl.transpose() * jia;
  }

  _Yia.setZero();
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j);
      _residual.segment(i * _nv, _nv).noalias() += tij * _Fia.segment(j * _nv, _nv);
      _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Jia->middleRows(j * _nv, _nv);
      if (i != j) {
        _residual.segment(j * _nv, _nv).noalias() += tij.transpose() * _Fia.segment(i * _nv, _nv);
        _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Jia->middleRows(i * _nv, _nv);
      }
    }
  }

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> yia(_Yia.data() + Q * _no * _nv, _nv, _no);
    res.noalias() -= yia * bij;
  }

  _residual.noalias() += this->getJ2GContribution(_Yia.data(), _Bij->data(), _singles, false);
} /* this->calculateResidual() restricted */

template<>
void CC2Sigmavector<UNRESTRICTED>::calculateResidual() {
  _Fia = Eigen::VectorXd::Zero(_nDim);
  _residual = Eigen::VectorXd::Zero(_nDim);

  Eigen::Map<Eigen::MatrixXd> sgl_a(_singles.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> res_a(_residual.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> fia_a(_Fia.data(), _nv_a, _no_a);

  Eigen::Map<Eigen::MatrixXd> sgl_b(_singles.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> res_b(_residual.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> fia_b(_Fia.data() + _alpha, _nv_b, _no_b);

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * _singles.head(_alpha) + _Jia->beta.transpose() * _singles.tail(_beta);
  _Fia.head(_alpha) += _Jia->alpha * XQ;
  _Fia.tail(_beta) += _Jia->beta * XQ;
  _residual.head(_alpha) += _Bia->alpha * XQ;
  _residual.tail(_beta) += _Bia->beta * XQ;

  res_a.noalias() += _e.alpha.segment(_no_a, _nv_a).asDiagonal() * sgl_a;
  res_a.noalias() -= sgl_a * _e.alpha.segment(0, _no_a).asDiagonal();

  res_b.noalias() += _e.beta.segment(_no_b, _nv_b).asDiagonal() * sgl_b;
  res_b.noalias() -= sgl_b * _e.beta.segment(0, _no_b).asDiagonal();

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _no_a * _nv_a, _nv_a, _no_a);
    fia_a.noalias() -= jia_a * sgl_a.transpose() * jia_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _no_b * _nv_b, _nv_b, _no_b);
    fia_b.noalias() -= jia_b * sgl_b.transpose() * jia_b;
  }
  _Yia.alpha.setZero();
  _Yia.beta.setZero();
  // alpha
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
      _residual.segment(i * _nv_a, _nv_a).noalias() += tij * _Fia.segment(j * _nv_a, _nv_a);
      _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
      if (i != j) {
        _residual.segment(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
        _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }

  // beta
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
      _residual.segment(_alpha + i * _nv_b, _nv_b).noalias() += tij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
      _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
      if (i != j) {
        _residual.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * _Fia.segment(_alpha + i * _nv_b, _nv_b);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
      }
    }
  }

  // mixed
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
      _residual.segment(i * _nv_a, _nv_a).noalias() += tij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
      _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Jia->beta.middleRows(j * _nv_b, _nv_b);

      _residual.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
      _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
    }
  }

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    res_a.noalias() -= yia_a * bij_a;

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    res_b.noalias() -= yia_b * bij_b;
  }

  _residual.head(_alpha).noalias() +=
      this->getJ2GContribution(_Yia.alpha.data(), _Bij->alpha.data(), _singles.head(_alpha), false, 1);
  _residual.tail(_beta).noalias() +=
      this->getJ2GContribution(_Yia.beta.data(), _Bij->beta.data(), _singles.tail(_beta), false, -1);
} /* this->calculateResidual() unrestricted */

template<>
Eigen::VectorXd CC2Sigmavector<RESTRICTED>::getRightXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  _trafoVector = guessVector;
  this->performTransformation(_Xia, guessVector, false);

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);
  Eigen::VectorXd Fia = Eigen::VectorXd::Zero(_nDim);

  Eigen::VectorXd XQ = _Jia->transpose() * guessVector;
  sigmaVector = 2 * _Bia->operator*(XQ);
  Fia = 2 * _Jia->operator*(XQ);

  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> fia(Fia.data(), _nv, _no);

  sigma.noalias() += _Eab * guess - guess * _Eij;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
    fia.noalias() -= jia * guess.transpose() * jia;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  Timings::takeTime("Exc. State WF -  Amp-Contraction");
  _Zia.setZero();
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j);
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue);
      sigmaVector.segment(i * _nv, _nv).noalias() += tij * Fia.segment(j * _nv, _nv);
      sigmaVector.segment(i * _nv, _nv).noalias() += rij * _Fia.segment(j * _nv, _nv);
      _Zia.middleRows(i * _nv, _nv).noalias() += rij * _Jia->middleRows(j * _nv, _nv);
      if (i != j) {
        sigmaVector.segment(j * _nv, _nv).noalias() += tij.transpose() * Fia.segment(i * _nv, _nv);
        sigmaVector.segment(j * _nv, _nv).noalias() += rij.transpose() * _Fia.segment(i * _nv, _nv);
        _Zia.middleRows(j * _nv, _nv).noalias() += rij.transpose() * _Jia->middleRows(i * _nv, _nv);
      }
    }
  }
  Timings::timeTaken("Exc. State WF -  Amp-Contraction");

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> zia(_Zia.data() + Q * _nv * _no, _nv, _no);
    sigma.noalias() -= zia * bij;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  sigmaVector += this->getJ2GContribution(_Zia.data(), _Bij->data(), guessVector, false);

  return sigmaVector;
} /* this->getRightXWFSigma() restricted */

template<>
Eigen::VectorXd CC2Sigmavector<UNRESTRICTED>::getRightXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  this->performTransformation(_Xia.alpha, guessVector.head(_alpha), false, 1);
  this->performTransformation(_Xia.beta, guessVector.tail(_beta), false, -1);

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);
  Eigen::VectorXd Fia = Eigen::VectorXd::Zero(_nDim);

  Eigen::VectorXd XQ = _Jia->alpha.transpose() * guessVector.head(_alpha) + _Jia->beta.transpose() * guessVector.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Bia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += _Bia->beta * XQ;
  Fia.head(_alpha).noalias() += _Jia->alpha * XQ;
  Fia.tail(_beta).noalias() += _Jia->beta * XQ;

  Eigen::Map<Eigen::MatrixXd> guess_a(guessVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> fia_a(Fia.data(), _nv_a, _no_a);

  Eigen::Map<Eigen::MatrixXd> guess_b(guessVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> fia_b(Fia.data() + _alpha, _nv_b, _no_b);

  sigma_a.noalias() += _Eab.alpha * guess_a - guess_a * _Eij.alpha;
  sigma_b.noalias() += _Eab.beta * guess_b - guess_b * _Eij.beta;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _no_a * _nv_a, _nv_a, _no_a);
    fia_a.noalias() -= jia_a * guess_a.transpose() * jia_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _no_b * _nv_b, _nv_b, _no_b);
    fia_b.noalias() -= jia_b * guess_b.transpose() * jia_b;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  Timings::takeTime("Exc. State WF -  Amp-Contraction");
  _Zia.alpha.setZero();
  _Zia.beta.setZero();
  // alpha
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 1);
      sigmaVector.segment(i * _nv_a, _nv_a).noalias() += tij * Fia.segment(j * _nv_a, _nv_a);
      sigmaVector.segment(i * _nv_a, _nv_a).noalias() += rij * _Fia.segment(j * _nv_a, _nv_a);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
      if (i != j) {
        sigmaVector.segment(j * _nv_a, _nv_a).noalias() += tij.transpose() * Fia.segment(i * _nv_a, _nv_a);
        sigmaVector.segment(j * _nv_a, _nv_a).noalias() += rij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
        _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }

  // beta
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, -1);
      sigmaVector.segment(_alpha + i * _nv_b, _nv_b).noalias() += tij * Fia.segment(_alpha + j * _nv_b, _nv_b);
      sigmaVector.segment(_alpha + i * _nv_b, _nv_b).noalias() += rij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
      _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
      if (i != j) {
        sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * Fia.segment(_alpha + i * _nv_b, _nv_b);
        sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += rij.transpose() * _Fia.segment(_alpha + i * _nv_b, _nv_b);
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
      }
    }
  }

  // mixed
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
      Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 0);

      sigmaVector.segment(i * _nv_a, _nv_a).noalias() += tij * Fia.segment(_alpha + j * _nv_b, _nv_b);
      sigmaVector.segment(i * _nv_a, _nv_a).noalias() += rij * _Fia.segment(_alpha + j * _nv_b, _nv_b);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);

      sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * Fia.segment(i * _nv_a, _nv_a);
      sigmaVector.segment(_alpha + j * _nv_b, _nv_b).noalias() += rij.transpose() * _Fia.segment(i * _nv_a, _nv_a);
      _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
    }
  }

  Timings::takeTime("Exc. State WF -    Q-Contraction");

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    sigma_a.noalias() -= zia_a * bij_a;

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    sigma_b.noalias() -= zia_b * bij_b;
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  sigmaVector.head(_alpha).noalias() +=
      this->getJ2GContribution(_Zia.alpha.data(), _Bij->alpha.data(), guessVector.head(_alpha), false, 1);
  sigmaVector.tail(_beta).noalias() +=
      this->getJ2GContribution(_Zia.beta.data(), _Bij->beta.data(), guessVector.tail(_beta), false, -1);

  return sigmaVector;
} /* this->getRightXWFSigma() unrestricted */

template<>
Eigen::VectorXd CC2Sigmavector<RESTRICTED>::getLeftXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  _trafoVector = guessVector;
  this->performTransformation(_Zia, guessVector, true);

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  Eigen::VectorXd sigmaVector = 2 * _Jia->operator*(_Bia->transpose() * guessVector);

  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);

  sigma.noalias() += _Eab.transpose() * guess - guess * _Eij.transpose();
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  Timings::takeTime("Exc. State WF -  Amp-Contraction");
  _Xia.setZero();
  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_no * _nv);
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j);
      Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue);
      Tia.segment(i * _nv, _nv).noalias() += tij * guessVector.segment(j * _nv, _nv);
      _Xia.middleRows(i * _nv, _nv).noalias() += lij * _Bia->middleRows(j * _nv, _nv);
      if (i != j) {
        Tia.segment(j * _nv, _nv).noalias() += tij.transpose() * guessVector.segment(i * _nv, _nv);
        _Xia.middleRows(j * _nv, _nv).noalias() += lij.transpose() * _Bia->middleRows(i * _nv, _nv);
      }
    }
  }
  Timings::timeTaken("Exc. State WF -  Amp-Contraction");

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  sigmaVector.noalias() += 2 * _Jia->operator*(_Jia->transpose() * Tia);
  Eigen::Map<Eigen::MatrixXd> tia(Tia.data(), _nv, _no);

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
    Eigen::Map<Eigen::MatrixXd> xia(_Xia.data() + Q * _nv * _no, _nv, _no);
    sigma.noalias() -= jia * tia.transpose() * jia;
    sigma.noalias() -= xia * bij.transpose();
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  sigmaVector += this->getJ2GContribution(_Xia.data(), _Bij->data(), guessVector, true);

  return sigmaVector;
} /* this->getLeftXWFSigma() restricted */

template<>
Eigen::VectorXd CC2Sigmavector<UNRESTRICTED>::getLeftXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) {
  _trafoVector = guessVector;
  this->performTransformation(_Zia.alpha, guessVector.head(_alpha), true, 1);
  this->performTransformation(_Zia.beta, guessVector.tail(_beta), true, -1);

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nDim);
  Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nDim);

  Eigen::VectorXd XQ = _Bia->alpha.transpose() * guessVector.head(_alpha) + _Bia->beta.transpose() * guessVector.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Jia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += _Jia->beta * XQ;

  Eigen::Map<Eigen::MatrixXd> guess_a(guessVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sigma_a(sigmaVector.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> tia_a(Tia.data(), _nv_a, _no_a);

  Eigen::Map<Eigen::MatrixXd> guess_b(guessVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> sigma_b(sigmaVector.data() + _alpha, _nv_b, _no_b);
  Eigen::Map<Eigen::MatrixXd> tia_b(Tia.data() + _alpha, _nv_b, _no_b);

  sigma_a.noalias() += _Eab.alpha.transpose() * guess_a - guess_a * _Eij.alpha.transpose();
  sigma_b.noalias() += _Eab.beta.transpose() * guess_b - guess_b * _Eij.beta.transpose();
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  Timings::takeTime("Exc. State WF -  Amp-Contraction");
  _Xia.alpha.setZero();
  _Xia.beta.setZero();
  // alpha
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, 1);
      Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 1);
      Tia.segment(i * _nv_a, _nv_a).noalias() += tij * guessVector.segment(j * _nv_a, _nv_a);
      _Xia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += lij * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
      if (i != j) {
        Tia.segment(j * _nv_a, _nv_a).noalias() += tij.transpose() * guessVector.segment(i * _nv_a, _nv_a);
        _Xia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += lij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }

  // beta
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, -1);
      Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, -1);
      Tia.segment(_alpha + i * _nv_b, _nv_b).noalias() += tij * guessVector.segment(_alpha + j * _nv_b, _nv_b);
      _Xia.beta.middleRows(i * _nv_b, _nv_b).noalias() += lij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
      if (i != j) {
        Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * guessVector.segment(_alpha + i * _nv_b, _nv_b);
        _Xia.beta.middleRows(j * _nv_b, _nv_b).noalias() += lij.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
      }
    }
  }

  // mixed
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      Eigen::MatrixXd tij = this->getAmplitudes(i, j, 0);
      Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 0);

      Tia.segment(i * _nv_a, _nv_a).noalias() += tij * guessVector.segment(_alpha + j * _nv_b, _nv_b);
      _Xia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += lij * _Bia->beta.middleRows(j * _nv_b, _nv_b);

      Tia.segment(_alpha + j * _nv_b, _nv_b).noalias() += tij.transpose() * guessVector.segment(i * _nv_a, _nv_a);
      _Xia.beta.middleRows(j * _nv_b, _nv_b).noalias() += lij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
    }
  }

  Timings::takeTime("Exc. State WF -    Q-Contraction");
  XQ = _Jia->alpha.transpose() * Tia.head(_alpha) + _Jia->beta.transpose() * Tia.tail(_beta);
  sigmaVector.head(_alpha).noalias() += _Jia->alpha * XQ;
  sigmaVector.tail(_beta).noalias() += _Jia->beta * XQ;

  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> xia_a(_Xia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    sigma_a.noalias() -= jia_a * tia_a.transpose() * jia_a;
    sigma_a.noalias() -= xia_a * bij_a.transpose();

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> xia_b(_Xia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    sigma_b.noalias() -= jia_b * tia_b.transpose() * jia_b;
    sigma_b.noalias() -= xia_b * bij_b.transpose();
  }
  Timings::timeTaken("Exc. State WF -    Q-Contraction");

  sigmaVector.head(_alpha).noalias() +=
      this->getJ2GContribution(_Xia.alpha.data(), _Bij->alpha.data(), guessVector.head(_alpha), true, 1);
  sigmaVector.tail(_beta).noalias() +=
      this->getJ2GContribution(_Xia.beta.data(), _Bij->beta.data(), guessVector.tail(_beta), true, -1);

  return sigmaVector;
} /* this->getLeftXWFSigma() unrestricted */

template<>
void CC2Sigmavector<RESTRICTED>::calcDensityMatrices(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues,
                                                     std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("Exc. State WF - Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned mo = _no + _nv;
  densityMatrices[0] = Eigen::MatrixXd::Zero(mo * mo, nEigen);
  densityMatrices[1] = Eigen::MatrixXd::Zero(mo * mo, nEigen);
  printf("  Calculating %3i eta and %3i xi density matrices.\n\n", nEigen, 2 * nEigen);

  // cc2 density matrix intermediate
  Eigen::MatrixXd CC2virt = Eigen::MatrixXd::Zero(_nv, _nv);
  Eigen::MatrixXd CC2occ = Eigen::MatrixXd::Zero(_no, _no);

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    // will need two density matrix objects
    // for left and right transition moments
    //  rightD = D^eta(R) + D^xi(M)
    //   leftD = D^xi(L)
    Eigen::Map<Eigen::MatrixXd> rightD(densityMatrices[0].col(iEigen).data(), mo, mo);
    Eigen::Map<Eigen::MatrixXd> leftD(densityMatrices[1].col(iEigen).data(), mo, mo);

    // reference blocks
    Eigen::Ref<Eigen::MatrixXd> rightDoo = rightD.topLeftCorner(_no, _no);
    Eigen::Ref<Eigen::MatrixXd> rightDov = rightD.topRightCorner(_no, _nv);
    Eigen::Ref<Eigen::MatrixXd> rightDvo = rightD.bottomLeftCorner(_nv, _no);
    Eigen::Ref<Eigen::MatrixXd> rightDvv = rightD.bottomRightCorner(_nv, _nv);

    Eigen::Ref<Eigen::MatrixXd> leftDoo = leftD.topLeftCorner(_no, _no);
    Eigen::Ref<Eigen::MatrixXd> leftDov = leftD.topRightCorner(_no, _nv);
    Eigen::Ref<Eigen::MatrixXd> leftDvo = leftD.bottomLeftCorner(_nv, _no);
    Eigen::Ref<Eigen::MatrixXd> leftDvv = leftD.bottomRightCorner(_nv, _nv);

    // reference eigenvectors and lagrangian multiplier
    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> excLagrange = eigenvectors[2].col(iEigen);
    double eigenvalue = eigenvalues(iEigen);

    // need to recalculate -Fia (stored in Fai) for each excited state
    _Fai = 2 * _Jia->operator*(_Jia->transpose() * rightEigenvector);
    Eigen::MatrixXd rev = Eigen::Map<Eigen::MatrixXd>(rightEigenvector.data(), _nv, _no);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
      _Fai -= jia * rev.transpose() * jia;
    }

    // remap to matrix objects
    Eigen::MatrixXd gsl = Eigen::Map<Eigen::MatrixXd>(_gsLagrange.data(), _nv, _no);
    Eigen::MatrixXd fia = Eigen::Map<Eigen::MatrixXd>(_Fia.data(), _nv, _no);

    Eigen::MatrixXd lev = Eigen::Map<Eigen::MatrixXd>(leftEigenvector.data(), _nv, _no);
    Eigen::MatrixXd exl = Eigen::Map<Eigen::MatrixXd>(excLagrange.data(), _nv, _no);
    Eigen::MatrixXd fai = Eigen::Map<Eigen::MatrixXd>(_Fai.data(), _nv, _no);

    // ij-wise right contributions [rightD]
    this->performTransformation(_Xia, rightEigenvector, false);
    this->performTransformation(_Zia, _gsLagrange, true);
    _trafoVector = _gsLagrange;
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j);
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue);
        Eigen::MatrixXd rrij = _soss * rij - _sss * rij.transpose();
        rightDvv.noalias() += tij * rij.transpose();
        rightDov.row(i).noalias() += rrij * gsl.col(j);
        if (iEigen == 0) {
          CC2virt.noalias() += ij * tij.transpose();
        }
        if (i != j) {
          rightDvv.noalias() += tij.transpose() * rij;
          rightDov.row(j).noalias() += rrij.transpose() * gsl.col(i);
          if (iEigen == 0) {
            CC2virt.noalias() += ij.transpose() * tij;
          }
        }
      }
    }
    rightDov.noalias() += rev.transpose();
    rightDov.noalias() -= (CC2virt * rev).transpose();
    rightDvv.noalias() += gsl * rev.transpose();
    /* ij-wise right contributions [rightD] */

    // here, transition-moment Lagrange multiplier are contracted, so this is only necessary for cc2
    if (_xwfModel == Options::LR_METHOD::CC2) {
      // ij-wise multiplier contributions [rightD]
      this->performTransformation(_Zia, excLagrange, true);
      _trafoVector = excLagrange;
      Eigen::MatrixXd rg = -1 * rev.transpose() * gsl;
      Eigen::MatrixXd gr = -1 * gsl * rev.transpose();
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
        Eigen::Map<Eigen::MatrixXd> xia(_Xia.data() + Q * _nv * _no, _nv, _no);
        xia = jia * rg + gr * jia;
      }
      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = i; j < _no; ++j) {
          Eigen::MatrixXd ij = this->getAmplitudesA(i, j);
          Eigen::MatrixXd iij = _soss * ij - _sss * ij.transpose();
          Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue);
          rightDvv.noalias() += mij * ij.transpose();
          rightDov.row(i).noalias() += iij * exl.col(j);
          if (i != j) {
            rightDvv.noalias() += mij.transpose() * ij;
            rightDov.row(j).noalias() += iij.transpose() * exl.col(i);
          }
        }
      }
      rightDvo.noalias() += exl;
      /* ij-wise multiplier contributions [rightD] */
    }

    // ij-wise left contributions [leftD]
    this->performTransformation(_Zia, leftEigenvector, true);
    _trafoVector = leftEigenvector;
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j);
        Eigen::MatrixXd iij = _soss * ij - _sss * ij.transpose();
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue);
        leftDvv.noalias() += lij * ij.transpose();
        leftDov.row(i).noalias() += iij * lev.col(j);
        if (i != j) {
          leftDvv.noalias() += lij.transpose() * ij;
          leftDov.row(j).noalias() += iij.transpose() * lev.col(i);
        }
      }
    }
    leftDvo.noalias() += lev;
    /* ij-wise left contributions [leftD] */

    // done with ij-wise amps, must reorder ints
    this->reorderTensor(*_Jia, _nv, _no, _nx);
    if (_Jia != _Bia) {
      this->reorderTensor(*_Bia, _nv, _no, _nx);
    }

    // ab-wise right contributions [rightD]
    this->performTransformation(_Xia, rightEigenvector, false);
    this->reorderTensor(_Xia, _nv, _no, _nx);
    this->performTransformation(_Zia, _gsLagrange, true);
    this->reorderTensor(_Zia, _nv, _no, _nx);
    this->reorderTensor(_Fia, _nv, _no);
    this->reorderTensor(_gsLagrange, _nv, _no);
    for (size_t a = 0; a < _nv; ++a) {
      for (size_t b = a; b < _nv; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b);
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue);
        rightDoo.noalias() -= rab * tab.transpose();
        if (iEigen == 0) {
          CC2occ.noalias() += ab * tab.transpose();
        }
        if (a != b) {
          rightDoo.noalias() -= rab.transpose() * tab;
          if (iEigen == 0) {
            CC2occ.noalias() += ab.transpose() * tab;
          }
        }
      }
    }
    rightDoo.noalias() -= rev.transpose() * gsl;
    rightDov.noalias() -= CC2occ * rev.transpose();
    /* ab-wise right contributions [rightD] */

    if (_xwfModel == Options::LR_METHOD::CC2) {
      // ab-wise multiplier contributions [rightD]
      Eigen::MatrixXd rg = -1 * rev * gsl.transpose();
      Eigen::MatrixXd gr = -1 * gsl.transpose() * rev;
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _no * _nv, _no, _nv);
        Eigen::Map<Eigen::MatrixXd> xia(_Xia.data() + Q * _no * _nv, _no, _nv);
        xia = gr * jia + jia * rg;
      }
      this->performTransformation(_Zia, excLagrange, true);
      _trafoVector = excLagrange;
      this->reorderTensor(_Zia, _nv, _no, _nx);
      this->reorderTensor(_Fai, _nv, _no);
      this->reorderTensor(_trafoVector, _nv, _no);
      for (size_t a = 0; a < _nv; ++a) {
        for (size_t b = a; b < _nv; ++b) {
          Eigen::MatrixXd ab = this->getAmplitudesAV(a, b);
          Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue);
          rightDoo.noalias() -= ab * mab.transpose();
          if (a != b) {
            rightDoo.noalias() -= ab.transpose() * mab;
          }
        }
      }
      /* ab-wise multiplier contributions [rightD] */
    }

    // ab-wise left contributions [leftD]
    this->performTransformation(_Zia, leftEigenvector, true);
    _trafoVector = leftEigenvector;
    this->reorderTensor(_Zia, _nv, _no, _nx);
    this->reorderTensor(_trafoVector, _nv, _no);
    for (size_t a = 0; a < _nv; ++a) {
      for (size_t b = a; b < _nv; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b);
        Eigen::MatrixXd lab = getLeftAmplitudesV(a, b, eigenvalue);
        leftDoo.noalias() -= ab * lab.transpose();
        if (a != b) {
          leftDoo.noalias() -= ab.transpose() * lab;
        }
      }
    }
    /* ab-wise left contributions [leftD] */

    // done with ab-wise amps, must order ints back
    this->reorderTensor(*_Jia, _no, _nv, _nx);
    if (_Jia != _Bia) {
      this->reorderTensor(*_Bia, _no, _nv, _nx);
    }

    // only these refer to the ground state
    this->reorderTensor(_Fia, _no, _nv);
    this->reorderTensor(_gsLagrange, _no, _nv);
  } /* iEigen */

  Timings::timeTaken("Exc. State WF - Density Matrices");
} /* this->calculateDensityMatrices() restricted */

template<>
void CC2Sigmavector<UNRESTRICTED>::calcDensityMatrices(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues,
                                                       std::vector<Eigen::MatrixXd>& densityMatrices) {
  Timings::takeTime("Exc. State WF - Density Matrices");

  unsigned nEigen = eigenvalues.size();
  unsigned alphaMO = _no_a + _nv_a;
  unsigned betaMO = _no_b + _nv_b;
  densityMatrices[0] = Eigen::MatrixXd::Zero(alphaMO * alphaMO + betaMO * betaMO, nEigen);
  densityMatrices[1] = Eigen::MatrixXd::Zero(alphaMO * alphaMO + betaMO * betaMO, nEigen);
  printf("  Calculating %3i eta and %3i xi density matrices.\n\n", nEigen, 2 * nEigen);

  // cc2 density matrix intermediate
  Eigen::MatrixXd CC2virt_a = Eigen::MatrixXd::Zero(_nv_a, _nv_a);
  Eigen::MatrixXd CC2occ_a = Eigen::MatrixXd::Zero(_no_a, _no_a);
  Eigen::MatrixXd CC2virt_b = Eigen::MatrixXd::Zero(_nv_b, _nv_b);
  Eigen::MatrixXd CC2occ_b = Eigen::MatrixXd::Zero(_no_b, _no_b);

  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    // will need two density matrix objects
    // for left and right transition moments
    //  rightD = D^eta(R) + D^xi(M)
    //   leftD = D^xi(L)
    Eigen::Map<Eigen::MatrixXd> rightD_a(densityMatrices[0].col(iEigen).data(), alphaMO, alphaMO);
    Eigen::Map<Eigen::MatrixXd> rightD_b(densityMatrices[0].col(iEigen).data() + alphaMO * alphaMO, betaMO, betaMO);
    Eigen::Map<Eigen::MatrixXd> leftD_a(densityMatrices[1].col(iEigen).data(), alphaMO, alphaMO);
    Eigen::Map<Eigen::MatrixXd> leftD_b(densityMatrices[1].col(iEigen).data() + alphaMO * alphaMO, betaMO, betaMO);

    // reference blocks
    Eigen::Ref<Eigen::MatrixXd> rightDoo_a = rightD_a.topLeftCorner(_no_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> rightDov_a = rightD_a.topRightCorner(_no_a, _nv_a);
    Eigen::Ref<Eigen::MatrixXd> rightDvo_a = rightD_a.bottomLeftCorner(_nv_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> rightDvv_a = rightD_a.bottomRightCorner(_nv_a, _nv_a);

    Eigen::Ref<Eigen::MatrixXd> leftDoo_a = leftD_a.topLeftCorner(_no_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> leftDov_a = leftD_a.topRightCorner(_no_a, _nv_a);
    Eigen::Ref<Eigen::MatrixXd> leftDvo_a = leftD_a.bottomLeftCorner(_nv_a, _no_a);
    Eigen::Ref<Eigen::MatrixXd> leftDvv_a = leftD_a.bottomRightCorner(_nv_a, _nv_a);

    Eigen::Ref<Eigen::MatrixXd> rightDoo_b = rightD_b.topLeftCorner(_no_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> rightDov_b = rightD_b.topRightCorner(_no_b, _nv_b);
    Eigen::Ref<Eigen::MatrixXd> rightDvo_b = rightD_b.bottomLeftCorner(_nv_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> rightDvv_b = rightD_b.bottomRightCorner(_nv_b, _nv_b);

    Eigen::Ref<Eigen::MatrixXd> leftDoo_b = leftD_b.topLeftCorner(_no_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> leftDov_b = leftD_b.topRightCorner(_no_b, _nv_b);
    Eigen::Ref<Eigen::MatrixXd> leftDvo_b = leftD_b.bottomLeftCorner(_nv_b, _no_b);
    Eigen::Ref<Eigen::MatrixXd> leftDvv_b = leftD_b.bottomRightCorner(_nv_b, _nv_b);

    // reference eigenvectors and lagrangian multiplier
    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> leftEigenvector = eigenvectors[1].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> excLagrange = eigenvectors[2].col(iEigen);
    double eigenvalue = eigenvalues(iEigen);

    // need to recalculate -Fia (stored in Fai) for each excited state
    Eigen::Map<Eigen::MatrixXd> mfai_a(_Fai.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> mfai_b(_Fai.data() + _alpha, _nv_b, _no_b);
    Eigen::VectorXd XQ =
        _Jia->alpha.transpose() * rightEigenvector.head(_alpha) + _Jia->beta.transpose() * rightEigenvector.tail(_beta);
    _Fai.head(_alpha) = _Jia->alpha * XQ;
    _Fai.tail(_beta) = _Jia->beta * XQ;
    Eigen::MatrixXd rev_a = Eigen::Map<Eigen::MatrixXd>(rightEigenvector.data(), _nv_a, _no_a);
    Eigen::MatrixXd rev_b = Eigen::Map<Eigen::MatrixXd>(rightEigenvector.data() + _alpha, _nv_b, _no_b);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
      mfai_a -= jia_a * rev_a.transpose() * jia_a;

      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
      mfai_b -= jia_b * rev_b.transpose() * jia_b;
    }

    // remap to matrix objects
    Eigen::MatrixXd gsl_a = Eigen::Map<Eigen::MatrixXd>(_gsLagrange.data(), _nv_a, _no_a);
    Eigen::MatrixXd fia_a = Eigen::Map<Eigen::MatrixXd>(_Fia.data(), _nv_a, _no_a);

    Eigen::MatrixXd lev_a = Eigen::Map<Eigen::MatrixXd>(leftEigenvector.data(), _nv_a, _no_a);
    Eigen::MatrixXd exl_a = Eigen::Map<Eigen::MatrixXd>(excLagrange.data(), _nv_a, _no_a);
    Eigen::MatrixXd fai_a = Eigen::Map<Eigen::MatrixXd>(_Fai.data(), _nv_a, _no_a);

    Eigen::MatrixXd gsl_b = Eigen::Map<Eigen::MatrixXd>(_gsLagrange.data() + _alpha, _nv_b, _no_b);
    Eigen::MatrixXd fia_b = Eigen::Map<Eigen::MatrixXd>(_Fia.data() + _alpha, _nv_b, _no_b);

    Eigen::MatrixXd lev_b = Eigen::Map<Eigen::MatrixXd>(leftEigenvector.data() + _alpha, _nv_b, _no_b);
    Eigen::MatrixXd exl_b = Eigen::Map<Eigen::MatrixXd>(excLagrange.data() + _alpha, _nv_b, _no_b);
    Eigen::MatrixXd fai_b = Eigen::Map<Eigen::MatrixXd>(_Fai.data() + _alpha, _nv_b, _no_b);

    // ij-wise right contributions [rightD]
    this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);
    this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 1);
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 1);
        Eigen::MatrixXd rrij = _sss * rij - _sss * rij.transpose();
        rightDvv_a.noalias() += tij * rij.transpose();
        rightDov_a.row(i).noalias() += rrij * gsl_a.col(j);
        if (iEigen == 0) {
          CC2virt_a.noalias() += ij * tij.transpose();
        }
        if (i != j) {
          rightDvv_a.noalias() += tij.transpose() * rij;
          rightDov_a.row(j).noalias() += rrij.transpose() * gsl_a.col(i);
          if (iEigen == 0) {
            CC2virt_a.noalias() += ij.transpose() * tij;
          }
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, -1);
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, -1);
        Eigen::MatrixXd rrij = _sss * rij - _sss * rij.transpose();
        rightDvv_b.noalias() += tij * rij.transpose();
        rightDov_b.row(i).noalias() += rrij * gsl_b.col(j);
        if (iEigen == 0) {
          CC2virt_b.noalias() += ij * tij.transpose();
        }
        if (i != j) {
          rightDvv_b.noalias() += tij.transpose() * rij;
          rightDov_b.row(j).noalias() += rrij.transpose() * gsl_b.col(i);
          if (iEigen == 0) {
            CC2virt_b.noalias() += ij.transpose() * tij;
          }
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 0);
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
        Eigen::MatrixXd rij = this->getRightAmplitudesA(i, j, eigenvalue, 0);
        rightDvv_a.noalias() += tij * rij.transpose();
        rightDov_a.row(i).noalias() += _oss * rij * gsl_b.col(j);
        rightDvv_b += tij.transpose() * rij;
        rightDov_b.row(j).noalias() += _oss * rij.transpose() * gsl_a.col(i);
        if (iEigen == 0) {
          CC2virt_a.noalias() += ij * tij.transpose();
          CC2virt_b.noalias() += ij.transpose() * tij;
        }
      }
    }
    rightDov_a.noalias() += rev_a.transpose();
    rightDov_a.noalias() -= (CC2virt_a * rev_a).transpose();
    rightDvv_a.noalias() += gsl_a * rev_a.transpose();
    rightDov_b.noalias() += rev_b.transpose();
    rightDov_b.noalias() -= (CC2virt_b * rev_b).transpose();
    rightDvv_b.noalias() += gsl_b * rev_b.transpose();
    /* ij-wise right contributions [rightD] */

    // here, transition-moment Lagrange multiplier are contracted, so this is only necessary for cc2
    if (_xwfModel == Options::LR_METHOD::CC2) {
      // ij-wise multiplier contributions [rightD]
      this->performTransformation(_Zia.alpha, excLagrange.head(_alpha), true, 1);
      this->performTransformation(_Zia.beta, excLagrange.tail(_beta), true, -1);
      _trafoVector = excLagrange;
      Eigen::MatrixXd rg_a = -1 * rev_a.transpose() * gsl_a;
      Eigen::MatrixXd gr_a = -1 * gsl_a * rev_a.transpose();
      Eigen::MatrixXd rg_b = -1 * rev_b.transpose() * gsl_b;
      Eigen::MatrixXd gr_b = -1 * gsl_b * rev_b.transpose();
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
        Eigen::Map<Eigen::MatrixXd> xia_a(_Xia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
        xia_a = jia_a * rg_a + gr_a * jia_a;

        Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
        Eigen::Map<Eigen::MatrixXd> xia_b(_Xia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
        xia_b = jia_b * rg_b + gr_b * jia_b;
      }
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 1);
          Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();
          Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue, 1);
          rightDvv_a.noalias() += mij * ij.transpose();
          rightDov_a.row(i).noalias() += iij * exl_a.col(j);
          if (i != j) {
            rightDvv_a.noalias() += mij.transpose() * ij;
            rightDov_a.row(j).noalias() += iij.transpose() * exl_a.col(i);
          }
        }
      }
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          Eigen::MatrixXd ij = this->getAmplitudesA(i, j, -1);
          Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();
          Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue, -1);
          rightDvv_b.noalias() += mij * ij.transpose();
          rightDov_b.row(i).noalias() += iij * exl_b.col(j);
          if (i != j) {
            rightDvv_b.noalias() += mij.transpose() * ij;
            rightDov_b.row(j).noalias() += iij.transpose() * exl_b.col(i);
          }
        }
      }
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 0);
          Eigen::MatrixXd mij = this->getELagrangeAmplitudes(i, j, eigenvalue, 0);
          rightDvv_a.noalias() += mij * ij.transpose();
          rightDov_a.row(i).noalias() += _oss * ij * exl_b.col(j);
          rightDvv_b.noalias() += mij.transpose() * ij;
          rightDov_b.row(j).noalias() += _oss * ij.transpose() * exl_a.col(i);
        }
      }
      rightDvo_a.noalias() += exl_a;
      rightDvo_b.noalias() += exl_b;
      /* ij-wise multiplier contributions [rightD] */
    }

    // ij-wise left contributions [leftD]
    this->performTransformation(_Zia.alpha, leftEigenvector.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, leftEigenvector.tail(_beta), true, -1);
    _trafoVector = leftEigenvector;
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 1);
        Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 1);
        leftDvv_a.noalias() += lij * ij.transpose();
        leftDov_a.row(i).noalias() += iij * lev_a.col(j);
        if (i != j) {
          leftDvv_a.noalias() += lij.transpose() * ij;
          leftDov_a.row(j).noalias() += iij.transpose() * lev_a.col(i);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, -1);
        Eigen::MatrixXd iij = _sss * ij - _sss * ij.transpose();
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, -1);
        leftDvv_b.noalias() += lij * ij.transpose();
        leftDov_b.row(i).noalias() += iij * lev_b.col(j);
        if (i != j) {
          leftDvv_b.noalias() += lij.transpose() * ij;
          leftDov_b.row(j).noalias() += iij.transpose() * lev_b.col(i);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd ij = this->getAmplitudesA(i, j, 0);
        Eigen::MatrixXd lij = this->getLeftAmplitudes(i, j, eigenvalue, 0);
        leftDvv_a.noalias() += lij * ij.transpose();
        leftDov_a.row(i).noalias() += _oss * ij * lev_b.col(j);
        leftDvv_b.noalias() += lij.transpose() * ij;
        leftDov_b.row(j).noalias() += _oss * ij.transpose() * lev_a.col(i);
      }
    }
    leftDvo_a.noalias() += lev_a;
    leftDvo_b.noalias() += lev_b;
    /* ij-wise left contributions [leftD] */

    // done with ij-wise amps, must reorder ints
    this->reorderTensor(_Jia->alpha, _nv_a, _no_a, _nx);
    this->reorderTensor(_Jia->beta, _nv_b, _no_b, _nx);
    if (_Jia != _Bia) {
      this->reorderTensor(_Bia->alpha, _nv_a, _no_a, _nx);
      this->reorderTensor(_Bia->beta, _nv_b, _no_b, _nx);
    }

    // ab-wise right contributions [rightD]
    this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);
    this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);
    this->reorderTensor(_Xia.alpha, _nv_a, _no_a, _nx);
    this->reorderTensor(_Xia.beta, _nv_b, _no_b, _nx);
    this->reorderTensor(_Zia.alpha, _nv_a, _no_a, _nx);
    this->reorderTensor(_Zia.beta, _nv_b, _no_b, _nx);

    this->reorderTensor(_Fia.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_Fia.tail(_beta), _nv_b, _no_b);
    this->reorderTensor(_gsLagrange.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_gsLagrange.tail(_beta), _nv_b, _no_b);

    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = a; b < _nv_a; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 1);
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, 1);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, 1);
        rightDoo_a.noalias() -= rab * tab.transpose();
        if (iEigen == 0) {
          CC2occ_a.noalias() += ab * tab.transpose();
        }
        if (a != b) {
          rightDoo_a.noalias() -= rab.transpose() * tab;
          if (iEigen == 0) {
            CC2occ_a.noalias() += ab.transpose() * tab;
          }
        }
      }
    }
    for (size_t a = 0; a < _nv_b; ++a) {
      for (size_t b = a; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, -1);
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, -1);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, -1);
        rightDoo_b.noalias() -= rab * tab.transpose();
        if (iEigen == 0) {
          CC2occ_b.noalias() += ab * tab.transpose();
        }
        if (a != b) {
          rightDoo_b.noalias() -= rab.transpose() * tab;
          if (iEigen == 0) {
            CC2occ_b.noalias() += ab.transpose() * tab;
          }
        }
      }
    }
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = 0; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 0);
        Eigen::MatrixXd tab = this->getGLagrangeAmplitudesV(a, b, 0);
        Eigen::MatrixXd rab = this->getRightAmplitudesAV(a, b, eigenvalue, 0);
        rightDoo_a.noalias() -= rab * tab.transpose();
        rightDoo_b.noalias() -= rab.transpose() * tab;
        if (iEigen == 0) {
          CC2occ_a.noalias() += ab * tab.transpose();
          CC2occ_b.noalias() += ab.transpose() * tab;
        }
      }
    }
    rightDoo_a.noalias() -= rev_a.transpose() * gsl_a;
    rightDov_a.noalias() -= CC2occ_a * rev_a.transpose();
    rightDoo_b.noalias() -= rev_b.transpose() * gsl_b;
    rightDov_b.noalias() -= CC2occ_b * rev_b.transpose();
    /* ab-wise right contributions [rightD] */

    if (_xwfModel == Options::LR_METHOD::CC2) {
      // ab-wise multiplier contributions [rightD]
      Eigen::MatrixXd rg_a = -1 * rev_a * gsl_a.transpose();
      Eigen::MatrixXd gr_a = -1 * gsl_a.transpose() * rev_a;
      Eigen::MatrixXd rg_b = -1 * rev_b * gsl_b.transpose();
      Eigen::MatrixXd gr_b = -1 * gsl_b.transpose() * rev_b;
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _no_a * _nv_a, _no_a, _nv_a);
        Eigen::Map<Eigen::MatrixXd> xia_a(_Xia.alpha.data() + Q * _no_a * _nv_a, _no_a, _nv_a);
        xia_a = gr_a * jia_a + jia_a * rg_a;

        Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _no_b * _nv_b, _no_b, _nv_b);
        Eigen::Map<Eigen::MatrixXd> xia_b(_Xia.beta.data() + Q * _no_b * _nv_b, _no_b, _nv_b);
        xia_b = gr_b * jia_b + jia_b * rg_b;
      }
      this->performTransformation(_Zia.alpha, excLagrange.head(_alpha), true, 1);
      this->performTransformation(_Zia.beta, excLagrange.tail(_beta), true, -1);
      _trafoVector = excLagrange;
      this->reorderTensor(_Zia.alpha, _nv_a, _no_a, _nx);
      this->reorderTensor(_Zia.beta, _nv_b, _no_b, _nx);
      this->reorderTensor(_Fai.head(_alpha), _nv_a, _no_a);
      this->reorderTensor(_Fai.tail(_beta), _nv_b, _no_b);
      this->reorderTensor(_trafoVector.head(_alpha), _nv_a, _no_a);
      this->reorderTensor(_trafoVector.tail(_beta), _nv_b, _no_b);
      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = a; b < _nv_a; ++b) {
          Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 1);
          Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue, 1);
          rightDoo_a.noalias() -= ab * mab.transpose();
          if (a != b) {
            rightDoo_a.noalias() -= ab.transpose() * mab;
          }
        }
      }
      for (size_t a = 0; a < _nv_b; ++a) {
        for (size_t b = a; b < _nv_b; ++b) {
          Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, -1);
          Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue, -1);
          rightDoo_b.noalias() -= ab * mab.transpose();
          if (a != b) {
            rightDoo_b.noalias() -= ab.transpose() * mab;
          }
        }
      }
      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = 0; b < _nv_b; ++b) {
          Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 0);
          Eigen::MatrixXd mab = this->getELagrangeAmplitudesV(a, b, eigenvalue, 0);
          rightDoo_a.noalias() -= ab * mab.transpose();
          rightDoo_b.noalias() -= ab.transpose() * mab;
        }
      }
      /* ab-wise multiplier contributions [rightD] */
    }

    // ab-wise left contributions [leftD]
    this->performTransformation(_Zia.alpha, leftEigenvector.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, leftEigenvector.tail(_beta), true, -1);
    _trafoVector = leftEigenvector;
    this->reorderTensor(_Zia.alpha, _nv_a, _no_a, _nx);
    this->reorderTensor(_Zia.beta, _nv_b, _no_b, _nx);
    this->reorderTensor(_trafoVector.head(_alpha), _nv_a, _no_a);
    this->reorderTensor(_trafoVector.tail(_beta), _nv_b, _no_b);
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = a; b < _nv_a; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 1);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, 1);
        leftDoo_a.noalias() -= ab * lab.transpose();
        if (a != b) {
          leftDoo_a.noalias() -= ab.transpose() * lab;
        }
      }
    }
    for (size_t a = 0; a < _nv_b; ++a) {
      for (size_t b = a; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, -1);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, -1);
        leftDoo_b.noalias() -= ab * lab.transpose();
        if (a != b) {
          leftDoo_b.noalias() -= ab.transpose() * lab;
        }
      }
    }
    for (size_t a = 0; a < _nv_a; ++a) {
      for (size_t b = 0; b < _nv_b; ++b) {
        Eigen::MatrixXd ab = this->getAmplitudesAV(a, b, 0);
        Eigen::MatrixXd lab = this->getLeftAmplitudesV(a, b, eigenvalue, 0);
        leftDoo_a.noalias() -= ab * lab.transpose();
        leftDoo_b.noalias() -= ab.transpose() * lab;
      }
    }
    /* ab-wise left contributions [leftD] */

    // done with ab-wise amps, must order ints back
    this->reorderTensor(_Jia->alpha, _no_a, _nv_a, _nx);
    this->reorderTensor(_Jia->beta, _no_b, _nv_b, _nx);
    if (_Jia != _Bia) {
      this->reorderTensor(_Bia->alpha, _no_a, _nv_a, _nx);
      this->reorderTensor(_Bia->beta, _no_b, _nv_b, _nx);
    }

    // only these refer to the ground state
    this->reorderTensor(_Fia.head(_alpha), _no_a, _nv_a);
    this->reorderTensor(_Fia.tail(_beta), _no_b, _nv_b);
    this->reorderTensor(_gsLagrange.head(_alpha), _no_a, _nv_a);
    this->reorderTensor(_gsLagrange.tail(_beta), _no_b, _nv_b);
  } /* iEigen */

  Timings::timeTaken("Exc. State WF - Density Matrices");
} /* this->calculateDensityMatrices() unrestricted */

template<>
void CC2Sigmavector<RESTRICTED>::calcTransitionMomentLagrangeMultiplier(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                        Eigen::VectorXd eigenvalues) {
  Timings::takeTime("Exc. State WF - EX Lagrange Mult");
  printf("\n    Setting up right-hand sides.\n");

  unsigned nEigen = eigenvectors[0].cols();

  // construct right-hand sides
  Eigen::MatrixXd rightHandSides = Eigen::MatrixXd::Zero(_no * _nv, nEigen);
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    // get eigenpair data
    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> rightHandSide = rightHandSides.col(iEigen);
    double eigenvalue = eigenvalues(iEigen);
    Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nv * _no);
    Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nv * _no);

    Eigen::Map<Eigen::MatrixXd> rev(rightEigenvector.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> rhs(rightHandSide.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> gsl(_gsLagrange.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> fia(_Fia.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> fai(_Fai.data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> tia(Tia.data(), _nv, _no);

    Eigen::MatrixXd Eij = Eigen::MatrixXd::Zero(_no, _no);
    Eigen::MatrixXd Eab = Eigen::MatrixXd::Zero(_nv, _nv);

    Eij += fia.transpose() * rev;
    Eab -= rev * fia.transpose();

    this->performTransformation(_Xia, rightEigenvector, false);

    _Zia.setZero();
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue);
        _Zia.middleRows(i * _nv, _nv).noalias() += rij * _Jia->middleRows(j * _nv, _nv);
        Tia.segment(i * _nv, _nv).noalias() += rij * _gsLagrange.segment(j * _nv, _nv);
        if (i != j) {
          _Zia.middleRows(j * _nv, _nv).noalias() += rij.transpose() * _Jia->middleRows(i * _nv, _nv);
          Tia.segment(j * _nv, _nv).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv, _nv);
        }
      }
    }

    Eigen::VectorXd XQ = _Jia->transpose() * rightEigenvector;
    _Fai = 2 * _Jia->operator*(XQ);

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
      Eigen::Map<Eigen::MatrixXd> zia(_Zia.data() + Q * _nv * _no, _nv, _no);
      // ^Fij and ^Fab
      Eij.noalias() += 2 * bij * XQ(Q);

      Eij.noalias() -= jia.transpose() * rev * bij;

      Eij.noalias() += jia.transpose() * zia;
      Eab.noalias() -= zia * jia.transpose();

      fai.noalias() -= jia * rev.transpose() * jia;
    }

    Eab += this->transformedABFock(rightEigenvector, XQ, _Xia);

    // 0 contribution
    rhs += Eab.transpose() * gsl - gsl * Eij.transpose();

    // I contribution
    // this is nowhere mentioned in turbomole-related literature
    Tia += rightEigenvector;

    rightHandSide += 2 * _Jia->operator*(_Jia->transpose() * Tia);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _no * _nv, _nv, _no);
      rhs -= jia * tia.transpose() * jia;
    }

    // GH contribution
    Eigen::MatrixXd rg = -1 * rev.transpose() * gsl;
    Eigen::MatrixXd gr = -1 * gsl * rev.transpose();
    _Xia.resize(_nv * _no, _nx);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
      Eigen::Map<Eigen::MatrixXd> xia(_Xia.data() + Q * _nv * _no, _nv, _no);
      xia = jia * rg + gr * jia;
    }

    // fock amplitude contribution
    _Yia = Eigen::MatrixXd::Zero(_nv * _no, _nx);
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue);
        _Yia.middleRows(i * _nv, _nv).noalias() += fij * _Bia->middleRows(j * _nv, _nv);
        if (i != j) {
          _Yia.middleRows(j * _nv, _nv).noalias() += fij.transpose() * _Bia->middleRows(i * _nv, _nv);
        }
      }
    }
    /* fock amplitude contribution */

    // ground state Lagrangian doubles contribution
    this->performTransformation(_Xia, rightEigenvector, false);
    this->performTransformation(_Zia, _gsLagrange, true);
    _trafoVector = _gsLagrange;

    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
        _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Xia.middleRows(j * _nv, _nv);
        if (i != j) {
          _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Xia.middleRows(i * _nv, _nv);
        }
      }
    }

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
      Eigen::Map<Eigen::MatrixXd> yia(_Yia.data() + Q * _no * _nv, _nv, _no);
      rhs -= yia * bij.transpose();
    }
    rightHandSide += this->getJ2GContribution(_Yia.data(), nullptr, nullVector, true);

    _Yia.setZero();
    for (size_t i = 0; i < _no; ++i) {
      for (size_t j = i; j < _no; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j);
        _Yia.middleRows(i * _nv, _nv).noalias() += tij * _Bia->middleRows(j * _nv, _nv);
        if (i != j) {
          _Yia.middleRows(j * _nv, _nv).noalias() += tij.transpose() * _Bia->middleRows(i * _nv, _nv);
        }
      }
    }
    rightHandSide += 2 * _Jia->operator*(_Xia.transpose() * _gsLagrange);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
      Eigen::Map<Eigen::MatrixXd> jia(_Jia->data() + Q * _nv * _no, _nv, _no);
      Eigen::Map<Eigen::MatrixXd> yia(_Yia.data() + Q * _nv * _no, _nv, _no);

      // right transformed Bab/Bij
      Eigen::MatrixXd _bab = -1 * rev * jia.transpose();
      Eigen::MatrixXd _bij = jia.transpose() * rev;

      rhs.noalias() += _bab.transpose() * yia;
      rhs.noalias() -= yia * _bij.transpose();
      rhs.noalias() -= _bab.transpose() * gsl * bij.transpose();

      yia.noalias() = gsl * _bij.transpose();
    }
    rightHandSide -= this->getJ2GContribution(_Yia.data(), nullptr, nullVector, true);
    /* ground state Lagrangian doubles contribution */
  }
  rightHandSides *= -1;
  _Yia.resize(0, 0);

  eigenvectors[2] = this->lagrangianSubspaceSolve(rightHandSides, eigenvalues);

  printf("\n  Norm of lagrange-multiplier singles (0th root is the ground state):\n\n");
  printTableHead(" root     norm ");
  printf("%6i %11.6f\n", 0, _gsLagrange.norm());
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    printf("%6i %11.6f\n", iEigen + 1, eigenvectors[2].colwise().norm()(iEigen));
  }
  printf("\n");
  Timings::timeTaken("Exc. State WF - EX Lagrange Mult");
} /* this->calcTransitionMomentLagrangeMultiplier() restricted */

template<>
void CC2Sigmavector<UNRESTRICTED>::calcTransitionMomentLagrangeMultiplier(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                                          Eigen::VectorXd eigenvalues) {
  Timings::takeTime("Exc. State WF - EX Lagrange Mult");
  printf("\n    Setting up right-hand sides.\n");

  unsigned nEigen = eigenvectors[0].cols();

  // construct right-hand sides
  Eigen::MatrixXd rightHandSides = Eigen::MatrixXd::Zero(_nDim, nEigen);
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    // get eigenpair data
    Eigen::Ref<Eigen::VectorXd> rightEigenvector = eigenvectors[0].col(iEigen);
    Eigen::Ref<Eigen::VectorXd> rightHandSide = rightHandSides.col(iEigen);
    double eigenvalue = eigenvalues(iEigen);
    Eigen::VectorXd Tia = Eigen::VectorXd::Zero(_nDim);
    Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nDim);
    _Fai = Eigen::VectorXd::Zero(_nDim);

    Eigen::Map<Eigen::MatrixXd> rev_a(rightEigenvector.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> rhs_a(rightHandSide.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> gsl_a(_gsLagrange.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> fia_a(_Fia.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> fai_a(_Fai.data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> tia_a(Tia.data(), _nv_a, _no_a);

    Eigen::Map<Eigen::MatrixXd> rev_b(rightEigenvector.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> rhs_b(rightHandSide.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> gsl_b(_gsLagrange.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> fia_b(_Fia.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> fai_b(_Fai.data() + _alpha, _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> tia_b(Tia.data() + _alpha, _nv_b, _no_b);

    Eigen::MatrixXd Eij_a = Eigen::MatrixXd::Zero(_no_a, _no_a);
    Eigen::MatrixXd Eab_a = Eigen::MatrixXd::Zero(_nv_a, _nv_a);
    Eigen::MatrixXd Eij_b = Eigen::MatrixXd::Zero(_no_b, _no_b);
    Eigen::MatrixXd Eab_b = Eigen::MatrixXd::Zero(_nv_b, _nv_b);

    Eij_a += fia_a.transpose() * rev_a;
    Eab_a -= rev_a * fia_a.transpose();
    Eij_b += fia_b.transpose() * rev_b;
    Eab_b -= rev_b * fia_b.transpose();

    this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);

    _Zia.alpha.setZero();
    _Zia.beta.setZero();
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 1);
        _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
        Tia.segment(i * _nv_a, _nv_a).noalias() += rij * _gsLagrange.segment(j * _nv_a, _nv_a);
        if (i != j) {
          _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
          Tia.segment(j * _nv_a, _nv_a).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv_a, _nv_a);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, -1);
        _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
        Tia.segment(i * _nv_b + _alpha, _nv_b).noalias() += rij * _gsLagrange.segment(j * _nv_b + _alpha, _nv_b);
        if (i != j) {
          _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
          Tia.segment(j * _nv_b + _alpha, _nv_b).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv_b + _alpha, _nv_b);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd rij = this->getRightAmplitudes(i, j, eigenvalue, 0);
        _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += rij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
        Tia.segment(i * _nv_a, _nv_a).noalias() += rij * _gsLagrange.segment(j * _nv_b + _alpha, _nv_b);
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += rij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
        Tia.segment(j * _nv_b + _alpha, _nv_b).noalias() += rij.transpose() * _gsLagrange.segment(i * _nv_a, _nv_a);
      }
    }

    Eigen::VectorXd XQ =
        _Jia->alpha.transpose() * rightEigenvector.head(_alpha) + _Jia->beta.transpose() * rightEigenvector.tail(_beta);
    _Fai.head(_alpha) = _Jia->alpha * XQ;
    _Fai.tail(_beta) = _Jia->beta * XQ;

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
      // ^Fij and ^Fab
      Eij_a.noalias() += bij_a * XQ(Q);

      Eij_a.noalias() -= jia_a.transpose() * rev_a * bij_a;

      Eij_a.noalias() += jia_a.transpose() * zia_a;
      Eab_a.noalias() -= zia_a * jia_a.transpose();

      fai_a.noalias() -= jia_a * rev_a.transpose() * jia_a;

      Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
      // ^Fij and ^Fab
      Eij_b.noalias() += bij_b * XQ(Q);

      Eij_b.noalias() -= jia_b.transpose() * rev_b * bij_b;

      Eij_b.noalias() += jia_b.transpose() * zia_b;
      Eab_b.noalias() -= zia_b * jia_b.transpose();

      fai_b.noalias() -= jia_b * rev_b.transpose() * jia_b;
    }

    Eab_a += this->transformedABFock(rightEigenvector.head(_alpha), XQ, _Xia.alpha, 1);
    Eab_b += this->transformedABFock(rightEigenvector.tail(_beta), XQ, _Xia.beta, -1);

    // 0 contribution
    rhs_a += Eab_a.transpose() * gsl_a - gsl_a * Eij_a.transpose();
    rhs_b += Eab_b.transpose() * gsl_b - gsl_b * Eij_b.transpose();

    // I contribution
    // this is nowhere mentioned in turbomole-related literature
    Tia += rightEigenvector;

    XQ = _Jia->alpha.transpose() * Tia.head(_alpha) + _Jia->beta.transpose() * Tia.tail(_beta);
    rightHandSide.head(_alpha) += _Jia->alpha * XQ;
    rightHandSide.tail(_beta) += _Jia->beta * XQ;
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _no_a * _nv_a, _nv_a, _no_a);
      rhs_a -= jia_a * tia_a.transpose() * jia_a;

      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _no_b * _nv_b, _nv_b, _no_b);
      rhs_b -= jia_b * tia_b.transpose() * jia_b;
    }

    // GH contribution
    Eigen::MatrixXd rg_a = -1 * rev_a.transpose() * gsl_a;
    Eigen::MatrixXd gr_a = -1 * gsl_a * rev_a.transpose();
    Eigen::MatrixXd rg_b = -1 * rev_b.transpose() * gsl_b;
    Eigen::MatrixXd gr_b = -1 * gsl_b * rev_b.transpose();
    _Xia.alpha.resize(_nv_a * _no_a, _nx);
    _Xia.beta.resize(_nv_b * _no_b, _nx);
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> xia_a(_Xia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
      xia_a = jia_a * rg_a + gr_a * jia_a;

      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> xia_b(_Xia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
      xia_b = jia_b * rg_b + gr_b * jia_b;
    }

    // fock amplitude contribution
    _Yia.alpha = Eigen::MatrixXd::Zero(_nv_a * _no_a, _nx);
    _Yia.beta = Eigen::MatrixXd::Zero(_nv_b * _no_b, _nx);
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue, 1);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += fij * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += fij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue, -1);
        _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += fij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += fij.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd fij = this->getFockAmplitudes(i, j, eigenvalue, 0);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += fij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += fij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
    /* fock amplitude contribution */

    // ground state Lagrangian doubles contribution
    this->performTransformation(_Xia.alpha, rightEigenvector.head(_alpha), false, 1);
    this->performTransformation(_Xia.beta, rightEigenvector.tail(_beta), false, -1);
    this->performTransformation(_Zia.alpha, _gsLagrange.head(_alpha), true, 1);
    this->performTransformation(_Zia.beta, _gsLagrange.tail(_beta), true, -1);
    _trafoVector = _gsLagrange;

    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Xia.alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Xia.alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
        _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Xia.beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Xia.beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Xia.beta.middleRows(j * _nv_b, _nv_b);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Xia.alpha.middleRows(i * _nv_a, _nv_a);
      }
    }

    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.data() + Q * _no_a * _nv_a, _nv_a, _no_a);
      rhs_a -= yia_a * bij_a.transpose();

      Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.data() + Q * _no_b * _nv_b, _nv_b, _no_b);
      rhs_b -= yia_b * bij_b.transpose();
    }
    rightHandSide.head(_alpha) += this->getJ2GContribution(_Yia.alpha.data(), nullptr, nullVector, true, 1);
    rightHandSide.tail(_beta) += this->getJ2GContribution(_Yia.beta.data(), nullptr, nullVector, true, -1);

    _Yia.alpha = Eigen::MatrixXd::Zero(_nv_a * _no_a, _nx);
    _Yia.beta = Eigen::MatrixXd::Zero(_nv_b * _no_b, _nx);
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = i; j < _no_a; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 1);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
        if (i != j) {
          _Yia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += tij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
        }
      }
    }
    for (size_t i = 0; i < _no_b; ++i) {
      for (size_t j = i; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, -1);
        _Yia.beta.middleRows(i * _nv_b, _nv_b).noalias() += tij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        if (i != j) {
          _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
        }
      }
    }
    for (size_t i = 0; i < _no_a; ++i) {
      for (size_t j = 0; j < _no_b; ++j) {
        Eigen::MatrixXd tij = this->getGLagrangeAmplitudes(i, j, 0);
        _Yia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += tij * _Bia->beta.middleRows(j * _nv_b, _nv_b);
        _Yia.beta.middleRows(j * _nv_b, _nv_b).noalias() += tij.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }

    XQ = _Xia.alpha.transpose() * _gsLagrange.head(_alpha) + _Xia.beta.transpose() * _gsLagrange.tail(_beta);
    rightHandSide.head(_alpha) += _Jia->alpha * XQ;
    rightHandSide.tail(_beta) += _Jia->beta * XQ;
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> yia_a(_Yia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);

      // right transformed Bab/Bij
      Eigen::MatrixXd _bab_a = -1 * rev_a * jia_a.transpose();
      Eigen::MatrixXd _bij_a = jia_a.transpose() * rev_a;

      rhs_a.noalias() += _bab_a.transpose() * yia_a;
      rhs_a.noalias() -= yia_a * _bij_a.transpose();
      rhs_a.noalias() -= _bab_a.transpose() * gsl_a * bij_a.transpose();

      yia_a.noalias() = gsl_a * _bij_a.transpose();

      Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
      Eigen::Map<Eigen::MatrixXd> yia_b(_Yia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);

      // right transformed Bab/Bij
      Eigen::MatrixXd _bab_b = -1 * rev_b * jia_b.transpose();
      Eigen::MatrixXd _bij_b = jia_b.transpose() * rev_b;

      rhs_b.noalias() += _bab_b.transpose() * yia_b;
      rhs_b.noalias() -= yia_b * _bij_b.transpose();
      rhs_b.noalias() -= _bab_b.transpose() * gsl_b * bij_b.transpose();

      yia_b.noalias() = gsl_b * _bij_b.transpose();
    }
    rightHandSide.head(_alpha) -= this->getJ2GContribution(_Yia.alpha.data(), nullptr, nullVector, true, 1);
    rightHandSide.tail(_beta) -= this->getJ2GContribution(_Yia.beta.data(), nullptr, nullVector, true, -1);
    // /* ground state Lagrangian doubles contribution */
  }
  rightHandSides *= -1;
  _Yia.alpha.resize(0, 0);
  _Yia.beta.resize(0, 0);

  eigenvectors[2] = this->lagrangianSubspaceSolve(rightHandSides, eigenvalues);

  printf("\n  Norm of lagrange-multiplier singles (0th root is the ground state):\n\n");
  printTableHead(" root     norm ");
  printf("%6i %11.6f\n", 0, _gsLagrange.norm());
  for (unsigned iEigen = 0; iEigen < nEigen; ++iEigen) {
    printf("%6i %11.6f\n", iEigen + 1, eigenvectors[2].colwise().norm()(iEigen));
  }
  printf("\n");
  Timings::timeTaken("Exc. State WF - EX Lagrange Mult");
} /* this->calcTransitionMomentLagrangeMultiplier() unrestricted */

template<>
void CC2Sigmavector<RESTRICTED>::calcGroundStateLagrangeMultiplier() {
  Timings::takeTime("Exc. State WF - GS Lagrange Mult");
  printf("\n    Setting up right-hand sides.\n");

  _Zia.setZero();
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      Eigen::MatrixXd iajb = _Jia->middleRows(i * _nv, _nv) * _Jia->middleRows(j * _nv, _nv).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv, _nv, _e(i) + _e(j)) - _eVirt;
      Eigen::MatrixXd amps = (_soss * iajb - _sss * iajb.transpose()).cwiseQuotient(denom);
      _Zia.middleRows(i * _nv, _nv).noalias() += amps * _Bia->middleRows(j * _nv, _nv);
      if (i != j) {
        _Zia.middleRows(j * _nv, _nv).noalias() += amps.transpose() * _Bia->middleRows(i * _nv, _nv);
      }
    }
  }

  Eigen::MatrixXd rightHandSide = _Fia;
  Eigen::Map<Eigen::MatrixXd> rhs(rightHandSide.data(), _nv, _no);
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->data() + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> zia(_Zia.data() + Q * _nv * _no, _nv, _no);
    rhs -= zia * bij.transpose();
  }
  Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nv * _no);
  rightHandSide.col(0) += this->getJ2GContribution(_Zia.data(), nullptr, nullVector, true);
  rightHandSide *= -1;

  Eigen::VectorXd zeroEigenvalue = Eigen::VectorXd::Zero(1);
  _gsLagrange = this->lagrangianSubspaceSolve(rightHandSide, zeroEigenvalue);
  Timings::timeTaken("Exc. State WF - GS Lagrange Mult");
} /* this->calcGroundStateLagrangeMultiplier() restricted */

template<>
void CC2Sigmavector<UNRESTRICTED>::calcGroundStateLagrangeMultiplier() {
  Timings::takeTime("Exc. State WF - GS Lagrange Mult");
  printf("\n    Setting up right-hand sides.\n");

  _Zia.alpha.setZero();
  _Zia.beta.setZero();
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      Eigen::MatrixXd iajb = _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Jia->alpha.middleRows(j * _nv_a, _nv_a).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_a, _e.alpha(i) + _e.alpha(j)) - _eVirt.alpha;
      Eigen::MatrixXd amps = (_sss * iajb - _sss * iajb.transpose()).cwiseQuotient(denom);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += amps * _Bia->alpha.middleRows(j * _nv_a, _nv_a);
      if (i != j) {
        _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += amps.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
      }
    }
  }
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      Eigen::MatrixXd iajb = _Jia->beta.middleRows(i * _nv_b, _nv_b) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_b, _nv_b, _e.beta(i) + _e.beta(j)) - _eVirt.beta;
      Eigen::MatrixXd amps = (_sss * iajb - _sss * iajb.transpose()).cwiseQuotient(denom);
      _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += amps * _Bia->beta.middleRows(j * _nv_b, _nv_b);
      if (i != j) {
        _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += amps.transpose() * _Bia->beta.middleRows(i * _nv_b, _nv_b);
      }
    }
  }
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      Eigen::MatrixXd iajb = _Jia->alpha.middleRows(i * _nv_a, _nv_a) * _Jia->beta.middleRows(j * _nv_b, _nv_b).transpose();
      Eigen::MatrixXd denom = Eigen::MatrixXd::Constant(_nv_a, _nv_b, _e.alpha(i) + _e.beta(j)) - _eVirtab;
      Eigen::MatrixXd amps = (_oss * iajb).cwiseQuotient(denom);
      _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += amps * _Bia->beta.middleRows(j * _nv_b, _nv_b);
      _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += amps.transpose() * _Bia->alpha.middleRows(i * _nv_a, _nv_a);
    }
  }

  Eigen::MatrixXd rightHandSide = _Fia;
  Eigen::Map<Eigen::MatrixXd> rhs_a(rightHandSide.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> rhs_b(rightHandSide.data() + _alpha, _nv_b, _no_b);
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.data() + Q * _no_a * _no_a, _no_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.data() + Q * _nv_a * _no_a, _nv_a, _no_a);
    rhs_a -= zia_a * bij_a.transpose();

    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.data() + Q * _no_b * _no_b, _no_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.data() + Q * _nv_b * _no_b, _nv_b, _no_b);
    rhs_b -= zia_b * bij_b.transpose();
  }
  Eigen::VectorXd nullVector = Eigen::VectorXd::Zero(_nDim);
  rightHandSide.col(0).head(_alpha) += this->getJ2GContribution(_Zia.alpha.data(), nullptr, nullVector, true, 1);
  rightHandSide.col(0).tail(_beta) += this->getJ2GContribution(_Zia.beta.data(), nullptr, nullVector, true, -1);
  rightHandSide *= -1;

  Eigen::VectorXd zeroEigenvalue = Eigen::VectorXd::Zero(1);
  _gsLagrange = this->lagrangianSubspaceSolve(rightHandSide, zeroEigenvalue);
  Timings::timeTaken("Exc. State WF - GS Lagrange Mult");
} /* this->calcGroundStateLagrangeMultiplier() unrestricted */

template class CC2Sigmavector<Options::SCF_MODES::RESTRICTED>;
template class CC2Sigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity