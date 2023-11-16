/**
 * @file CC2IntDirectFunctions.cpp
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

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void CC2Controller<SCFMode>::reorderTensor(Eigen::Ref<Eigen::MatrixXd> Yia, unsigned n, unsigned m) {
  Timings::takeTime("CC2 -          Tensor Reorder");
  for (long j = 0; j < Yia.cols(); ++j) {
    Eigen::Map<Eigen::MatrixXd> yian(Yia.col(j).data(), n, m);
    Eigen::Map<Eigen::MatrixXd> yiat(Yia.col(j).data(), m, n);
    yiat = yian.transpose().eval();
  }
  Timings::timeTaken("CC2 -          Tensor Reorder");
} /* this->reorderTensor() */

template<Options::SCF_MODES SCFMode>
void CC2Controller<SCFMode>::performRITransformation(Eigen::MatrixXd& Yia, unsigned no, bool transposeRITrafo) {
  Timings::takeTime("CC2 -            RI Transform");
  unsigned blockSize = Yia.rows() / no;
  bool naf = _nx < _nxb;

  if (transposeRITrafo && naf) {
    Yia.conservativeResize(Yia.rows(), _nxb);
  }

  for (size_t i = 0; i < no; ++i) {
    if (transposeRITrafo) {
      Yia.block(i * blockSize, 0, blockSize, _nxb) =
          (Yia.block(i * blockSize, 0, blockSize, _nx) * _riTrafo.transpose()).eval();
    }
    else {
      Yia.block(i * blockSize, 0, blockSize, _nx) = (Yia.block(i * blockSize, 0, blockSize, _nxb) * _riTrafo).eval();
    }
  }

  if (!transposeRITrafo && naf) {
    Yia.conservativeResize(Yia.rows(), _nx);
  }
  Timings::timeTaken("CC2 -            RI Transform");
} /* this->performRITransformation() restricted/unrestricted */

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd CC2Controller<SCFMode>::getJ2GContribution(double* YiaPtr, double* JijPtr,
                                                           Eigen::Ref<Eigen::VectorXd> guessVector, bool fromLeft, int iSpin) {
  Timings::takeTime("CC2 -            J2G Contrib.");
  // declarations
  unsigned _nv = (iSpin > 0) ? _nv_a : _nv_b;
  unsigned _no = (iSpin > 0) ? _no_a : _no_b;

  Eigen::MatrixXd coeffA;
  Eigen::MatrixXd coeffB;
  int fSpin = 1;
  for_spin(_P, _H) {
    if (iSpin * fSpin > 0) {
      coeffA = fromLeft ? _H_spin.middleCols(_no, _nv) : _P_spin.middleCols(_no, _nv);
      coeffB = fromLeft ? _P_spin.middleCols(_no, _nv) : _H_spin.middleCols(_no, _nv);
    }
    fSpin *= -1;
  };

  Eigen::VectorXd sigmaVector = Eigen::VectorXd::Zero(_nv * _no);
  Eigen::Map<Eigen::MatrixXd> sigma(sigmaVector.data(), _nv, _no);
  Eigen::Map<Eigen::MatrixXd> guess(guessVector.data(), _nv, _no);

  unsigned nThreads = omp_get_max_threads();
  Eigen::MatrixXd Gia = Eigen::MatrixXd::Zero(_nb * _no, _nxb);
  Eigen::MatrixXd Eta = Eigen::MatrixXd::Zero(_nb * _no, nThreads);

  // contractions
  if (YiaPtr) {
    for (size_t Q = 0; Q < _nx; ++Q) {
      Eigen::Map<Eigen::MatrixXd> gia(Gia.col(Q).data(), _nb, _no);
      Eigen::Map<Eigen::MatrixXd> yia(YiaPtr + Q * _nv * _no, _nv, _no);
      gia.noalias() += coeffB * yia;
    }
  }
  if (JijPtr) {
    if (fromLeft) {
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> gia(Gia.col(Q).data(), _nb, _no);
        Eigen::Map<Eigen::MatrixXd> jij(JijPtr + Q * _no * _no, _no, _no);
        gia.noalias() -= coeffB * (guess * jij.transpose());
      }
    }
    else {
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> gia(Gia.col(Q).data(), _nb, _no);
        Eigen::Map<Eigen::MatrixXd> jij(JijPtr + Q * _no * _no, _no, _no);
        gia.noalias() -= coeffB * guess * jij;
      }
    }
  }

  Timings::timeTaken("CC2 -            J2G Contrib.");
  this->performRITransformation(Gia, _no, true);
  Timings::takeTime("CC2 -            J2G Contrib.");

  // loop over the rest
  auto distribute = [&](Eigen::Map<Eigen::MatrixXd> AO, size_t P, unsigned iThread) {
    Eigen::Map<Eigen::MatrixXd> eta(Eta.col(iThread).data(), _nb, _no);
    Eigen::Map<Eigen::MatrixXd> gia(Gia.col(P).data(), _nb, _no);
    eta.noalias() += AO * gia;
  };

  _integrals->loop(distribute);

  // sum over threads
  for (size_t iThread = 0; iThread < nThreads; ++iThread) {
    Eigen::Map<Eigen::MatrixXd> eta(Eta.col(iThread).data(), _nb, _no);
    sigma.noalias() += coeffA.transpose() * eta;
  }

  Timings::timeTaken("CC2 -            J2G Contrib.");
  return sigmaVector;
} /* this->getJ2GContribution() restricted/unrestricted */

template<Options::SCF_MODES SCFMode>
void CC2Controller<SCFMode>::performTransformation(Eigen::MatrixXd& Xia, Eigen::Ref<Eigen::VectorXd> trafoVector,
                                                   bool fromLeft, int iSpin, bool singlesTransformed) {
  Timings::takeTime("CC2 -           Int-Transform");

  // declarations
  unsigned _nv = (iSpin > 0) ? _nv_a : _nv_b;
  unsigned _no = (iSpin > 0) ? _no_a : _no_b;

  double* JijPtr;
  double* coeffAPtr;
  double* coeffBPtr;

  int fSpin = 1;
  auto& Jij = *_Jij;
  auto& Bij = *_Bij;
  for_spin(_C, _P, _H, Jij, Bij) {
    if (iSpin * fSpin > 0) {
      if (singlesTransformed) {
        coeffAPtr = fromLeft ? _H_spin.data() : _P_spin.data();
        coeffBPtr = fromLeft ? _P_spin.data() : _H_spin.data();
        JijPtr = Bij_spin.data();
      }
      else {
        coeffAPtr = _C_spin.data();
        coeffBPtr = _C_spin.data();
        JijPtr = Jij_spin.data();
      }
    }
    fSpin *= -1;
  };

  Eigen::Map<Eigen::MatrixXd> coeffA(coeffAPtr + _nb * _no, _nb, _nv);
  Eigen::Map<Eigen::MatrixXd> coeffB(coeffBPtr + _nb * _no, _nb, _nv);

  Eigen::Map<Eigen::MatrixXd> trafo(trafoVector.data(), _nv, _no);
  Eigen::MatrixXd eta = coeffB * trafo;

  Xia.resize(_nb * _no, _nxb);

  auto distribute = [&](Eigen::Map<Eigen::MatrixXd> AO, unsigned long P, unsigned) {
    Eigen::Map<Eigen::MatrixXd> xia(Xia.col(P).data(), _nb, _no);
    xia.noalias() = AO * eta;
  };

  _integrals->loop(distribute);
  this->performRITransformation(Xia, _no);

  // Jij contribution
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jij(JijPtr + Q * _no * _no, _no, _no);
    Eigen::Map<Eigen::MatrixXd> xia_b(Xia.col(Q).data(), _nb, _no);
    Eigen::Map<Eigen::MatrixXd> xia_v(Xia.col(Q).data(), _nv, _no);
    xia_v = coeffA.rightCols(_nv).transpose() * xia_b;

    if (fromLeft) {
      xia_v.noalias() -= trafo * jij.transpose();
    }
    else {
      xia_v.noalias() -= trafo * jij;
    }
  }

  Xia.conservativeResize(_nv * _no, _nx);
  Timings::timeTaken("CC2 -           Int-Transform");
} /* this->performTransformation() restricted/unrestricted */

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd CC2Controller<SCFMode>::transformedABFock(Eigen::Ref<Eigen::VectorXd> trafovector,
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
  Eigen::MatrixXd eta = coeffB * trafo;
  Yia = Eigen::MatrixXd::Zero(this->_nb * _no, this->_nxb);

  Eigen::VectorXd XP = this->_riTrafo * XQ;
  unsigned nThreads = omp_get_max_threads();
  Eigen::MatrixXd pseudoFock = Eigen::MatrixXd::Zero(this->_nb * this->_nb, nThreads);

  auto distribute = [&](Eigen::Map<Eigen::MatrixXd> AO, size_t P, unsigned iThread) {
    Eigen::Map<Eigen::MatrixXd> psf(pseudoFock.col(iThread).data(), this->_nb, this->_nb);
    Eigen::Map<Eigen::MatrixXd> yia(Yia.col(P).data(), this->_nb, _no);
    psf.noalias() += AO * XP(P);
    yia.noalias() = AO * eta;
  };
  _integrals->loop(distribute);

  this->performRITransformation(Yia, _no);

  unsigned spinFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;

  // sum over all threads
  for (size_t iThread = 0; iThread < nThreads; ++iThread) {
    Eigen::Map<Eigen::MatrixXd> psf(pseudoFock.col(iThread).data(), this->_nb, this->_nb);
    Fab.noalias() += spinFactor * coeffA.transpose() * psf * coeffB;
  }

  for (size_t Q = 0; Q < this->_nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(JiaPtr + Q * _nv * _no, _nv, _no);
    Eigen::Map<Eigen::MatrixXd> yia(Yia.col(Q).data(), this->_nb, _no);
    Fab.noalias() -= coeffA.transpose() * yia * jia.transpose();
  }
  Yia.resize(_nv * _no, _nx);

  return Fab;
} /* this->transformedABFock() restricted / unrestricted */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity