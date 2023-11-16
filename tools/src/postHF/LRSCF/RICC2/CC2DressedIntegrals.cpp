/**
 * @file CC2DressedIntegrals.cpp
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

namespace Serenity {

template<>
void CC2Controller<RESTRICTED>::transformIntegrals() {
  unsigned mo = _no + _nv;
  Eigen::MatrixXd t1 = Eigen::MatrixXd::Zero(mo, mo);
  t1.bottomLeftCorner(_nv, _no) = Eigen::Map<Eigen::MatrixXd>(_singles.data(), _nv, _no);
  _P = _C;
  _H = _C;
  _P.leftCols(mo) = _C.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) - t1.transpose());
  _H.leftCols(mo) = _C.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) + t1);

  // here I'm misusing the transformation function since basically everything it does is:
  // _Bia->col(Q) = jab(Q) * sgl - sgl * jij(Q), which just happens to be part of Bia
  Eigen::Map<Eigen::MatrixXd> sgl(_singles.data(), _nv, _no);
  this->performTransformation(*_Bia, _singles, false, 1, false);

  (*_Bij) = (*_Jij);
  (*_Bia) += (*_Jia);
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia(_Jia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> bia(_Bia->col(Q).data(), _nv, _no);
    Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
    bij.noalias() += jia.transpose() * sgl;
    bia.noalias() -= sgl * jia.transpose() * sgl;
  }
} /* this->transformIntegrals() */

template<>
void CC2Controller<UNRESTRICTED>::transformIntegrals() {
  unsigned iSpin = 0;
  for_spin(_P, _H, _C, _nv, _no) {
    unsigned mo = _no_spin + _nv_spin;
    Eigen::MatrixXd t1 = Eigen::MatrixXd::Zero(mo, mo);
    t1.bottomLeftCorner(_nv_spin, _no_spin) =
        Eigen::Map<Eigen::MatrixXd>(_singles.data() + iSpin * _alpha, _nv_spin, _no_spin);
    _P_spin = _C_spin;
    _H_spin = _C_spin;
    _P_spin.leftCols(mo) = _C_spin.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) - t1.transpose());
    _H_spin.leftCols(mo) = _C_spin.leftCols(mo) * Eigen::MatrixXd(Eigen::MatrixXd::Identity(mo, mo) + t1);
    ++iSpin;
  };

  // here I'm misusing the transformation function since basically everything it does is:
  // _Bia->col(Q) = jab(Q) * sgl - sgl * jij(Q), which just happens to be part of Bia
  Eigen::Map<Eigen::MatrixXd> sgl_a(_singles.data(), _nv_a, _no_a);
  Eigen::Map<Eigen::MatrixXd> sgl_b(_singles.data() + _alpha, _nv_b, _no_b);
  this->performTransformation(_Bia->alpha, _singles.head(_alpha), false, 1, false);
  this->performTransformation(_Bia->beta, _singles.tail(_beta), false, -1, false);
  _Bij->alpha = _Jij->alpha;
  _Bia->alpha += _Jia->alpha;
  _Bij->beta = _Jij->beta;
  _Bia->beta += _Jia->beta;
  for (size_t Q = 0; Q < _nx; ++Q) {
    Eigen::Map<Eigen::MatrixXd> jia_a(_Jia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> bia_a(_Bia->alpha.col(Q).data(), _nv_a, _no_a);
    Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
    bij_a.noalias() += jia_a.transpose() * sgl_a;
    bia_a.noalias() -= sgl_a * jia_a.transpose() * sgl_a;

    Eigen::Map<Eigen::MatrixXd> jia_b(_Jia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> bia_b(_Bia->beta.col(Q).data(), _nv_b, _no_b);
    Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
    bij_b.noalias() += jia_b.transpose() * sgl_b;
    bia_b.noalias() -= sgl_b * jia_b.transpose() * sgl_b;
  }
} /* this->transformIntegrals() */

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity