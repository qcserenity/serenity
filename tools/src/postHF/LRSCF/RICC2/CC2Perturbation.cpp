/**
 * @file CC2Perturbation.cpp
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
std::vector<Eigen::MatrixXd> CC2Controller<RESTRICTED>::getPerturbation(std::vector<double> frequencies,
                                                                        Eigen::Ref<Eigen::MatrixXd> dipoleIntegrals) {
  Timings::takeTime("CC2 -            Perturbation");

  unsigned mo = _no + _nv;
  unsigned nFreq = frequencies.size();
  unsigned nCart = dipoleIntegrals.cols();
  printf("  Computing %3i perturbation vectors.\n\n", nFreq * nCart);
  printf("  root  time(min)\n");
  printf(" ------------------");

  std::vector<Eigen::MatrixXd> rightHandSides(nFreq, Eigen::MatrixXd(_nv * _no, nCart));

  std::vector<SPMatrix<RESTRICTED>> dip_oo(nCart);
  std::vector<SPMatrix<RESTRICTED>> dip_vo(nCart);
  std::vector<SPMatrix<RESTRICTED>> dip_ov(nCart);
  std::vector<SPMatrix<RESTRICTED>> dip_vv(nCart);

  for (unsigned iCart = 0; iCart < nCart; ++iCart) {
    Eigen::Map<Eigen::MatrixXd> MO(dipoleIntegrals.col(iCart).data(), mo, mo);
    dip_oo[iCart] = MO.topLeftCorner(_no, _no);
    dip_vo[iCart] = MO.bottomLeftCorner(_nv, _no);
    dip_ov[iCart] = MO.topRightCorner(_no, _nv).transpose();
    dip_vv[iCart] = MO.bottomRightCorner(_nv, _nv);
  }

  // First part.
  for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
    for (unsigned iCart = 0; iCart < nCart; ++iCart) {
      rightHandSides[iFreq].col(iCart) = Eigen::Map<Eigen::VectorXd>(dip_vo[iCart].data(), dip_vo[iCart].size());
    }
  }

  // Second part.
  for (size_t i = 0; i < _no; ++i) {
    for (size_t j = i; j < _no; ++j) {
      auto tij = this->getAmplitudes(i, j);
      for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
        for (unsigned iCart = 0; iCart < nCart; ++iCart) {
          rightHandSides[iFreq].col(iCart).segment(i * _nv, _nv) += tij * dip_ov[iCart].col(j);
          if (i != j) {
            rightHandSides[iFreq].col(iCart).segment(j * _nv, _nv) += tij.transpose() * dip_ov[iCart].col(i);
          }
        }
      }
    }
  }

  // Third part.
  for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
    double frequency = frequencies[iFreq];
    for (unsigned iCart = 0; iCart < nCart; ++iCart) {
      auto itStart = std::chrono::steady_clock::now();
      _Zia.setZero();
      Eigen::Ref<Eigen::VectorXd> rightHandSide = rightHandSides[iFreq].col(iCart);

      for (size_t i = 0; i < _no; ++i) {
        for (size_t j = i; j < _no; ++j) {
          auto pij = this->getXiDoubles(i, j, frequency, dip_vv[iCart]);
          _Zia.middleRows(i * _nv, _nv).noalias() += pij * _Jia->middleRows(j * _nv, _nv);
          rightHandSide.segment(i * _nv, _nv).noalias() += pij * _Fia.segment(j * _nv, _nv);
          if (i != j) {
            _Zia.middleRows(j * _nv, _nv).noalias() += pij.transpose() * _Jia->middleRows(i * _nv, _nv);
            rightHandSide.segment(j * _nv, _nv).noalias() += pij.transpose() * _Fia.segment(i * _nv, _nv);
          }
        }
      }

      this->reorderTensor(_Fia, _nv, _no);
      this->reorderTensor(rightHandSide, _nv, _no);
      this->reorderTensor(*_Jia, _nv, _no);
      this->reorderTensor(_Zia, _nv, _no);
      if (_Bia != _Jia) {
        this->reorderTensor(*_Bia, _nv, _no);
      }

      for (size_t a = 0; a < _nv; ++a) {
        for (size_t b = a; b < _nv; ++b) {
          auto pab = this->getXiDoublesV(a, b, frequency, dip_oo[iCart]);
          _Zia.middleRows(a * _no, _no).noalias() += pab * _Jia->middleRows(b * _no, _no);
          rightHandSide.middleRows(a * _no, _no).noalias() += pab * _Fia.middleRows(b * _no, _no);
          if (a != b) {
            _Zia.middleRows(b * _no, _no).noalias() += pab.transpose() * _Jia->middleRows(a * _no, _no);
            rightHandSide.middleRows(b * _no, _no).noalias() += pab.transpose() * _Fia.middleRows(a * _no, _no);
          }
        }
      }

      this->reorderTensor(_Fia, _no, _nv);
      this->reorderTensor(rightHandSide, _no, _nv);
      this->reorderTensor(*_Jia, _no, _nv);
      this->reorderTensor(_Zia, _no, _nv);
      if (_Bia != _Jia) {
        this->reorderTensor(*_Bia, _no, _nv);
      }

      Eigen::Map<Eigen::MatrixXd> rhs(rightHandSide.data(), _nv, _no);
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> bij(_Bij->col(Q).data(), _no, _no);
        Eigen::Map<Eigen::MatrixXd> zia(_Zia.col(Q).data(), _nv, _no);
        rhs.noalias() -= zia * bij;
      }

      Eigen::VectorXd nullvec = Eigen::VectorXd::Zero(_nv * _no);
      rightHandSide += this->getJ2GContribution(_Zia.data(), nullptr, nullvec, false);

      auto itEnd = std::chrono::steady_clock::now();
      double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();
      auto str = ((iFreq * nCart + iCart) % 5 == 0) ? "\n %5i %9.3f" : " %5i %9.3f";
      printf(str, iFreq * nCart + iCart + 1, duration / 60.0);
    }

    // Account for the negative sign.
    rightHandSides[iFreq] *= -1.0;
  }
  printf("\n\n");

  Timings::timeTaken("CC2 -            Perturbation");
  return rightHandSides;
}

template<>
std::vector<Eigen::MatrixXd> CC2Controller<UNRESTRICTED>::getPerturbation(std::vector<double> frequencies,
                                                                          Eigen::Ref<Eigen::MatrixXd> dipoleIntegrals) {
  Timings::takeTime("CC2 -            Perturbation");
  size_t mo_a = _no_a + _nv_a;
  size_t mo_b = _no_b + _nv_b;
  unsigned nFreq = frequencies.size();
  unsigned nCart = dipoleIntegrals.cols();
  printf("  Computing %3i perturbation vectors:\n\n", nFreq * nCart);
  printf("  root  time(min)\n");
  printf(" ------------------");

  std::vector<Eigen::MatrixXd> rightHandSides(nFreq, Eigen::MatrixXd(_nDim, nCart));

  std::vector<SPMatrix<UNRESTRICTED>> dip_oo(nCart);
  std::vector<SPMatrix<UNRESTRICTED>> dip_vo(nCart);
  std::vector<SPMatrix<UNRESTRICTED>> dip_ov(nCart);
  std::vector<SPMatrix<UNRESTRICTED>> dip_vv(nCart);

  for (unsigned iCart = 0; iCart < nCart; ++iCart) {
    Eigen::Map<Eigen::MatrixXd> MO_a(dipoleIntegrals.col(iCart).data(), mo_a, mo_a);
    dip_oo[iCart].alpha = MO_a.topLeftCorner(_no_a, _no_a);
    dip_vo[iCart].alpha = MO_a.bottomLeftCorner(_nv_a, _no_a);
    dip_ov[iCart].alpha = MO_a.topRightCorner(_no_a, _nv_a).transpose();
    dip_vv[iCart].alpha = MO_a.bottomRightCorner(_nv_a, _nv_a);

    Eigen::Map<Eigen::MatrixXd> MO_b(dipoleIntegrals.col(iCart).data() + mo_a * mo_a, mo_b, mo_b);
    dip_oo[iCart].beta = MO_b.topLeftCorner(_no_b, _no_b);
    dip_vo[iCart].beta = MO_b.bottomLeftCorner(_nv_b, _no_b);
    dip_ov[iCart].beta = MO_b.topRightCorner(_no_b, _nv_b).transpose();
    dip_vv[iCart].beta = MO_b.bottomRightCorner(_nv_b, _nv_b);
  }

  // First part.
  for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
    for (unsigned iCart = 0; iCart < nCart; ++iCart) {
      rightHandSides[iFreq].col(iCart).head(_alpha) =
          Eigen::Map<Eigen::VectorXd>(dip_vo[iCart].alpha.data(), dip_vo[iCart].alpha.size());
      rightHandSides[iFreq].col(iCart).tail(_beta) =
          Eigen::Map<Eigen::VectorXd>(dip_vo[iCart].beta.data(), dip_vo[iCart].beta.size());
    }
  }

  // Second part.
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = i; j < _no_a; ++j) {
      auto tij = this->getAmplitudes(i, j, 1);
      for (unsigned iCart = 0; iCart < nCart; ++iCart) {
        for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
          rightHandSides[iFreq].block(i * _nv_a, iCart, _nv_a, 1).noalias() += tij * dip_ov[iCart].alpha.col(j);
          if (i != j) {
            rightHandSides[iFreq].block(j * _nv_a, iCart, _nv_a, 1).noalias() +=
                tij.transpose() * dip_ov[iCart].alpha.col(i);
          }
        }
      }
    }
  }
  for (size_t i = 0; i < _no_b; ++i) {
    for (size_t j = i; j < _no_b; ++j) {
      auto tij = this->getAmplitudes(i, j, -1);
      for (unsigned iCart = 0; iCart < nCart; ++iCart) {
        for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
          rightHandSides[iFreq].block(_alpha + i * _nv_b, iCart, _nv_b, 1).noalias() += tij * dip_ov[iCart].beta.col(j);
          if (i != j) {
            rightHandSides[iFreq].block(_alpha + j * _nv_b, iCart, _nv_b, 1).noalias() +=
                tij.transpose() * dip_ov[iCart].beta.col(i);
          }
        }
      }
    }
  }
  for (size_t i = 0; i < _no_a; ++i) {
    for (size_t j = 0; j < _no_b; ++j) {
      auto tij = this->getAmplitudes(i, j, 0);
      for (unsigned iCart = 0; iCart < nCart; ++iCart) {
        for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
          rightHandSides[iFreq].block(i * _nv_a, iCart, _nv_a, 1).noalias() += tij * dip_ov[iCart].beta.col(j);
          rightHandSides[iFreq].block(_alpha + j * _nv_b, iCart, _nv_b, 1).noalias() +=
              tij.transpose() * dip_ov[iCart].alpha.col(i);
        }
      }
    }
  }

  // Third part.
  for (unsigned iFreq = 0; iFreq < nFreq; ++iFreq) {
    double frequency = frequencies[iFreq];
    for (unsigned iCart = 0; iCart < nCart; ++iCart) {
      auto itStart = std::chrono::steady_clock::now();
      _Zia.alpha.setZero();
      _Zia.beta.setZero();
      Eigen::Ref<Eigen::VectorXd> rightHandSide_a = rightHandSides[iFreq].col(iCart).head(_alpha);
      Eigen::Ref<Eigen::VectorXd> rightHandSide_b = rightHandSides[iFreq].col(iCart).tail(_beta);

      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = i; j < _no_a; ++j) {
          Eigen::MatrixXd pij = this->getXiDoubles(i, j, frequency, dip_vv[iCart], 1);
          _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += pij * _Jia->alpha.middleRows(j * _nv_a, _nv_a);
          rightHandSide_a.middleRows(i * _nv_a, _nv_a).noalias() += pij * _Fia.middleRows(j * _nv_a, _nv_a);
          if (i != j) {
            _Zia.alpha.middleRows(j * _nv_a, _nv_a).noalias() += pij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
            rightHandSide_a.middleRows(j * _nv_a, _nv_a).noalias() += pij.transpose() * _Fia.middleRows(i * _nv_a, _nv_a);
          }
        }
      }
      for (size_t i = 0; i < _no_b; ++i) {
        for (size_t j = i; j < _no_b; ++j) {
          Eigen::MatrixXd pij = this->getXiDoubles(i, j, frequency, dip_vv[iCart], -1);
          _Zia.beta.middleRows(i * _nv_b, _nv_b).noalias() += pij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
          rightHandSide_b.middleRows(i * _nv_b, _nv_b).noalias() += pij * _Fia.middleRows(_alpha + j * _nv_b, _nv_b);
          if (i != j) {
            _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += pij.transpose() * _Jia->beta.middleRows(i * _nv_b, _nv_b);
            rightHandSide_b.middleRows(j * _nv_b, _nv_b).noalias() +=
                pij.transpose() * _Fia.middleRows(_alpha + i * _nv_b, _nv_b);
          }
        }
      }
      for (size_t i = 0; i < _no_a; ++i) {
        for (size_t j = 0; j < _no_b; ++j) {
          auto pij = this->getXiDoubles(i, j, frequency, dip_vv[iCart], 0);
          _Zia.alpha.middleRows(i * _nv_a, _nv_a).noalias() += pij * _Jia->beta.middleRows(j * _nv_b, _nv_b);
          rightHandSide_a.middleRows(i * _nv_a, _nv_a).noalias() += pij * _Fia.middleRows(_alpha + j * _nv_b, _nv_b);

          _Zia.beta.middleRows(j * _nv_b, _nv_b).noalias() += pij.transpose() * _Jia->alpha.middleRows(i * _nv_a, _nv_a);
          rightHandSide_b.middleRows(j * _nv_b, _nv_b).noalias() += pij.transpose() * _Fia.middleRows(i * _nv_a, _nv_a);
        }
      }

      this->reorderTensor(_Fia.head(_alpha), _nv_a, _no_a);
      this->reorderTensor(_Fia.tail(_beta), _nv_b, _no_b);
      this->reorderTensor(rightHandSide_a, _nv_a, _no_a);
      this->reorderTensor(rightHandSide_b, _nv_b, _no_b);
      this->reorderTensor(_Jia->alpha, _nv_a, _no_a);
      this->reorderTensor(_Jia->beta, _nv_b, _no_b);
      this->reorderTensor(_Zia.alpha, _nv_a, _no_a);
      this->reorderTensor(_Zia.beta, _nv_b, _no_b);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->alpha, _nv_a, _no_a);
        this->reorderTensor(_Bia->beta, _nv_b, _no_b);
      }

      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = a; b < _nv_a; ++b) {
          Eigen::MatrixXd pab = this->getXiDoublesV(a, b, frequency, dip_oo[iCart], 1);
          _Zia.alpha.middleRows(a * _no_a, _no_a).noalias() += pab * _Jia->alpha.middleRows(b * _no_a, _no_a);
          rightHandSide_a.middleRows(a * _no_a, _no_a).noalias() += pab * _Fia.middleRows(b * _no_a, _no_a);
          if (a != b) {
            _Zia.alpha.middleRows(b * _no_a, _no_a).noalias() += pab.transpose() * _Jia->alpha.middleRows(a * _no_a, _no_a);
            rightHandSide_a.middleRows(b * _no_a, _no_a).noalias() += pab.transpose() * _Fia.middleRows(a * _no_a, _no_a);
          }
        }
      }
      for (size_t a = 0; a < _nv_b; ++a) {
        for (size_t b = a; b < _nv_b; ++b) {
          Eigen::MatrixXd pab = this->getXiDoublesV(a, b, frequency, dip_oo[iCart], -1);
          _Zia.beta.middleRows(a * _no_b, _no_b).noalias() += pab * _Jia->beta.middleRows(b * _no_b, _no_b);
          rightHandSide_b.middleRows(a * _no_b, _no_b).noalias() += pab * _Fia.middleRows(_alpha + b * _no_b, _no_b);
          if (a != b) {
            _Zia.beta.middleRows(b * _no_b, _no_b).noalias() += pab.transpose() * _Jia->beta.middleRows(a * _no_b, _no_b);
            rightHandSide_b.middleRows(b * _no_b, _no_b).noalias() +=
                pab.transpose() * _Fia.middleRows(_alpha + a * _no_b, _no_b);
          }
        }
      }
      for (size_t a = 0; a < _nv_a; ++a) {
        for (size_t b = 0; b < _nv_b; ++b) {
          auto pab = this->getXiDoublesV(a, b, frequency, dip_oo[iCart], 0);
          _Zia.alpha.middleRows(a * _no_a, _no_a).noalias() += pab * _Jia->beta.middleRows(b * _no_b, _no_b);
          rightHandSide_a.middleRows(a * _no_a, _no_a).noalias() += pab * _Fia.middleRows(_alpha + b * _no_b, _no_b);

          _Zia.beta.middleRows(b * _no_b, _no_b).noalias() += pab.transpose() * _Jia->alpha.middleRows(a * _no_a, _no_a);
          rightHandSide_b.middleRows(b * _no_b, _no_b).noalias() += pab.transpose() * _Fia.middleRows(a * _no_a, _no_a);
        }
      }

      this->reorderTensor(_Fia.head(_alpha), _no_a, _nv_a);
      this->reorderTensor(_Fia.tail(_beta), _no_b, _nv_b);
      this->reorderTensor(rightHandSide_a, _no_a, _nv_a);
      this->reorderTensor(rightHandSide_b, _no_b, _nv_b);
      this->reorderTensor(_Jia->alpha, _no_a, _nv_a);
      this->reorderTensor(_Jia->beta, _no_b, _nv_b);
      this->reorderTensor(_Zia.alpha, _no_a, _nv_a);
      this->reorderTensor(_Zia.beta, _no_b, _nv_b);
      if (_Bia != _Jia) {
        this->reorderTensor(_Bia->alpha, _no_a, _nv_a);
        this->reorderTensor(_Bia->beta, _no_b, _nv_b);
      }

      Eigen::Map<Eigen::MatrixXd> rhs_a(rightHandSide_a.data(), _nv_a, _no_a);
      Eigen::Map<Eigen::MatrixXd> rhs_b(rightHandSide_b.data(), _nv_b, _no_b);
      for (size_t Q = 0; Q < _nx; ++Q) {
        Eigen::Map<Eigen::MatrixXd> bij_a(_Bij->alpha.col(Q).data(), _no_a, _no_a);
        Eigen::Map<Eigen::MatrixXd> zia_a(_Zia.alpha.col(Q).data(), _nv_a, _no_a);
        rhs_a.noalias() -= zia_a * bij_a;

        Eigen::Map<Eigen::MatrixXd> bij_b(_Bij->beta.col(Q).data(), _no_b, _no_b);
        Eigen::Map<Eigen::MatrixXd> zia_b(_Zia.beta.col(Q).data(), _nv_b, _no_b);
        rhs_b.noalias() -= zia_b * bij_b;
      }

      Eigen::VectorXd nullvec = Eigen::VectorXd::Zero(_alpha + _beta);
      rightHandSide_a += this->getJ2GContribution(_Zia.alpha.data(), nullptr, nullvec.head(_alpha), false, 1);
      rightHandSide_b += this->getJ2GContribution(_Zia.beta.data(), nullptr, nullvec.tail(_beta), false, -1);

      auto itEnd = std::chrono::steady_clock::now();
      double duration = std::chrono::duration_cast<std::chrono::duration<double>>(itEnd - itStart).count();
      auto str = ((iFreq * nCart + iCart) % 5 == 0) ? "\n %5i %9.3f" : " %5i %9.3f";
      printf(str, iFreq * nCart + iCart + 1, duration / 60.0);
    }

    // Account for the negative sign.
    rightHandSides[iFreq] *= -1.0;
  }
  printf("\n");

  Timings::timeTaken("CC2 -            Perturbation");
  return rightHandSides;
}

template class CC2Controller<Options::SCF_MODES::RESTRICTED>;
template class CC2Controller<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity