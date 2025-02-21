/**
 * @file DeltaSpinSquared.cpp
 *
 * @date Jul. 7, 2021
 * @author Johannes Toelle
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
#include "postHF/LRSCF/Analysis/DeltaSpinSquared.h"
/* Include Serenity Internal Headers */
#include "integrals/OneElectronIntegralController.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
DeltaSpinSquared<SCFMode>::DeltaSpinSquared(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                            unsigned int nEigen, Options::LR_METHOD method, Options::LRSCF_TYPE type)
  : _lrscf(lrscf), _nEigen(nEigen), _method(method), _type(type) {
}

template<>
double DeltaSpinSquared<Options::SCF_MODES::UNRESTRICTED>::calculateSpinSquared(unsigned int iState) {
  double deltaSsquared = 0.0;
  for (auto lrscf : _lrscf) {
    auto excVec = *(lrscf->getExcitationVectors(_type));
    auto nOcc = lrscf->getNOccupied();
    auto nVirt = lrscf->getNVirtual();
    auto coef = lrscf->getCoefficients();
    auto sys = lrscf->getSys();
    auto overlap = sys->getOneElectronIntegralController()->getOverlapIntegrals();
    unsigned ia_a = nOcc.alpha * nVirt.alpha;

    // (X+Y)_alpha
    Eigen::MatrixXd xpy_a = (Eigen::Map<Eigen::MatrixXd>(excVec[0].col(iState).data(), nVirt.alpha, nOcc.alpha)).transpose();
    // (X+Y)_beta
    Eigen::MatrixXd xpy_b =
        (Eigen::Map<Eigen::MatrixXd>(excVec[0].col(iState).data() + ia_a, nVirt.beta, nOcc.beta)).transpose();

    // (X-Y)_alpha
    Eigen::MatrixXd xmy_a = (Eigen::Map<Eigen::MatrixXd>(excVec[1].col(iState).data(), nVirt.alpha, nOcc.alpha)).transpose();
    // (X-Y)_beta
    Eigen::MatrixXd xmy_b =
        (Eigen::Map<Eigen::MatrixXd>(excVec[1].col(iState).data() + ia_a, nVirt.beta, nOcc.beta)).transpose();

    Eigen::MatrixXd x_a = 0.5 * (xpy_a + xmy_a);
    Eigen::MatrixXd x_b = 0.5 * (xpy_b + xmy_b);

    Eigen::MatrixXd y_a = 0.5 * (xpy_a - xmy_a);
    Eigen::MatrixXd y_b = 0.5 * (xpy_b - xmy_b);

    Eigen::MatrixXd hole_aa = x_a * x_a.transpose();
    Eigen::MatrixXd hole_bb = x_b * x_b.transpose();
    Eigen::MatrixXd particle_aa = x_a.transpose() * x_a;
    Eigen::MatrixXd particle_bb = x_b.transpose() * x_b;
    Eigen::MatrixXd overlapMOab = coef.alpha.transpose() * overlap * coef.beta;

    Eigen::MatrixXd overlap_ba_virt = overlapMOab.transpose().block(0, nOcc.alpha, nOcc.beta, nVirt.alpha);
    Eigen::MatrixXd overlap_ba_occ = overlapMOab.transpose().block(0, 0, nOcc.beta, nOcc.alpha);
    Eigen::MatrixXd overlap_ab_virt = overlapMOab.block(0, nOcc.beta, nOcc.alpha, nVirt.beta);
    Eigen::MatrixXd overlap_ab_occ = overlapMOab.block(0, 0, nOcc.alpha, nOcc.beta);

    double contr1 = (overlap_ba_occ * hole_aa * overlap_ba_occ.transpose()).trace();
    double contr2 = (overlap_ba_virt * particle_aa * overlap_ba_virt.transpose()).trace();
    double contr3 = (overlap_ab_occ * hole_bb * overlap_ab_occ.transpose()).trace();
    double contr4 = (overlap_ab_virt * particle_bb * overlap_ab_virt.transpose()).trace();
    double contr5 = 2.0 * ((overlapMOab.transpose().leftCols(nOcc.alpha) * x_a *
                            overlapMOab.block(nOcc.alpha, 0, nVirt.alpha, overlapMOab.cols()))
                               .block(0, nOcc.beta, nOcc.beta, nVirt.beta) *
                           x_b.transpose())
                              .trace();
    deltaSsquared = contr1 - contr2 + contr3 - contr4 - contr5;
    if (_method == Options::LR_METHOD::TDDFT) {
      hole_aa = y_a * y_a.transpose();
      hole_bb = y_b * y_b.transpose();
      particle_aa = y_a.transpose() * y_a;
      particle_bb = y_b.transpose() * y_b;
      contr1 = (overlap_ba_occ * hole_aa * overlap_ba_occ.transpose()).trace();
      contr2 = (overlap_ba_virt * particle_aa * overlap_ba_virt.transpose()).trace();
      contr3 = (overlap_ab_occ * hole_bb * overlap_ab_occ.transpose()).trace();
      contr4 = (overlap_ab_virt * particle_bb * overlap_ab_virt.transpose()).trace();
      contr5 = 2.0 * ((overlapMOab.transpose().leftCols(nOcc.alpha) * y_a *
                       overlapMOab.block(nOcc.alpha, 0, nVirt.alpha, overlapMOab.cols()))
                          .block(0, nOcc.beta, nOcc.beta, nVirt.beta) *
                      y_b.transpose())
                         .trace();
      deltaSsquared = deltaSsquared - contr1 + contr2 - contr3 + contr4 + contr5;
      // mixed terms
      double mixed1 = 2.0 * ((overlapMOab.leftCols(nOcc.beta) * y_b *
                              overlapMOab.transpose().block(nOcc.beta, 0, nVirt.beta, overlapMOab.cols()))
                                 .block(0, nOcc.alpha, nOcc.alpha, nVirt.alpha) *
                             x_a.transpose())
                                .trace();
      double mixed2 = 2.0 * ((overlapMOab.transpose().leftCols(nOcc.alpha) * y_a *
                              overlapMOab.block(nOcc.alpha, 0, nVirt.alpha, overlapMOab.cols()))
                                 .block(0, nOcc.beta, nOcc.beta, nVirt.beta) *
                             x_b.transpose())
                                .trace();
      deltaSsquared += mixed1 - mixed2;
    }
  }
  return deltaSsquared;
}

template<>
double DeltaSpinSquared<Options::SCF_MODES::RESTRICTED>::calculateSpinSquared(unsigned int iState) {
  (void)iState;
  std::cout << " --- Delta S^2 not supported for Restricted! --- " << std::endl;
  return 0.0;
}

template<Options::SCF_MODES SCFMode>
void DeltaSpinSquared<SCFMode>::print() {
  _spinSquared = Eigen::VectorXd::Zero(_nEigen);
  printf("                                    Delta <S*S>                                        \n");
  printf(" %5s %12s \n", "state", "Delta <S*S>");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned int i = 0; i < _nEigen; i++) {
    _spinSquared(i) = this->calculateSpinSquared(i);
    printf(" %4i %11.6f ", i + 1, _spinSquared(i));
    printf("\n");
  }
  printf("---------------------------------------------------------------------------------------\n");
}

template class DeltaSpinSquared<Options::SCF_MODES::RESTRICTED>;
template class DeltaSpinSquared<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
