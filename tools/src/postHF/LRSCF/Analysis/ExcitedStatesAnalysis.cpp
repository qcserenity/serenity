/**
 * @file ExcitedStatesAnalysis.cpp
 *
 * @date Apr 21, 2021
 * @author Anton Rikus
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
#include "postHF/LRSCF/Analysis/ExcitedStatesAnalysis.h"
/* Include Serenity Internal Headers */
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ExcitedStatesAnalysis<SCFMode>::ExcitedStatesAnalysis(std::vector<unsigned int> excitations,
                                                      std::shared_ptr<LRSCFController<SCFMode>> lrscfc)
  : _excitations(excitations), _lrscfcontroller(lrscfc) {
  for (unsigned i = 0; i < _excitations.size(); i++) {
    if (_lrscfcontroller->getLRSCFSettings().nEigen < _excitations[i])
      throw SerenityError("You need to perform an LRSCFTask with neigen greater or equal to the PlotTask's highest "
                          "excitationnumber!");
    _transdensmatrix.push_back(nullptr);
    _holedensmatrix.push_back(nullptr);
    _particledensmatrix.push_back(nullptr);
  }
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<MatrixInBasis<SCFMode>> ExcitedStatesAnalysis<SCFMode>::getTransitionDensityMatrix(unsigned int n) {
  if (!_transdensmatrix[n]) {
    this->calculateTransitionDensityMatrix(n);
  }
  return (_transdensmatrix[n]);
};

template<Options::SCF_MODES SCFMode>
std::shared_ptr<MatrixInBasis<SCFMode>> ExcitedStatesAnalysis<SCFMode>::getHoleDensityMatrix(unsigned int n) {
  if (!_holedensmatrix[n]) {
    this->calculateHoleDensityMatrix(n);
  }
  return (_holedensmatrix[n]);
};

template<Options::SCF_MODES SCFMode>
std::shared_ptr<MatrixInBasis<SCFMode>> ExcitedStatesAnalysis<SCFMode>::getParticleDensityMatrix(unsigned int n) {
  if (!_particledensmatrix[n]) {
    this->calculateParticleDensityMatrix(n);
  }
  return (_particledensmatrix[n]);
};

template<Options::SCF_MODES SCFMode>
void ExcitedStatesAnalysis<SCFMode>::calculateTransitionDensityMatrix(unsigned int n) {
  auto nocc = _lrscfcontroller->getNOccupied();
  auto nvirt = _lrscfcontroller->getNVirtual();
  auto& coeffmat = _lrscfcontroller->getCoefficients();
  const auto& excvectors = _lrscfcontroller->getExcitationVectors(Options::LRSCF_TYPE::ISOLATED);
  unsigned ia_start = 0;
  _transdensmatrix[n] = std::make_shared<MatrixInBasis<SCFMode>>(_lrscfcontroller->getBasisController());
  auto& temp = *_transdensmatrix[n];
  unsigned exc = _excitations[n] - 1;
  Eigen::VectorXd XpY = (*excvectors)[0].col(exc) + (*excvectors)[1].col(exc);
  for_spin(nocc, nvirt, coeffmat, temp) {
    Eigen::MatrixXd xpy = Eigen::Map<Eigen::MatrixXd>(XpY.data() + ia_start, nvirt_spin, nocc_spin);
    temp_spin =
        coeffmat_spin.leftCols(nocc_spin) * xpy.transpose() * coeffmat_spin.middleCols(nocc_spin, nvirt_spin).transpose();
    ia_start += nocc_spin * nvirt_spin;
  };
}

template<Options::SCF_MODES SCFMode>
void ExcitedStatesAnalysis<SCFMode>::calculateHoleDensityMatrix(unsigned int n) {
  auto nocc = _lrscfcontroller->getNOccupied();
  auto nvirt = _lrscfcontroller->getNVirtual();
  const auto& coeffmat = _lrscfcontroller->getCoefficients();
  const auto& excvectors = _lrscfcontroller->getExcitationVectors(Options::LRSCF_TYPE::ISOLATED);
  unsigned ia_start = 0;
  _holedensmatrix[n] = std::make_shared<MatrixInBasis<SCFMode>>(_lrscfcontroller->getBasisController());
  auto& temp = *_holedensmatrix[n];
  unsigned exc = _excitations[n] - 1;
  Eigen::VectorXd XpY = (*excvectors)[0].col(exc) + (*excvectors)[1].col(exc);
  Eigen::VectorXd XmY = (*excvectors)[0].col(exc) - (*excvectors)[1].col(exc);
  for_spin(nocc, nvirt, coeffmat, temp) {
    Eigen::MatrixXd xpy = Eigen::Map<Eigen::MatrixXd>(XpY.data() + ia_start, nvirt_spin, nocc_spin);
    Eigen::MatrixXd xmy = Eigen::Map<Eigen::MatrixXd>(XmY.data() + ia_start, nvirt_spin, nocc_spin);
    temp_spin = -0.5 * coeffmat_spin.leftCols(nocc_spin) * (xpy.transpose() * xmy + xmy.transpose() * xpy) *
                coeffmat_spin.leftCols(nocc_spin).transpose();
    ia_start += nocc_spin * nvirt_spin;
  };
}

template<Options::SCF_MODES SCFMode>
void ExcitedStatesAnalysis<SCFMode>::calculateParticleDensityMatrix(unsigned int n) {
  auto nocc = _lrscfcontroller->getNOccupied();
  auto nvirt = _lrscfcontroller->getNVirtual();
  const auto& coeffmat = _lrscfcontroller->getCoefficients();
  const auto& excvectors = _lrscfcontroller->getExcitationVectors(Options::LRSCF_TYPE::ISOLATED);
  unsigned ia_start = 0;
  _particledensmatrix[n] = std::make_shared<MatrixInBasis<SCFMode>>(_lrscfcontroller->getBasisController());
  auto& temp = *_particledensmatrix[n];
  unsigned exc = _excitations[n] - 1;
  Eigen::VectorXd XpY = (*excvectors)[0].col(exc) + (*excvectors)[1].col(exc);
  Eigen::VectorXd XmY = (*excvectors)[0].col(exc) - (*excvectors)[1].col(exc);
  for_spin(nocc, nvirt, coeffmat, temp) {
    Eigen::MatrixXd xpy = Eigen::Map<Eigen::MatrixXd>(XpY.data() + ia_start, nvirt_spin, nocc_spin);
    Eigen::MatrixXd xmy = Eigen::Map<Eigen::MatrixXd>(XmY.data() + ia_start, nvirt_spin, nocc_spin);
    temp_spin = 0.5 * coeffmat_spin.middleCols(nocc_spin, nvirt_spin) * (xmy * xpy.transpose() + xpy * xmy.transpose()) *
                coeffmat_spin.middleCols(nocc_spin, nvirt_spin).transpose();
    ia_start += nocc_spin * nvirt_spin;
  };
}

template class ExcitedStatesAnalysis<Options::SCF_MODES::RESTRICTED>;
template class ExcitedStatesAnalysis<Options::SCF_MODES::UNRESTRICTED>;

} // namespace Serenity