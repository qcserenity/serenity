/**
 * @file OrbitalPairDIISWrapper.cpp
 *
 * @date Jul 30, 2019
 * @author Moritz Bensberg
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
#include "postHF/LocalCorrelation/OrbitalPairDIISWrapper.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalPair.h"        //Orbital pair definition.
#include "data/SingleSubstitution.h" //Singles definition.

namespace Serenity {

OrbitalPairDIISWrapper::OrbitalPairDIISWrapper(unsigned int maxStore) : _diis(maxStore) {
}

void OrbitalPairDIISWrapper::optimize(std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs,
                                      std::vector<std::shared_ptr<SingleSubstitution>> singles) {
  unsigned int nAmplitudes = 0;
  for (const auto& pair : orbitalPairs)
    nAmplitudes += pair->k_ij.cols() * pair->k_ij.cols();
  for (const auto& single : singles)
    nAmplitudes += single->t_i.rows();
  Eigen::VectorXd amplitudes(nAmplitudes);
  Eigen::VectorXd residuals(nAmplitudes);
  unsigned int row = 0;
  for (const auto& pair : orbitalPairs) {
    unsigned int nAmpl = pair->k_ij.cols() * pair->k_ij.cols();
    amplitudes.segment(row, nAmpl) = Eigen::Map<Eigen::VectorXd>(pair->t_ij.data(), nAmpl);
    Eigen::MatrixXd update = Eigen::MatrixXd(-pair->residual.array() / pair->uncoupledTerm.array());
    residuals.segment(row, nAmpl) = Eigen::Map<Eigen::VectorXd>(update.data(), nAmpl);
    row += nAmpl;
  }
  for (const auto& single : singles) {
    unsigned int nAmpl = single->t_i.rows();
    amplitudes.segment(row, nAmpl) = single->t_i;
    Eigen::VectorXd update = Eigen::MatrixXd(-single->residual.array() / single->epsMinusF.array());
    residuals.segment(row, nAmpl) = update;
    row += nAmpl;
  }
  _diis.optimize(amplitudes, residuals);
  row = 0;
  for (auto& pair : orbitalPairs) {
    unsigned int nPNOs = pair->k_ij.cols();
    unsigned int nAmpl = nPNOs * nPNOs;
    pair->t_ij = Eigen::Map<Eigen::MatrixXd>(amplitudes.segment(row, nAmpl).data(), nPNOs, nPNOs);
    row += nAmpl;
  }
  for (auto& single : singles) {
    unsigned int nAmpl = single->t_i.rows();
    single->t_i = amplitudes.segment(row, nAmpl);
    row += nAmpl;
  }
}

} /* namespace Serenity */
