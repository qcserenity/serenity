/**
 * @file NROCalculator.h
 *
 * @date Apr 23, 2021
 * @author Anton Rikus
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
#ifndef NROCALCULATOR_H_
#define NROCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class  NROCalculator NROCalculator.h
 * @brief  Class to calculate natural response orbitals as pairs of particle and hole orbitals.
 *
 * The theory can be found in: J. Chem. Theory Comput, 16, 7709-7720 (2020) "A Unified Strategy for the Chemically
 * Intuitive Interpretation of Molecular Optical Response Properties"
 */
template<Options::SCF_MODES SCFMode>
class NROCalculator {
 public:
  /**
   * @brief Constructor.
   * @param solutionvectors        A vector containing X+Y as first and X-Y as second entry, thus representing the
   *                               solution to the response problenm
   * @param lrscfController        An LRSCFController which can be used to gain information on the LRSCFTask
   */
  NROCalculator(std::vector<Eigen::MatrixXd> solutionvectors, std::shared_ptr<LRSCFController<SCFMode>> lrscfController);

  /**
   * @brief Default destructor.
   */
  virtual ~NROCalculator() = default;

  /**
   * @brief A function to calculate the natural response orbitals.
   * @param iFreq The index of the perturbance frequency (same index as in the input-File)
   * @return  A spinpolarized vector with six entries, where the first entry contains the hole-NROs in x-direction, the
   *          second entry describes the particle-NROs in x-direction, the third entry contains the hole-NROs in y and
   *          so on. The NROs are returned in AO basis.
   */
  SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> getNROs(unsigned int iFreq);

  /**
   * @brief Get computed singular values of the NROs calculated by this class.
   * @return Singular values.
   */
  SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> getSingularValues();

 private:
  /// @brief An LRSCFController which holds information on the LRSCFTask
  std::shared_ptr<LRSCFController<SCFMode>> _lrscfcontroller;
  /// @brief A vector of matrices that contains the matrix X+Y as the first entry and X-Y as the second entry.
  std::vector<Eigen::MatrixXd> _XY;
  /// @brief The singular values of the response density matrix.
  SpinPolarizedData<SCFMode, std::vector<Eigen::MatrixXd>> _singularValues;
};

} /* namespace Serenity */

#endif /* NROCALCULATOR_H_ */