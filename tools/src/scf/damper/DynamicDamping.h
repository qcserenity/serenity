/**
 * @file   DynamicDamping.h
 *
 * @date   13 August 2020
 * @author Niklas Niemeyer
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
#ifndef DYNAMICDAMPING_H
#define DYNAMICDAMPING_H
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "scf/damper/Damper.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace Serenity {
/**
 * @class DynamicDamping DynamicDamping.h
 * @brief Employs the 'Optimal Damping Algorithm' [Cancés et. al, Int. J. Quantum Chem. 79, 82-90 (2000)].
 */
template<Options::SCF_MODES SCFMode>
class DynamicDamping : public Damper<SCFMode> {
 public:
  /**
   * @brief Dynamic damping based on an energy-like criterion.
   */
  explicit DynamicDamping();

  /**
   * @brief Default destructor.
   */
  virtual ~DynamicDamping() = default;

  /**
   * @brief Calculates the damping factor and updates the Fock matrix accordingly.
   * @param newFock     A reference of the Fock matrix to be damped. Is updated and stored as _oldFock.
   * @param newDensity  A copy of the density matrix to be damped. Is also updated, so it can be stored as _oldDensity.
   */
  void dynamicDamp(FockMatrix<SCFMode>& newFock, DensityMatrix<SCFMode> newDensity) override;

  ///@brief Not used here, just override.
  virtual void damp(FockMatrix<SCFMode>&) override final{};

  ///@brief Not used here, just override.
  virtual void damp(SpinPolarizedData<SCFMode, Eigen::MatrixXd>&) override final{};

 private:
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _oldFock;
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _oldDensity;

  bool _initialized;
};

} /* namespace Serenity */
#endif /* DYNAMICDAMPING_H */
