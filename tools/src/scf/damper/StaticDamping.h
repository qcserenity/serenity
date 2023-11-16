/**
 * @file   StaticDamping.h
 * @date   29. Dezember 2013, 14:13
 * @author Thomas Dresselhaus, Jan Unsleber
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
#ifndef STATICDAMPING_H
#define STATICDAMPING_H
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "scf/damper/Damper.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace Serenity {
/* Forward declarations */
/**
 * @class StaticDamping StaticDamping.h
 * @brief This class uses the simplest way of improving convergence of an SCF.
 */
template<Options::SCF_MODES T>
class StaticDamping : public Damper<T> {
 public:
  /**
   * @param dampingFactor The amount of 'old matrix' that will be mixed in. 0.0 would
   *                      have no effect on the newMatrix, 1.0 would completely overwrite
   *                      the newMatrix with the old one. A negative dampingFactor is also
   *                      possible but risky. It would mean that the changes (difference
   *                      between old and new matrix) are enhanced instead of being reduced.
   *
   * @brief               Perform static damping with dampingFactor for maxCycles.
   *                      See Damper.h for more information.
   */
  explicit StaticDamping(const double dampingFactor);
  virtual ~StaticDamping() = default;

  /**
   * @brief see above. Only with Eigen-matrix.
   *
   * @param newMatrix     the matrix to be damped, a part of the 'old matrix' is mixed
   *                      into this. The result will also become the new 'old matrix'.
   */
  virtual void damp(FockMatrix<T>& newMatrix) override final;
  virtual void damp(SpinPolarizedData<T, Eigen::MatrixXd>& newMatrix) override final;

  ///@brief Just override.
  virtual void dynamicDamp(FockMatrix<T>&, DensityMatrix<T>) override final{};

 private:
  SpinPolarizedData<T, Eigen::MatrixXd> _oldMatrix;
  const double _dampingFactor;

  bool _initialized;
};

} /* namespace Serenity */
#endif /* STATIC_DAMPING_H */
