/**
 * @file ArithmeticSeriesDamping.h
 *
 * @date Oct 28, 2016
 * @author M. Boeckers
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

#ifndef ELECTRONICSTRUCTURECALCULATIONS_SCF_DAMPER_ARITHMETICSERIESDAMPING_H_
#define ELECTRONICSTRUCTURECALCULATIONS_SCF_DAMPER_ARITHMETICSERIESDAMPING_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "scf/damper/Damper.h"

namespace Serenity {
template<Options::SCF_MODES T>
class ArithmeticSeriesDamping : public Damper<T> {
 public:
  /**
   * @brief           Damping with varying damping factor in a given range.
   *
   *                  Damping with varying damping factor in range dStart to dEnd. Note that
   *                  the end value must not be reached. The damping will stop when the damping
   *                  factor falls below dEnd. Negative damping factors are not possible.
   *                  If iStartUp > 0, it will perform iStartUp static damping steps at the
   *                  beginning of the SCF calculation.  It can make sense to use a higher
   *                  damping factor multiple times at the beginning because the Fock matrix
   *                  varies strongly in the first few iterations.
   *                  See Damper.h for more information.
   *
   * @param dStart    Start damping factor (dStart > 0.0)
   * @param dStep     Step width (dStep > 0.0)
   * @param dEnd      End damping factor (dEnd > 0.0)
   * @param iStartUp  Number of static damping steps with dStart at the beginning of the SCF
   *
   */
  ArithmeticSeriesDamping(const double dStart, const double dStep, const double dEnd, int iStartUp);

  virtual ~ArithmeticSeriesDamping() = default;

  /**
   *
   * @param newMatrix  The matrix to be damped. a part of the 'old matrix' is mixed into this.
   *                   The result will also become the new 'old matrix'.
   * @return           Returns true if damping is finished
   */
  virtual void damp(SpinPolarizedData<T, Eigen::MatrixXd>& newMatrix) override final;
  /**
   *
   * @param newMatrix  The matrix to be damped. a part of the 'old matrix' is mixed into this.
   *                   The result will also become the new 'old matrix'.
   * @return           Returns true if damping is finished
   */
  virtual void damp(FockMatrix<T>& newMatrix) override final;

  ///@brief Just override.
  virtual void dynamicDamp(FockMatrix<T>&, DensityMatrix<T>) override final{};

 private:
  const double _dStart;
  const double _dStep;
  const double _dEnd;
  int _iStartUp;
  double _dampingFactor;
  bool _initialized;
  SpinPolarizedData<T, Eigen::MatrixXd> _oldMatrix;
};
} /* namespace Serenity */
#endif /* ELECTRONICSTRUCTURECALCULATIONS_SCF_DAMPER_ARITHMETICSERIESDAMPING_H_ */
