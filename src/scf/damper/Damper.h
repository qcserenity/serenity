/**
 * @file   Damper.h
 *
 * @date   Oct 10, 2016, reworked on Feb 07, 2024
 * @author Jan Unsleber, rework by Lukas Paetow
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
#ifndef MATHS_DAMPER_DAMPER_H_
#define MATHS_DAMPER_DAMPER_H_
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "settings/Options.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
/**
 * @brief Damping is one of the earliest and most used methods to obtain SCF convergence.
 *        A new matrix is constructed from a weighted average of the matrix of the last iteration and
 *        the new one
 *        \f[
 *          M^{(n+1)\prime} = (1-\lambda) M^{(n+1)} + \lambda M^{(n)}
 *        \f]
 *        where \f$ 0 \leq \lambda < 1 \f$. In the SCF, damping can be applied to the the density matrix,
 *        the Fock matrix or the MO coefficients. If applied to the MO coefficients, one must take care of
 *        the signs (phases) and the mixing of degenerate orbitals.
 *        During the first iterations of the SCF, damping has the effect to reduce the oscillations
 *        in the matrices. In the first iterations, \f$ \lambda \f$ can be constant or determined
 *        dynamically until the convergence falls below some predefined threshold after which
 *        \f$ \lambda \f$ is set to zero.
 *        \n
 *        Note that damping does not preserve the orthonormality of the MOs or the indempotency of the
 *        density matrix. However, this is not a severe problem since the damping is switched off in later
 *        SCF iterations. Furthermore, damping may decrease the speed of convergence, but slow convergence
 *        is better than no convergence if no better technique is available.
 *        \n
 *        Other quantum chemistry programs may use the slightly different damping relation
 *        \f[
 *          M^{(n+1)\prime} = M^{(n+1)} + \lambda M^{(n)}
 *        \f]
 *        with \f$ \lambda \in \mathrm{R} \f$ .
 *        \n
 *        For further information about damping, see actual implementations.
 *        \n
 *        Ref.: See e.g. H. B. Schlegel and J. J. W. McDouall, Do you have SCF Stability and Convergence
 *              Problems? in Computational Advances in Organic Chemistry: Molecular Structure and Reactivity,
 *              167-185, 1991 or some basic text books.
 *
 */
class Damper {
 public:
  /**
   * @brief Default constructor. Sufficient for dynamic damping.
   */
  Damper();

  /**
   * @brief Constructor for static damping.
   *
   * @param dampingFactor Damping factor that determines the ratio of previous and current Fock matrix that is being
   * mixed in.
   */
  Damper(const double dampingFactor);

  /**
   * @brief           Constructor for arithmetic series damping with varying damping factor in a given range.
   *
   *                  Damping with varying damping factor in range dStart to dEnd. Note that
   *                  the end value must not be reached. The damping will stop when the damping
   *                  factor falls below dEnd. Negative damping factors are not possible.
   *                  If iStartUp > 0, it will perform iStartUp static damping steps at the
   *                  beginning of the SCF calculation.  It can make sense to use a higher
   *                  damping factor multiple times at the beginning because the Fock matrix
   *                  varies strongly in the first few iterations.
   *
   * @param dStart    Start damping factor (dStart > 0.0)
   * @param dStep     Step width (dStep > 0.0)
   * @param dEnd      End damping factor (dEnd > 0.0)
   * @param iStartUp  Number of static damping steps with dStart at the beginning of the SCF
   *
   */
  Damper(const double dStart, const double dStep, const double dEnd, int iStartUp);

  /**
   * @brief Default destructor.
   */
  virtual ~Damper() = default;

  /**
   * @param newFock The new Fock matrix to be damped.
   */
  void staticDamp(FockMatrix<SCFMode>& newFock);

  /**
   * @param newFock The new Fock matrix to be damped.
   */
  void arithmeticSeriesDamp(FockMatrix<SCFMode>& newFock);

  /**
   * @param newFock The new Fock matrix to be damped.
   * @param newDensity The new Density matrix to be damped.
   */
  void dynamicDamp(FockMatrix<SCFMode>& newFock, DensityMatrix<SCFMode> newDensity);

 private:
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _oldFock;
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _oldDensity;
  double _dampingFactor = 0.0;
  bool _initialized;
  double _dStart = 0.7;
  double _dStep = 0.05;
  double _dEnd = 0.2;
  int _iStartUp = 2;
};

} /* namespace Serenity */

#endif /* MATHS_DAMPER_DAMPER_H_ */
