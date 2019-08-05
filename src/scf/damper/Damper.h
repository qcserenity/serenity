/**
 * @file   Damper.h
 *
 * @date   Oct 10, 2016
 * @author Jan Unsleber
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
#ifndef MATHS_DAMPER_DAMPER_H_
#define MATHS_DAMPER_DAMPER_H_
/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h"
#include "settings/Options.h"
#include "data/SpinPolarizedData.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <cassert>
#include <memory>
#include <vector>

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
 *        the Fock matrix ore the MO coefficients. If applied to the MO coefficients, one must take care of
 *        the signs (phases) and the mixing of degenerate orbitals.
 *        During the first iterations of the SCF, damping has the effect to reduce the oscillations
 *        in the matrices. In the first iterations, \f$ \lambda \f$ can be constant or determined
 *        dynamically unit the convergence falls below some predefined threshold after which
 *        \f$ \lambda \f$ is set to zero.
 *        \n
 *        Note that damping does not preserve the orthonormality of the MOs or the indempotency of the
 *        density matrix. However, this is not a severe problem since the damping switched off in later
 *        SCF iterations. Furthermore, damping may decrease the speed of convergence, but slow convergence
 *        is better than no convergence if no better technique is available.
 *        \n
 *        Other quantum chemestry programs may use the slightly different damping relation
 *        \f[
 *          M^{(n+1)\prime} = M^{(n+1)} + \lambda M^{(n)}
 *        \f]
 *        with \f$ \lambda \in \mathrm{R} \f$ .
 *        \n
 *        For further information about damping see actual implementations.
 *        \n
 *        Ref.: See e.g. H. B. Schlegel and J. J. W. McDouall, Do you have SCF Stability and Convergence
 *              Problems? in Computational Advances in Organic Chemistry: Molecular Structure and Reactivity,
 *              167-185, 1991 or some basic text books.
 *
 */
class Damper {
public:
  /**
   * @param Default constructor
   */
  Damper() = default;
  /**
   * @brief Default destructor.
   */
  virtual ~Damper() = default;

  /**
   *
   * @param newMatrix The new Matrix to be damped.
   */
  virtual void damp(SpinPolarizedData<SCFMode, Eigen::MatrixXd>& newMatrix) = 0;
  /**
   *
   * @param newMatrix The new Matrix to be damped.
   */
  virtual void damp(FockMatrix<SCFMode>& newMatrix) = 0;

};


} /* namespace Serenity */

#endif /* MATHS_DAMPER_DAMPER_H_ */
