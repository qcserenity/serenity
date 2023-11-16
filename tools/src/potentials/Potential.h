/**
 * @file Potential.h
 *
 * @date Nov 22, 2016
 * @author Jan Unsleber
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

#ifndef POTENTIALS_POTENTIAL_H_
#define POTENTIALS_POTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"

namespace Serenity {

/* Forward Declarations */
class BasisController;

/**
 * @class Potential Potential.h
 * @brief An interface for all the potentials.
 *
 * All of the potentials interfaced here are in their matrix
 *   representation thus they are specific for one basis.
 * Each of the potentials needs to implement the getMatrix()
 *   function which returns the potential in matrix form in
 *   its entirety.
 */
template<Options::SCF_MODES SCFMode>
class Potential {
 public:
  /**
   * @brief Constructor
   * @param basis The basis this potential is defined in.
   */
  Potential(std::shared_ptr<BasisController> basis) : _basis(basis){};
  /// @brief Default destructor.
  virtual ~Potential() = default;
  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() = 0;

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  virtual double getEnergy(const DensityMatrix<SCFMode>& P) = 0;

  /**
   * @brief Getter for the gradient associated with this potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */

  virtual Eigen::MatrixXd getGeomGradients() = 0;

  /**
   * @param Small function for sanity checks.
   * @param basis The basis to compare with potential basis with.
   * @return Boolean, true if the basis sets match.
   */
  virtual bool compareBasis(std::shared_ptr<BasisController> basis) {
    return _basis == basis;
  };

 protected:
  ///@brief The basis this potential is defined in.
  std::shared_ptr<BasisController> _basis;
};

} /* namespace Serenity */

#endif /* POTENTIALS_POTENTIAL_H_ */
