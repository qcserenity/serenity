/**
 * @file ZeroPotential.h
 *
 * @date Nov 24, 2016
 * @author: Jan Unsleber
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

#ifndef POTENTIALS_ZEROPOTENTIAL_H_
#define POTENTIALS_ZEROPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "potentials/Potential.h"
#include "settings/Options.h"

namespace Serenity {
/**
 * @class ZeroPotential ZeroPotential.h
 *
 * @brief A potential dummy, with a potential of zero everywhere
 */
template<Options::SCF_MODES SCFMode>
class ZeroPotential : public Potential<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param basis The basis this potential is defined in.
   */
  ZeroPotential(std::shared_ptr<BasisController> basis) : Potential<SCFMode>(basis){};
  /// @brief Default destructor.
  virtual ~ZeroPotential() = default;
  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns a potential that is zero.
   */
  FockMatrix<SCFMode>& getMatrix() override final {
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin.setZero();
    };
    return *_potential;
  }

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override final {
    auto bascont = std::dynamic_pointer_cast<AtomCenteredBasisController>(this->_basis);
    auto natoms = bascont->getBasisIndices().size();
    Eigen::MatrixXd gradientContr(natoms, 3);
    gradientContr.setZero();
    return gradientContr;
  }

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
  virtual double getEnergy(const DensityMatrix<SCFMode>& P) override final {
    return 0.0;
  }
#pragma GCC diagnostic pop
 private:
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_ZEROPOTENTIAL_H_ */
