/**
 * @file EffectiveCorePotential.h
 *
 * @date May 24, 2018
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

#ifndef POTENTIALS_EFFECTIVECOREPOTENTIAL_H_
#define POTENTIALS_EFFECTIVECOREPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "geometry/Atom.h"
#include "potentials/Potential.h"

namespace Serenity {

class SystemController;

template<Options::SCF_MODES SCFMode>
class EffectiveCorePotential : public Potential<SCFMode>, public ObjectSensitiveClass<Basis> {
 public:
  /**
   * @brief Constructor
   * @param atoms The atoms with potential ECPs.
   * @param basis The basis.
   */
  EffectiveCorePotential(std::shared_ptr<SystemController> system, std::vector<std::shared_ptr<Atom>> atoms,
                         std::shared_ptr<BasisController> basis);
  /**
   * @brief Default destructor.
   */
  virtual ~EffectiveCorePotential() = default;
  /**
   * @brief Getter for the actual potential.
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;
  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;
  /**
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override final;
  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  ///@brief Flag if the the contribution has to be calculated/ ECPs are used.
  bool _notZero;
  ///@brief The system controller (gradients only).
  std::weak_ptr<SystemController> _system;
  ///@brief The atoms.
  std::vector<std::shared_ptr<Atom>> _atoms;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_EFFECTIVECOREPOTENTIAL_H_ */
