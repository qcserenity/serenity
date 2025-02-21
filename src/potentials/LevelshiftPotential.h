/**
 * @file LevelshiftPotential.h
 *
 * @date Nov 24, 2016
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

#ifndef POTENTIALS_LEVELSHIFTPOTENTIAL_H_
#define POTENTIALS_LEVELSHIFTPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/matrices/DensityMatrixController.h"
#include "notification/ObjectSensitiveClass.h"
#include "potentials/Potential.h"
#include "settings/Options.h"

namespace Serenity {

/**
 * @class LevelshiftPotential LevelshiftPotential.h
 *
 * @brief A class for any levelshift and projection potential.
 */
template<Options::SCF_MODES SCFMode>
class LevelshiftPotential : public Potential<SCFMode>,
                            public ObjectSensitiveClass<Basis>,
                            public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param actBasis
   * @param envDensityMatrixController
   * @param levelShiftParameter The levelshift parameter in Hartree.
   */
  LevelshiftPotential(const std::shared_ptr<BasisController> actBasis,
                      std::shared_ptr<DensityMatrixController<SCFMode>> envDensityMatrixController,
                      const double levelShiftParameter);
  /**
   * @brief Default destructor.
   */
  virtual ~LevelshiftPotential() = default;

  /**
   * @brief Getter for the matrix representation of the level-shifted potential.
   *
   * This function makes use of the RI approximation.
   *
   * @return Returns the active Systems potential in matrix form.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return Throws a SerenityError, should never be called.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Getter for the energy associated with this potential.
   *
   * In this case the energy is a linear correction to the projection
   * operator.
   *
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;
  /**
   * @brief Potential is linked to the basis it is defines in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  ///@brief Environment densities.
  std::shared_ptr<DensityMatrixController<SCFMode>> _envDMatController;
  ///@brief The levelshift parameter in Hartree.
  const double _levelShiftParameter;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
};

} /* namespace Serenity */

#endif /* POTENTIALS_LEVELSHIFTPOTENTIAL_H_ */
