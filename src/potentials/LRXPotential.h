/**
 * @file LRXPotential.h
 *
 * @date Mar 31, 2017
 * @author M. Boeckers
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

#ifndef POTENTIALS_LRXPOTENTIAL_H_
#define POTENTIALS_LRXPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "settings/Options.h"
#include "potentials/Potential.h"


namespace Serenity {
/**
 * @class LRXPotential LRXPotential.h
 * @brief Calculates the matrix representation of the long range contribution to the
 *        exchange potential in range-separated exchange--correlation functionals.
 */
template<Options::SCF_MODES SCFMode>
class LRXPotential : public Potential<SCFMode>,
                     public ObjectSensitiveClass<Basis>,
                     public ObjectSensitiveClass<DensityMatrix<SCFMode> >{
public:
	/**
	 * @brief Constructor.
	 * @param dMat The density matrix controller.
	 * @param exchangeRatio The exchange ratio.
	 * @param prescreeningThreshold The schwartz prescreening threshold.
	 * @param mu The range seperation parameter.
	 */
  LRXPotential(
      std::shared_ptr<DensityMatrixController<SCFMode> > dMat,
      const double exchangeRatio,
      const double prescreeningThreshold,
      const double mu);

  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Adds an increment to an existing potential in matrix form.
   *
   * This function will prescreen based on the delta-density matrix.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP An increment of the density matrix.
   */
  void addToMatrix(
      FockMatrix<SCFMode>& F,
      const DensityMatrix<SCFMode>& deltaP);

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;


  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Potential is linked to the basis it is defined in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final{
    _potential.reset(nullptr);
  };

  virtual ~LRXPotential() = default;

private:
  ///@brief Threshold for the integral prescreening.
  const double _prescreeningThreshold;
  ///@brief The exchange ratio.
  const double _exc;
  ///@brief Density matrix controller for this potential
  std::shared_ptr<DensityMatrixController<SCFMode> > _dMatController;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode> >_potential;
  ///@brief The range separation parameter
  const double _mu;
};

} /* namespace Serenity */

#endif /* POTENTIALS_LRXPOTENTIAL_H_ */
