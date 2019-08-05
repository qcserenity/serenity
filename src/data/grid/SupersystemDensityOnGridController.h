/**
 * @file SupersystemDensityOnGridController.h
 *
 * @date Dec 30, 2014
 * @author Thomas Dresselhaus
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
#ifndef SUPERSYSTEMDENSITYONGRIDCONTROLLER_H
#define	SUPERSYSTEMDENSITYONGRIDCONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGridController.h"
#include "notification/ObjectSensitiveClass.h"


namespace Serenity {
/* Forward Declarations */

/**
 * @class SupersystemDensityOnGridController SupersystemDensityOnGridController.h
 * @brief Manages the density and derivatives on a grid for a supersystem consisting of subsystems.
 * 
 * Make sure when using it that all subsystem controllers know about the produced object by using
 * their member functions addSensitiveObject()
 */
template<Options::SCF_MODES T>class SupersystemDensityOnGridController :
          public DensityOnGridController<T>, public ObjectSensitiveClass<DensityOnGrid<T> > {
public:
  /**
   * @brief Assumes that the subsystem controllers all work on the same grid.
   * @param subsystemDensOnGridControllers the supersystem controller will take (shared) ownership
   *                                       of each of them
   */
  SupersystemDensityOnGridController(
      const std::vector<std::shared_ptr<DensityOnGridController<T> > >& subsystemDensOnGridControllers);
  virtual ~SupersystemDensityOnGridController();

  virtual const DensityOnGrid<T>& getDensityOnGrid() override final;

  virtual const Gradient<DensityOnGrid<T> >& getDensityGradientOnGrid() override final;
  
  virtual const Hessian<DensityOnGrid<T> >& getDensityHessianOnGrid() override final;

  virtual void notify() override final;

  virtual void setHighestDerivative(unsigned int newHighestDerivative) override final;
private:
  void updateData();
  
  const std::vector<std::shared_ptr<DensityOnGridController<T> > > _subsystemDensOnGridControllers;

  bool _upToDate;
};

} /* namespace Serenity */
#endif	/* SUPERSYSTEMDENSITYONGRIDCONTROLLER_H */

