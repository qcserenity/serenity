/**
 * @file   Kernel.h
 *
 * @date   Mar 19, 2017
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

#ifndef POSTHF_LRSCF_KERNEL_H_
#define POSTHF_LRSCF_KERNEL_H_

/* Include Serenity Internal Headers */
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "dft/functionals/wrappers/PartialDerivatives.h"
#include "system/SystemController.h"


namespace Serenity {
/**
 * @brief Class that calculates and stores objects to calculate the
 *        exchange-correlation kernel contribution to the LRSCF Hamiltonian. \n
 *        Warning: The getter functions will return nullptr if the quantity
 *        cannot be calculated (e.g. when using a LDA functional and requesting
 *        derivatives w.r.t. the gradient of the density or when requesting non-
 *        additive objects when no FDE calculation is performed).
 */
template<Options::SCF_MODES T> class Kernel {
public:
  Kernel(
      std::shared_ptr<SystemController> activeSystem,
      std::vector<std::shared_ptr<SystemController> > environmentSystems,
      bool superSystemGrid,
      bool noNaddKernel,
      Options::FUNCTIONALS func,
      Options::FUNCTIONALS naddKinFunc,
      Options::FUNCTIONALS naddXCFunc);
  virtual ~Kernel() = default;

  /**
   *
   * @return Returns the gridController which controls the grid on which
   *         the exchange correlation kernel is defined
   */
  std::shared_ptr<GridController> getGridController();


  std::shared_ptr<d2F_dRho2<T> > getD2F_dRho2();
  std::shared_ptr<dF_dSigma<T> > getDF_dSigma();
  std::shared_ptr<d2F_dSigma2<T> > getD2F_dSigma2();
  std::shared_ptr<d2F_dRhodSigma<T> > getD2F_dRhodSigma();

  std::shared_ptr<d2F_dSigma2<T> > getTotalNaddD2F_dSigma2();
  std::shared_ptr<d2F_dRhodSigma<T> > getTotalNaddD2F_dRhodSigma();

  std::shared_ptr<Gradient<DensityOnGrid<T> > > getActiveDensityGradient();

  std::shared_ptr<Gradient<DensityOnGrid<T> > > getTotalDensityGradient();

private:
  //The system controller of the active system
  std::shared_ptr<SystemController> _activeSystem;
  //The system controller of the environment systems
  std::vector<std::shared_ptr<SystemController> > _environmentSystems;
  //If true, use super-system grid
  bool _superSystemGrid;
  //If true, neglect non-additive kernel
  bool _noNaddKernel;
  //Functionals
  Options::FUNCTIONALS _func;
  Options::FUNCTIONALS _naddKinFunc;
  Options::FUNCTIONALS _naddXCFunc;


  //True if at least one environment system is present. If true, perform FDE-TDDFT calculation.
  bool _fde;

  //Data objects:
  std::shared_ptr<d2F_dRho2<T> > _d2FdRho2;
  std::shared_ptr<dF_dSigma<T> > _dFdSigma;
  std::shared_ptr<d2F_dSigma2<T> > _d2FdSigma2;
  std::shared_ptr<d2F_dRhodSigma<T> > _d2FdRhodSigma;

  //For some of the GGA parts of the total density, it is not possible to directly add the non-additive
  //contribution to the data objects above.
  std::shared_ptr<d2F_dSigma2<T> > _totalNaddD2FdSigma2;
  std::shared_ptr<d2F_dRhodSigma<T> > _totalNaddD2FdRhodSigma;

  //Density gradients
  std::shared_ptr<Gradient<DensityOnGrid<T> > > _activeDensityGradient;
  std::shared_ptr<Gradient<DensityOnGrid<T> > > _totalDensityGradient;


  //True if gga functional is used
  bool _gga;

  //True if non-additive gga functional is used
  bool _naddGGA;

  //True if kernel objects have been calculated
  bool _hasBeenCalculated;



  //The grid controller controlling the grid on which the kernel is calculated and stored
  std::shared_ptr<GridController> _gridController;

  //Calculates the different parts of the exchange correlation kernel on a grid of points.
  void calculateKernelOnGrid();

  //Calculates the (unrestricted) density of a system on the chosen grid of points (see LRSCFTask for
  //more information about the grids) and initializes the DensityOnGridController for this density
  std::unique_ptr<DensityOnGrid<T> >
    getDensity(
        std::shared_ptr<SystemController> systemController,
        std::shared_ptr<DensityOnGridController<T> >& densityOnGridController);


};

} /* namespace Serenity */

#endif /* POSTHF_LRSCF_KERNEL_H_ */
