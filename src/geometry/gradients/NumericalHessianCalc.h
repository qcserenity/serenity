/**
 * @file NumericalHessianCalc.h
 *
 * @date Feb 23, 2017
 * @author Kevin Klahr
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
#ifndef GEOMETRY_GRADIENTS_NUMERICALHESSIANCALC_H_
#define GEOMETRY_GRADIENTS_NUMERICALHESSIANCALC_H_

/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "geometry/gradients/HessianCalculator.h"
#include "math/Matrix.h"
#include "settings/Options.h"
#include "data/OrbitalController.h"
#include "settings/Settings.h"


namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES T>class ElectronicStructureCalculatorFactory;
/**
 * @class NumericalGeomGradCalc NumericalGeomGradCalc.h
 * @brief A GeometryGradientCalculator that calculates the gradients numerically.
 */
template<Options::SCF_MODES T>class NumericalHessianCalc :
	public HessianCalculator{
public:
  /**
   * @brief Constructor.
   * @param stepsize The step size for each step in the gradient calculations (in bohr).
   */
  NumericalHessianCalc(double stepsizeGrad, double stepsizeHess, bool printToFile);
  virtual ~NumericalHessianCalc() = default;
  /**
   * @brief Calculates the gradients for the given geometry numerically.
   *
   * This function calculates the gradients of all atoms in all
   * three dimensions of space numerically.
   * Calculated geometry gradients are stored in each atom.
   *
   * @param system   The system holding the geometry of interest.
   *
   */

  void calcHessian(std::shared_ptr<SystemController> systemController);

  void calcFaTHessian(std::vector<std::shared_ptr<SystemController> > activeSystems,
		  std::vector<std::shared_ptr<SystemController> > passiveSystems,
		  Options::KINFUNCTIONALS FaTnaddKinFunc,
		  Options::XCFUNCTIONALS FaTnaddXCFunc,
		  int FatmaxCycles,
		  double FaTenergyConvThresh,
		  double FaTgridCutOff,
		  Options::DFT_DISPERSION_CORRECTIONS dispersion);

//  void calcActiveSystemHessian(std::shared_ptr<SystemController> activeSystem,
//		  std::vector<std::shared_ptr<SystemController> > passiveSystems,
//		  Options::FUNCTIONALS FaTnaddKinFunc,
//		  Options::FUNCTIONALS FaTnaddXCFunc,
//		  bool FaTexactNaddKin,
//		  double FaTgridCutOff);

  void frequencyCalculation(Matrix<double> hessian,
                                    std::shared_ptr<Geometry> geometry,
                                    const std::vector<Settings> settings);

private:
  const double _deltaGrad;
  const double _deltaHess;
  const bool _printToFile;

};

} /* namespace Serenity */



#endif /* GEOMETRY_GRADIENTS_NUMERICALHESSIANCALC_H_ */
