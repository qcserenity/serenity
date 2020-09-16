/**
 * @file   DispersionCorrectionCalculator.h
 *
 * @date   Nov 26, 2015
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

#ifndef DISPERSIONCORRECTIONCALCULATOR_H_
#define DISPERSIONCORRECTIONCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "dft/dispersionCorrection/DispersionData.h"
#include "geometry/Atom.h"
#include "geometry/Geometry.h"
#include "math/Matrix.h"
#include "settings/DFTOptions.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/**
 * @class DispersionCorrectionCalculator DispersionCorrectionCalculator.h
 * @brief A class for the calculation of Stefan Grimme's dispersion correction.
 *
 * This class includes static functions that can calculate the DFT-D3(X) dispersion correction.
 * All functions include empty return values in case the DFT_DISPERSION_CORRECTIONS enum is set to NONE.
 *
 * For further references on the QC behind this class see:
 * D3: S.Grimme, J.Antony, S.Ehrlich and H.Krieg, J.Chem.Phys., 132, (2010), 154104
 *
 * Further references:
 * D1: S.Grimme, J.Comput.Chem., 25, (2004), 1463-1476
 * D2: S.Grimme, J.Comput.Chem., 27, (2006), 1787-1799
 *
 */
class DispersionCorrectionCalculator {
 public:
  DispersionCorrectionCalculator() = default;
  virtual ~DispersionCorrectionCalculator() = default;

  /**
   * @brief Analytically calculates the DFT dispersion correction to the energy.
   *
   * @param dispType   Which correction to use.
   * @param geometry   The geometry.
   * @param functional The functional for which the correction shall be calculated.
   * @return The DFT dispersion correction to the energy in a.u.
   */
  static double calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS dispType,
                                               std::shared_ptr<const Geometry> geometry,
                                               const CompositeFunctionals::XCFUNCTIONALS functional);
  /**
   * @brief Calculates the DFT dispersion correction to the interaction energy between two geometries.
   *
   * @param dispType            Which correction to use.
   * @param activeGeometry      The active geometry.
   * @param environmentGeometry The other geometry.
   * @param functional          The exchange--correlation functional.
   */
  static double calcDispersionEnergyInteractionCorrection(Options::DFT_DISPERSION_CORRECTIONS dispType,
                                                          std::shared_ptr<const Geometry> activeGeometry,
                                                          std::shared_ptr<const Geometry> environmentGeometry,
                                                          const CompositeFunctionals::XCFUNCTIONALS functional);
  /**
   * @brief Analytically calculates the DFT dispersion correction to the energy.
   *
   * @param geometry   The geometry.
   * @param functional The functional for which the correction shall be calculated.
   * @return The DFT dispersion correction to the energy in a.u.
   */
  template<Options::DFT_DISPERSION_CORRECTIONS>
  static double calcDispersionEnergyCorrection(std::shared_ptr<const Geometry> geometry,
                                               const CompositeFunctionals::XCFUNCTIONALS functional);
  /**
   * @brief Calculates the DFT dispersion correction to the interaction energy between two geometries.
   *
   * @param activeGeometry The active geometry.
   * @param environmentGeometry The other geometry.
   * @param functional The exchange--correlation functional.
   */
  template<Options::DFT_DISPERSION_CORRECTIONS>
  static double calcDispersionEnergyInteractionCorrection(std::shared_ptr<const Geometry> activeGeometry,
                                                          std::shared_ptr<const Geometry> environmentGeometry,
                                                          const CompositeFunctionals::XCFUNCTIONALS functional);

  /**
   * @brief Analytically calculates the DFT dispersion correction to the geometry gradient, adds it to the geometry.
   *
   * The scaling is N^2 in atoms.
   *
   * @param dispType   Which correction to use.
   * @param geometry   The geometry to which the correction shall  be added.
   * @param functional The functional for which the correction shall be calculated.
   * @return The DFT dispersion correction to the gradient (natoms,3).
   */
  static Eigen::MatrixXd calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS dispType,
                                                          std::shared_ptr<const Geometry> geometry,
                                                          const CompositeFunctionals::XCFUNCTIONALS functional);

  /**
   * @brief Analytically calculates the DFT dispersion correction to the geometry gradient.
   *
   * The scaling is N^2 in atoms.
   *
   * @param geometry   The geometry to which the correction shall  be added.
   * @param functional The functional for which the correction shall be calculated.
   * @return The DFT dispersion correction to the gradient (natoms,3).
   */
  template<Options::DFT_DISPERSION_CORRECTIONS>
  static Eigen::MatrixXd calcDispersionGradientCorrection(std::shared_ptr<const Geometry> geometry,
                                                          const CompositeFunctionals::XCFUNCTIONALS functional);

  /**
   * @brief Numerically calculates the DFT dispersion correction to the Hessian.
   *
   * @param dispType   Which correction to use.
   * @param geometry   The geometry.
   * @param functional The functional for which the correction shall be calculated.
   * @return The DFT dispersion correction to the Hessian.
   */
  static std::vector<Eigen::MatrixXd> calcDispersionHessianCorrection(Options::DFT_DISPERSION_CORRECTIONS dispType,
                                                                      std::shared_ptr<const Geometry> geometry,
                                                                      const CompositeFunctionals::XCFUNCTIONALS functional);

  /**
   * @brief Numerically calculates the DFT dispersion correction to the Hessian.
   *
   *
   * @param geometry   The geometry.
   * @param functional The functional for which the correction shall be calculated.
   * @return The DFT dispersion correction to the Hessian.
   */
  template<Options::DFT_DISPERSION_CORRECTIONS>
  static std::vector<Eigen::MatrixXd> calcDispersionHessianCorrection(std::shared_ptr<const Geometry> geometry,
                                                                      const CompositeFunctionals::XCFUNCTIONALS functional);

 private:
  /**
   * @brief Returns the C6 parameter for one atom pair.
   * @param atomI   The first atom.
   * @param atomJ   The second atom.
   * @param nCoordI The coordination number of the first atom.
   * @param nCoordI The coordination number of the second atom.
   * @return The C6 parameter.
   */
  static double getC6(std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ, const double& nCoordI, const double& nCoordJ);

  /**
   * @brief Returns the gradient of the C6 parameter for one atom pair.
   *
   * Derivative of C6 w.r.t. coordination number dC_6/d_nCoordI, dC_6/d_nCoordJ
   *
   * @param atomI   The first atom.
   * @param atomJ   The second atom.
   * @param nCoordI The coordination number of the first atom.
   * @param nCoordJ The coordination number of the second atom.
   * @return The C6 parameter derivatives against movement of both atoms
   */
  static std::pair<double, double> getDeltaC6(std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ, double& nCoordI,
                                              double& nCoordJ);

  /**
   * @brief Calculates the coordination numbers for a given geometry.
   * @param geometry The geometry.
   * @return A vector containing the coordination number for each atom in the given Geometry.
   */
  static std::vector<double> calcCoordNumbers(std::shared_ptr<const Geometry> geometry,
                                              std::shared_ptr<const Geometry> environmentGeometry = nullptr);

  static void calculateD3Term(std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ, const double& coordI,
                              const double& coordJ, const double& rs6, const double& rs18, const double& alp,
                              double& e6, double& e8);

  static void calculateD3BJTerm(std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ, const double& coordI,
                                const double& coordJ, const double& rs6, const double& rs18, double& e6, double& e8);
};

} /* namespace Serenity */

#endif /* DISPERSIONCORRECTIONCALCULATOR_H_ */
