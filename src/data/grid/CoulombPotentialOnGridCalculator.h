/**
 * @file CoulombPotentialOnGridCalculator.h
 *
 * @date Mar 31, 2016
 * @author David Schnieders
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

#ifndef BASICS_GRID_DATAONGRID_COULOMBPOTENTIALONGRIDCALCULATOR_H_
#define BASICS_GRID_DATAONGRID_COULOMBPOTENTIALONGRIDCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "data/grid/GridPotential.h"
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {

class Atom;
class BasisFunctionOnGridController;
class GridController;
template<Options::SCF_MODES SCFMode>
class DensityOnGridController;

/**
 * @class CoulombPotentialOnGridCalculator CoulombPotentialOnGridCalculator.h
 * This class calculates a real-space representation of the Coulomb potentials
 * (electron-nuclei and electron-electron). The electron-electron part is
 * calculated by "hacking" the one-electron integrals for the electron-nuclei
 * interaction: At every point r in real space, the Coulomb potential can be
 * described as
 *
 * \f$\int \sum_i^{n_\text{occ}} \frac{|\phi_i(r')|^2}{|r-r'|} \text{d}r'\f$
 *
 * The key part, i.e.
 *
 * \f$\int  \frac{|\phi_i(r')|^2}{|r-r'|} \text{d}r'\f$
 *
 * or, using a different notation
 *
 * \f$\langle \phi_i(r')|\frac{1}{r-r'}|\phi_i(r')\rangle \f$
 *
 * can be calculated using the one-electron integrals for the electron-nuclei
 * interaction and using a charge of -1 at position r. This has to be performed
 * on every point r.
 * This may be computationally costly, but circumvents the numerical problems
 * arising when trying to evaluate \f$\frac{|\phi_i(r')|^2}{|r-r'|}\f$ for r=r'.
 */
class CoulombPotentialOnGridCalculator {
 private:
  /**
   * @brief Constructor. Purely static.
   */
  CoulombPotentialOnGridCalculator() = default;

 public:
  virtual ~CoulombPotentialOnGridCalculator() = default;

  /**
   * @brief Calculates the real-space representation of the CLASSICAL electron-electron
   *        potential and adds it to a GridPotential.
   * @param result The GridPotential to which the Coulomb potential will be added.
   * @param densMat The DensityMatrix of the MOs to be used.
   */
  template<Options::SCF_MODES SCF_MODE>
  static void calculateElectronElectron(GridPotential<RESTRICTED>& result, const DensityMatrix<SCF_MODE>& densMat);
  /**
   * @brief Calculates the real-space representation of the electron-nuclei
   *        potential and adds it to a GridPotential.
   * @param result The GridPotential to which the Coulomb potential will be added.
   * @param atoms The atoms of the system.
   *
   * Even though the spin polarization is not needed (since the potentials will be
   * the same, of course) it is used to allow for an easy summation with the
   * electron-electron potential.
   */
  static void calculateElectronNuclei(GridPotential<RESTRICTED>& result, const std::vector<std::shared_ptr<Atom>>& atoms);
};

} /* namespace Serenity */

#endif /* BASICS_GRID_DATAONGRID_COULOMBPOTENTIALONGRIDCALCULATOR_H_ */
