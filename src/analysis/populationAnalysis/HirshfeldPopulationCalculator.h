/**
 * @file HirshfeldPopulationCalculator.h
 *
 * @date Mar 11, 2016
 * @author: Moritz Bensberg
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

#ifndef POSTSCF_ANALYSIS_HIRSHFELDPOPULATIONCALCULATOR_H_
#define POSTSCF_ANALYSIS_HIRSHFELDPOPULATIONCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGrid.h"
#include "data/matrices/DensityMatrix.h"
/* Include Std and External Headers */

namespace Serenity {
/* Forward declarations */

template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
template<Options::SCF_MODES SCFMode>
class DensityOnGridController;

class SystemController;
/**
 * @class HirshfeldPopulationCalculator HirshfeldPopulationCalculator.h
 * @brief A class for the calculation of the Hirshfeld population.
 *
 * The Hirshfeld population analysis is a grid-based tool. Based on atomic
 * reference densities the grid points are weighted and the molecular density
 * of each point is assigned to certain atoms based on these weights. After
 * grid integration the atomic population is obtained.
 * Reference: F. L. Hirshfeld, Theoret. Chim. Acta 44 (1977), 129-138.
 * "Bonded-Atom Fragments for Describing Molecular Charge Densities"\n\n
 *
 * The Hirshfeld population on atom \f$ A \f$ is given by\n
 *    \f$ Q_A = \int \mathrm{d}r \rho_A^{b.a.}(r) \f$\n ,
 * where \f$ \rho_A^{b.a.}(r) \f$  is the charge density of the bonded atom defined as\n
 *    \f$ \rho_A^{b.a.}(r) = w_A(r)\rho^\mathrm{mol}(r)
 *                         = \frac{\rho^\mathrm{at}_A(r)}{\rho^\mathrm{pro}(r)}\rho^\mathrm{mol}(r) \f$ \n
 * and \f$ \rho^\mathrm{at}_A(r) \f$, \f$ \rho^\mathrm{pro}(r) \f$ and \f$ \rho^\mathrm{mol}(r) \f$
 * denote the atomic density of \f$ A \f$, the promolecular density of the molecule (as a superposition
 * of atomic densities) and the converged molecular density, respectively.
 *
 * The promolecular density is generated from spherically averaged atom densities, obtained from
 * scf calculations using the same method as specified for the target system. These scf calculations will
 * be done in-place using quite weak convergence criteria. This may lead to some numerical instabilities
 * in the resulting atom populations (in the 4th decimal). If a higher accuracy is needed, the atom scf
 * convergence criteria need to be tightened (take a look at the AtomDensityGuessCalculator for this).
 *
 * In case of an UNRESTRICTED calculation, the Hirshfeld populations are evaluated separately for alpha and
 * beta densities using the same weighting factor \f$ \frac{\rho^\mathrm{at}_A(r)}{\rho^\mathrm{pro}(r)} \f$.
 */
template<Options::SCF_MODES SCFMode>
class HirshfeldPopulationCalculator : public ObjectSensitiveClass<DensityOnGrid<SCFMode>> {
 public:
  /**
   * @brief Constructor
   * @param systemController The SystemController of the system to be calculated.
   * @param densOnGrid The electron density on a grid.
   */
  HirshfeldPopulationCalculator(std::shared_ptr<SystemController> systemController,
                                std::shared_ptr<DensityOnGridController<SCFMode>> densOnGrid);
  ///@brief Default destructor.
  virtual ~HirshfeldPopulationCalculator() = default;

  /**
   * @brief Calculates the atom populations (if necessary) and returns them.
   * @return The atom populations.
   */
  const SpinPolarizedData<SCFMode, Eigen::VectorXd>& getAtomPopulations() {
    if (!_atomPopulations) {
      calculateHirshfeldAtomPopulations();
    }
    return *_atomPopulations;
  }

  /**
   * @brief The data is connected to the density of the system.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _atomPopulations.reset(nullptr);
  }

 private:
  /*
   * The actual calculation is implemented here.
   */
  void calculateHirshfeldAtomPopulations();

  std::shared_ptr<SystemController> _system;
  std::shared_ptr<DensityOnGridController<SCFMode>> _densOnGrid;
  std::shared_ptr<BasisController> _basis;
  std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>> _atomPopulations;
};
} // namespace Serenity
#endif /* POSTSCF_ANALYSIS_HIRSHFELDPOPULATIONCALCULATOR_H_ */
