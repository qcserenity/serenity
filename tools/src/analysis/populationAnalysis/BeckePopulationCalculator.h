/**
 * @file BeckePopulationCalculator.h
 *
 * @date Jun 16, 2020
 * @author: Patrick Eschenbach
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
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef POSTSCF_ANALYSIS_BECKEPOPULATIONCALCULATOR_H_
#define POSTSCF_ANALYSIS_BECKEPOPULATIONCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "settings/ElectronicStructureOptions.h"

namespace Serenity {
/* Forward declarations */
class SystemController;
/**
 * @class  BeckePopulationCalculator BeckePopulationCalculator.h
 * @brief  Performs a Becke Population Analysis for a given system
 *         [1] J. Chem. Phys. 88, 2547 (1988); https://doi.org/10.1063/1.454033\n
 *         [2] J. Comput. Chem. 2013, 34, 1819–1827; https://doi.org/10.1002/jcc.23323
 */
template<Options::SCF_MODES SCFMode>
class BeckePopulationCalculator {
 public:
  /**
   * @brief Constructor
   * @param systemController: The system we want to calculate populations for
   * @param density: The density matrix of the system
   */
  BeckePopulationCalculator(std::shared_ptr<SystemController> systemController, std::shared_ptr<DensityMatrix<SCFMode>> density);
  /**
   * @brief Default destructor.
   */
  virtual ~BeckePopulationCalculator() = default;
  /**
   * @brief Calculates the atom populations (if necessary) and returns them
   * @return The atom populations.
   */
  const Eigen::VectorXd& getAtomPopulations() {
    if (!_atomPopulations) {
      calculateBeckeAtomPopulations();
    }
    return *_atomPopulations;
  }
  /**
   * @brief Calculates the atom spin populations (if necessary) and returns them
   * @return The atom populations.
   */
  const Eigen::VectorXd& getSpinPopulations() {
    if (!_spinPopulations) {
      calculateBeckeSpinPopulations();
    }
    return *_spinPopulations;
  }

 private:
  /// @brief Calculates atomic populations
  void calculateBeckeAtomPopulations();
  /// @brief Calculates atomic spin populations
  void calculateBeckeSpinPopulations();
  /// @brief The system
  std::shared_ptr<SystemController> _system;
  /// @brief The density matrix of the system
  std::shared_ptr<DensityMatrix<SCFMode>> _density;
  /// @ brief A vector of atompopulations
  std::shared_ptr<Eigen::VectorXd> _atomPopulations;
  /// @ brief A vector of spinpopulations
  std::shared_ptr<Eigen::VectorXd> _spinPopulations;
};
} // namespace Serenity
#endif /* POSTSCF_ANALYSIS_BECKEPOPULATIONCALCULATOR_H_ */
