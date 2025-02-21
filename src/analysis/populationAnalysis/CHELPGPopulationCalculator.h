/**
 * @file CHELPGPopulationCalculator.h
 *
 * @date October 7, 2024
 * @author: Thorben Wiegmann
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
#ifndef ANALYSIS_CHELPGPOPULATIONCALCULATOR_H_
#define ANALYSIS_CHELPGPOPULATIONCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

/* Forward declarations */
class SystemController;
class Atom;
class GridController;

/**
 * @class CHELPGPopulationCalculator CHELPGPopulationCalculator.h
 * @brief Performs a CHELPG Population Analysis for a given system.
 *        [1] J. Comput. Chem. 1990, 3, 11, 361-373;  https://doi.org/10.1002/jcc.540110311
 */
template<Options::SCF_MODES SCFMode>
class CHELPGPopulationCalculator {
 public:
  /**
   * @brief Constructor
   * @param system The SystemController of the calculated system.
   */
  CHELPGPopulationCalculator(std::shared_ptr<SystemController> system);
  /**
   * @brief Default destructor
   */
  virtual ~CHELPGPopulationCalculator() = default;
  /**
   * @brief Calculates the atom populations (if necessary) and returns them.
   * @return The atom populations.
   */
  const Eigen::VectorXd& getAtomPopulations() {
    if (!_atomPopulations) {
      calculateCHELPGPopulations();
    }
    return *_atomPopulations;
  }

 private:
  void calculateCHELPGPopulations();
  Eigen::VectorXd calculateCharges();
  std::shared_ptr<SystemController> _system;
  std::shared_ptr<Eigen::VectorXd> _atomPopulations;
  std::vector<std::shared_ptr<Atom>> _atoms;
  unsigned int _nAtoms;
  double _headspace;
  double _pointDistance;
};
} // namespace Serenity

#endif /* ANALYSIS_CHELPGPOPULATIONCALCULATOR_H_ */