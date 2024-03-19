/**
 * @file SCFAnalysis.h
 *
 * @date Dec 1, 2016
 * @author M. Boeckers
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

#ifndef SCF_SCFANALYSIS_H_
#define SCF_SCFANALYSIS_H_

/* Include Serenity Internal Headers */
#include "energies/EnergyComponentController.h"
#include "integrals/OneElectronIntegralController.h"

namespace Serenity {
/* Forward declaration */
class SystemController;
class GridController;

/**
 * @class SCFAnalysis SCFAnalysis.h
 * @param systemController
 * @param supersystemGrid
 */
template<Options::SCF_MODES SCFMode>
class SCFAnalysis {
 public:
  SCFAnalysis(std::vector<std::shared_ptr<SystemController>> systemController,
              std::shared_ptr<GridController> supersystemGrid = nullptr);
  virtual ~SCFAnalysis() = default;

  /**
   * @brief Evaluates the expectation value of the S2 operator.
   *        When DFT is used, the expectation value is calculated
   *        as functional of the density.
   * @param useUHForbitals
   * @return The S2 expectation value
   */
  double getS2(bool useUHForbitals = false);

  /**
   * @brief Calculates -<V>/<T>
   * @return returns -<V>/<T>
   */
  double getVirialRatio();

 private:
  std::vector<std::shared_ptr<SystemController>> _systemController;
  std::shared_ptr<GridController> _supersystemGrid;
};

} /* namespace Serenity */

#endif /* SCF_SCFANALYSIS_H_ */
