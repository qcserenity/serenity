/**
 * @file EnergyComponentPrinter.h
 *
 * @date Jan 5, 2015
 * @author Thomas Dresselhaus
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
#ifndef ENERGYCOMPONENTPRINTER_H
#define ENERGYCOMPONENTPRINTER_H
/* Include Std and External Headers */
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Serenity {
/* Forward declarations */
/**
 * These types of energy can appear in Serenity. They are related to each other via the
 * ENERGY_CONTRIBUTIONS_CHILDREN_MAP and ENERGY_CONTRIBUTIONS_PARENTS_MAP. For each entry these
 * relations are clearly defined. If you need a different combination of energy values, you may
 * think about extending this list (and the maps accordingly).
 */
enum class ENERGY_CONTRIBUTIONS;

typedef std::vector<ENERGY_CONTRIBUTIONS> EnergyContributionChildren;
typedef std::vector<ENERGY_CONTRIBUTIONS> EnergyContributionParents;

/**
 * @class EnergyComponentPrinter EnergyComponentPrinter.h
 * @brief Manages how energy contributions are printed. Purely static, cannot be instanciated.
 */
class EnergyComponentPrinter {
 private:
  EnergyComponentPrinter() = default;

 public:
  virtual ~EnergyComponentPrinter() = default;
  /**
   * @brief Properly prints a certain energy contribution
   * @param dataToPrint
   */
  static void printEnergyComponent(const std::pair<ENERGY_CONTRIBUTIONS, double>& dataToPrint);

  // The redundancy is a bit annoying, but one cannot trivially loop over an enum, so we introduce a
  // vector as well.
  static const std::vector<ENERGY_CONTRIBUTIONS> ENERGY_CONTRIBUTIONS_VECTOR;
  /**
   * This map holds the relations of ENERGY_CONTRIBUTIONS to each other into one direction. A
   * non-empty vector for an entry means the contribution can be split up into smaller contributions.
   */
  static const std::map<ENERGY_CONTRIBUTIONS, EnergyContributionChildren> ENERGY_CONTRIBUTIONS_CHILDREN_MAP;

  /**
   * @returns the ENERGY_CONTRIBUTIONS_PARENTS_MAP, never call but use that object straight away.
   */
  static const std::map<ENERGY_CONTRIBUTIONS, EnergyContributionParents> determineEnergyContributionParentMap();
  /**
   * Analog to ENERGY_CONTRIBUTIONS_CHILDREN_MAP. A non-empty vector for an entry means the
   * contribution is part of another 'higher-level' contribution.
   */
  static const std::map<ENERGY_CONTRIBUTIONS, EnergyContributionParents> ENERGY_CONTRIBUTIONS_PARENTS_MAP;
#ifndef NDEBUG
  /**
   * @returns true iff the _energyContributionMap is complete.
   */
  static bool checkEnergyContributionMap() {
    for (const auto& energyType : ENERGY_CONTRIBUTIONS_VECTOR) {
      if (_energyContributionMap.find(energyType) == _energyContributionMap.end())
        return false;
    }
    return true;
  };
  /**
   * @returns true iff the ENERGY_CONTRIBUTIONS_CHILDREN_MAP is complete.
   */
  static bool checkEnergyContributionChildrenMap() {
    for (const auto& energyType : ENERGY_CONTRIBUTIONS_VECTOR) {
      if (ENERGY_CONTRIBUTIONS_CHILDREN_MAP.find(energyType) == ENERGY_CONTRIBUTIONS_CHILDREN_MAP.end())
        return false;
    }
    return true;
  };
#endif
  /**
   * TODO move somewhere else
   * @brief checks whether descendant is a child, grandchild, great grandchild, ... of ancestor.
   *
   * @param descendant
   * @param ancestor
   * @returns true iff descendant is actually a descendant of ancestor (or they are the same)
   */
  static bool isAncestorOf(const ENERGY_CONTRIBUTIONS& descendant, const ENERGY_CONTRIBUTIONS& ancestor);

 private:
  static const std::map<ENERGY_CONTRIBUTIONS, std::string> _energyContributionMap;
};

} /* namespace Serenity */
#endif /* ENERGYCOMPONENTPRINTER_H */
