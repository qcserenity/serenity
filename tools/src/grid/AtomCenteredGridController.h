/**
 * @file   AtomCenteredGridController.h
 *
 * @date   Mar 09, 2016
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
#ifndef ATOMCENTEREDGRIDCONTROLLER_H_
#define ATOMCENTEREDGRIDCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "grid/GridController.h"
#include "notification/ObjectSensitiveClass.h"

namespace Serenity {
/* Forward declarations */
class Geometry;
class GridFactory;
class Atom;
class AtomCenteredGrid;
/**
 * @class GridController GridController.h
 * @brief Atom-centered version of a GridController
 *
 * The grid is guaranteed to be up to date if the notification system is properly used. As a
 * consequence any call may lead to the construction of a new integration grid (if the held grid
 * is not up to date any more because an atom was moved).
 */
class AtomCenteredGridController : public GridController, ObjectSensitiveClass<Geometry> {
 public:
  /**
   * @brief Creates a grid with the given tools which it will then manage.
   * @param geometry The geometry of the System the grid shall be constructed for
   * @param gridFactory The GridFactory producing the AtomCenteredGrid
   * @param gridPointSorting Specifies whether the HilbertRTreeSorting of the grid points shall be performed or not.
   */
  AtomCenteredGridController(std::shared_ptr<const Geometry> geometry, std::shared_ptr<GridFactory> gridFactory,
                             const bool gridPointSorting);

  virtual ~AtomCenteredGridController() = default;

  /**
   * @brief Constructs a new grid that is not sorted by position in space, but by
   *          parent atom.
   * @return An atom centered grid including the information which point belongs to which atom.
   */
  std::unique_ptr<AtomCenteredGrid> getAtomGrid();

  void notify() override final;

 private:
  void produceGrid() override final;

  const std::shared_ptr<const Geometry> _geometry;
  const std::shared_ptr<GridFactory> _gridFactory;
  const bool _gridPointSorting;
};

} // namespace Serenity
#endif /* ATOMCENTEREDGRIDCONTROLLER_H_ */
