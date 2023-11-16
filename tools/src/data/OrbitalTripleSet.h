/**
 * @file OrbitalTripleSet.h
 *
 * @date Feb. 18, 2021
 * @author Moritz Bensberg
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

#ifndef SERENITY_ORBITALTRIPLESET_H
#define SERENITY_ORBITALTRIPLESET_H
/* Include Std and External Headers */
#include <Eigen/SparseCore> //Sparse matrices.
#include <memory>           //smrt_ptr.
#include <vector>           //std::vector.

namespace H5 {
class H5File;
} // namespace H5
namespace Serenity {

/* Forward Declarations */
class OrbitalTriple;
/**
 * @class
 * @brief Handles a set of orbital triples.
 *
 * This provides a convenient interface to handle a list of orbital triples simultaneously.
 */
class OrbitalTripleSet : public std::vector<std::shared_ptr<OrbitalTriple>> {
 public:
  /**
   * @brief Default constructor.
   */
  OrbitalTripleSet();
  /**
   * @brief Default destructor.
   */
  ~OrbitalTripleSet();
  /**
   * @brief Getter for the union of the fitting domain.
   * @return The union of the fitting domains.
   */
  Eigen::SparseVector<int> getTotalFittingDomain() const;
};

} /* namespace Serenity */

#endif // SERENITY_ORBITALTRIPLESET_H
