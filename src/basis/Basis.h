/**
 * @file   Basis.h
 *
 * @date   28. Juli 2013, 18:03
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
#ifndef BASIS_H
#define BASIS_H
/* Include Serenity Internal Headers */
#include "basis/Shell.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/**
 * @class Basis Basis.h
 * @brief A vector of basis functions
 *
 * which is wrapped up like this to have a handy object to use inside the calculations for e.g.
 * a whole molecule. A lot of data is pretty much coupled to the ordering of the held basis
 * functions. This meta-information can be received through the BasisController.
 */
class Basis : public std::vector<std::shared_ptr<const Shell>> {
 public:
  using std::vector<std::shared_ptr<const Shell>>::vector;
  virtual ~Basis() = default;
};

} /* namespace Serenity */
#endif /* BASIS_H */
