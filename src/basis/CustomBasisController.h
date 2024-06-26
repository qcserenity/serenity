/**
 * @file CustomBasisController.h
 *
 * @date Sep 11, 2015
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

#ifndef CONTROLLER_BASIS_CUSTOMBASISCONTROLLER_H_
#define CONTROLLER_BASIS_CUSTOMBASISCONTROLLER_H_

/* Include Serenity Internal Headers */
#include "basis/BasisController.h"

namespace Serenity {

class Shell;
class Basis;

class CustomBasisController : public BasisController {
 public:
  CustomBasisController(std::vector<std::shared_ptr<const Shell>>& basisFunctions, const std::string& basisString);
  virtual ~CustomBasisController() = default;

 private:
  virtual std::unique_ptr<Basis> produceBasisFunctionVector() override final;
  virtual void postConstruction() override final;
};
} // namespace Serenity
#endif /* CONTROLLER_BASIS_CUSTOMBASISCONTROLLER_H_ */
