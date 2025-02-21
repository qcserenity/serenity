/**
 * @file   AtomTypeFactory.h
 *
 * @date   Mar 19, 2013
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
#ifndef ATOMTYPEFACTORY_H_
#define ATOMTYPEFACTORY_H_
/* Include Serenity Internal Headers */
#include "geometry/AtomType.h"
/* Include Std and External Headers */
#include <map>
#include <memory>
#include <string>

namespace Serenity {
/**
 * @class AtomTypeFactory AtomTypeFactory.h
 * @brief Generates @ref AtomType "AtomTypes" solely based on a label.
 *
 * The label usually the atom symbol as found in the PSE, but
 *   does include isotiopes such as D and T.
 *
 * This is a purely static class for global use.
 */
class AtomTypeFactory {
 public:
  /**
   * @brief The default destructor.
   */
  virtual ~AtomTypeFactory() = default;
  /**
   * @param   name A unique identifier for an atom type.
   * @returns Returns the atom type which corresponds to name as stored in the data object.
   */
  static std::shared_ptr<const AtomType> getAtomType(std::string name);
  /**
   * @returns Returns the only allowed instance of this class
   */
  static AtomTypeFactory& getInstance();

 private:
  /*
   * A map of all atom types created in the past, to ensure that none
   *   of them are built twice.
   */
  typedef std::map<std::string, std::shared_ptr<const AtomType>> AtomTypeMap;
  static AtomTypeMap _atomTypes;

  /*
   * The private function to actually generate the AtomTypes that
   *   are not yet present in the AtomTypeMap if they are requested.
   */
  static void generateAtomType(std::string name);

  // Private Constructor; this class must never be instanciated
  AtomTypeFactory() = default;
  // Stop the compiler generating methods to copy the object
  AtomTypeFactory(AtomTypeFactory const& copy) = delete;
  AtomTypeFactory& operator=(AtomTypeFactory const& copy) = delete;
};

} /* namespace Serenity */
#endif /* ATOMTYPEFACTORY_H_ */
