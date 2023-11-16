/**
 * @file EnergyComponentController.h
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
#ifndef ENERGYCOMPONENTCONTROLLER_H
#define ENERGYCOMPONENTCONTROLLER_H
/* Include Serenity Internal Headers */
#include "energies/EnergyComponentPrinter.h"
/* Include Std and External Headers */
#include <assert.h>
#include <map>
#include <memory>

namespace Serenity {
/* Forward Declarations */
enum class ENERGY_CONTRIBUTIONS;

/**
 * @class EnergyComponentController EnergyComponentController.h
 * @brief Manages different energy contributions for a specific system/calculation/...
 */
class EnergyComponentController {
 public:
  EnergyComponentController();
  virtual ~EnergyComponentController() = default;
  /**
   * @brief Erases all stored energy values.
   */
  void clear() {
    _data.clear();
  }
  /**
   * @brief Prints all currently held data and the value of the managed energy type.
   *
   * This method will fail if no value for the managed energy type can be calculated, i.e. if the
   * data is incomplete.
   */
  void printAllComponents() const;
  /**
   * @brief adds an energy component or adds to an existing one
   * @param newData which will be added to the list of held energy values
   * @param isUnique if true: The value may not be added to an existing field as double counting
   *                 will happen. If false: The added value is e.g. just a correction, so adding up
   *                 is ok.
   */
  void addComponent(const std::pair<ENERGY_CONTRIBUTIONS, double>& newData, bool isUnique = true);
  /**
   * @brief replaces an energy component
   * @param newData New energy contribution
   */
  void replaceComponent(const std::pair<ENERGY_CONTRIBUTIONS, double>& newData);
  /**
   * @param   energyType the kind of energy for which the value is requested
   * @returns the energy value of energyType
   */
  double getEnergyComponent(ENERGY_CONTRIBUTIONS energyType) const;
  /**
   * @returns the summed energy value of all energy components controlled by the current instance
   */
  double getTotalEnergy() const;
  /**
   * @brief adds an energy component or adds to an existing one
   * @param newData which will be added or replaced in the list of held energy values
   */
  void addOrReplaceComponent(const std::pair<ENERGY_CONTRIBUTIONS, double>& newData);
  /**
   * @brief adds an energy component or adds to an existing one
   * @param contr Energy type.
   * @param newData Energy value.
   */
  void addOrReplaceComponent(const ENERGY_CONTRIBUTIONS contr, const double newData);
  /**
   * @brief Checks if all the data is there for the reference energy.
   * @param energyType Reference energy.
   */
  bool checkEnergyComponentFromChildren(const ENERGY_CONTRIBUTIONS energyType) const;

  /**
   * @brief Writes the data to file.
   * @param fBaseName The filename with path.
   * @param id The identifier corresponding to the system.
   */
  void toFile(std::string fBaseName, std::string id);

  /**
   * @brief Reads data from file.
   * @param fBaseName The filename with path.
   * @param id The identifier corresponding to the system.
   */
  void fromFile(std::string fBaseName, std::string id);

 private:
  /*
   * @brief calculates an energy component that is not directly available
   *
   * , but can be obtained through the 'children' of this energy contribution
   */
  double getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS energyType) const;

  /**
   * @brief Prints all parts of the reference energies that are not stored .
   * @param energyType Reference energy.
   */
  void printChildrenIfNotStored(ENERGY_CONTRIBUTIONS energyType) const;

  /**
   * @brief Prints all energies available that are in part composed of the reference energy.
   * @param energyType Reference energy.
   */
  void printDerived(const ENERGY_CONTRIBUTIONS energyType) const;

  std::map<ENERGY_CONTRIBUTIONS, double> _data;

  /* ===========================================
   *   Special Fields for FaT/FDE Calculations
   * ===========================================*/
  std::unique_ptr<double> _naddXCSub;
  std::unique_ptr<double> _naddKinSub;

 public:
  /**
   * @brief Returns the kinetic energy calculated for this system on a supersystem grid.
   *
   * Note: This energy is not necessarily present, please use checkNAddKinSub() to check first.
   *       This function is present to reduce the number of evaluation of the subsystem part of
   *       the non additive kinetic energy functional in the supersystem grid
   *
   * @return The kinetic energy of this subsystem generated on a supersystem grid.
   */
  double getNAddKinSub() {
    assert((_naddKinSub != nullptr) && "No _naddKinSub set in the EnergyComponentController.");
    return *_naddKinSub;
  }
  /**
   * @brief Returns the XC energy calculated for this system on a supersystem grid.
   *
   * Note: This energy is not necessarily present, please use checkNAddXCSub() to check first.
   *       This function is present to reduce the number of evaluation of the subsystem part of
   *       the non additive XC energy functional in the supersystem grid
   *
   * @return The XC energy of this subsystem generated on a supersystem grid.
   */
  double getNAddXCSub() {
    assert((_naddXCSub != nullptr) && "No _naddKinSub set in the EnergyComponentController.");
    return *_naddXCSub;
  }

  /**
   * @brief Checks if kinetic energy of this subsystem generated on a supersystem grid is present.
   * @return Returns true is the energy is available.
   */
  bool checkNAddKinSub() {
    return (_naddKinSub != nullptr);
  };
  /**
   * @brief Checks if XC energy of this subsystem generated on a supersystem grid is present.
   * @return Returns true is the energy is available.
   */
  bool checkNAddXCSub() {
    return (_naddXCSub != nullptr);
  };

  /**
   * @brief Sets the kinetic energy of this subsystem generated on a supersystem grid.
   * @param val The energy.
   */
  void setNAddKinSub(double val) {
    this->_naddKinSub.reset(new double(val));
  }
  /**
   * @brief Sets the XC energy of this subsystem generated on a supersystem grid.
   * @param val The energy.
   */
  void setNAddXCSub(double val) {
    this->_naddXCSub.reset(new double(val));
  }
  /**
   * @brief Clears the kinetic energy of this subsystem generated on a supersystem grid.
   */
  void clearNAddKinSub() {
    this->_naddKinSub.reset(nullptr);
  }
  /**
   * @brief Clears the XC energy of this subsystem generated on a supersystem grid.
   */
  void clearNAddXCSub() {
    this->_naddXCSub.reset(nullptr);
  }
};

} /* namespace Serenity */
#endif /* ENERGYCOMPONENTCONTROLLER_H */
