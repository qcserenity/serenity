/**
 * @file   EnergyComponentController.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   5. Januar 2015, 14:33
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
/* Include Class Header*/
#include "energies/EnergyComponentController.h"
/* Include Serenity Internal Headers */
#include "energies/EnergyComponentPrinter.h"
#include "energies/EnergyContributions.h"
#include "io/HDF5.h"
#include "misc/HelperFunctions.h"
/* Include Std and External Headers */
#include <cassert>
#include <fstream>
#include <iomanip> //std::setw()
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Serenity {

EnergyComponentController::EnergyComponentController() : _naddXCSub(nullptr), _naddKinSub(nullptr) {
}

void EnergyComponentController::addOrReplaceComponent(const std::pair<ENERGY_CONTRIBUTIONS, double>& newData) {
  if (_data.find(newData.first) != _data.end()) {
    _data.at(newData.first) = newData.second;
  }
  else {
    _data.insert(newData);
  }
}

void EnergyComponentController::addOrReplaceComponent(const ENERGY_CONTRIBUTIONS contr, const double newData) {
  if (_data.find(contr) != _data.end()) {
    _data.at(contr) = newData;
  }
  else {
    _data.insert(std::make_pair(contr, newData));
  }
}

double EnergyComponentController::getEnergyComponent(const ENERGY_CONTRIBUTIONS energyType) const {
  if (_data.find(energyType) != _data.end()) {
    // Data is there as is
    return _data.at(energyType);
  }
  else {
    // Check whether it's a compound energy which can be constructed from known data
    return getEnergyComponentFromChildren(energyType);
  }
}

double EnergyComponentController::getTotalEnergy() const {
  double totalE = 0.0;
  for (auto const& energy : _data) {
    totalE += energy.second;
  }
  return totalE;
}

double EnergyComponentController::getEnergyComponentFromChildren(const ENERGY_CONTRIBUTIONS energyType) const {
  if (EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.find(energyType) ==
      EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.end())
    throw SerenityError(
        (std::string) "ERROR: The energy contribution was not fully integrated into the ENERGY_CONTRIBUTION "
                      "framework.\n" +
        "       It is missing in the ENERGY_CONTRIBUTIONS_CHILDREN_MAP map. This may be an implementation error\n" +
        "       or may be caused by loading corrupted/outdated files.");
  // False if the present data is insufficient to calculate the requested energy type.
  if (not(EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType).size() > 0))
    throw SerenityError("ERROR: The energy contribution " + std::to_string((int)energyType) + " was not set before!");
  double energy = 0.0;
  for (const auto& child : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType)) {
    if (_data.find(child) == _data.end()) {
      energy += getEnergyComponentFromChildren(child);
    }
    else {
      energy += _data.at(child);
    }
  }
  return energy;
}

bool EnergyComponentController::checkEnergyComponentFromChildren(const ENERGY_CONTRIBUTIONS energyType) const {
  bool inMap = not(EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.find(energyType) ==
                   EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.end());
  if (not inMap)
    throw SerenityError(
        (std::string) "ERROR: The energy contribution was not fully integrated into the ENERGY_CONTRIBUTION "
                      "framework.\n" +
        "       It is missing in the ENERGY_CONTRIBUTIONS_CHILDREN_MAP map. This may be an implementation error\n" +
        "       or may be caused by loading corrupted/outdated files.");
  // False if the present data is insufficient to calculate the requested energy type.
  if (not(EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType).size() > 0)) {
    return false;
  }
  bool result = true;
  for (const auto& child : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType)) {
    if (_data.find(child) == _data.end()) {
      result = result && checkEnergyComponentFromChildren(child);
    }
  }
  return result;
}

void EnergyComponentController::printAllComponents() const {
  // Just print data as is
  std::cout << "Energy components as stored:" << std::endl;
  for (const auto& field : _data)
    EnergyComponentPrinter::printEnergyComponent(field);
  std::cout << std::string(100, '-') << std::endl;
  std::cout << "Derived energy components:" << std::endl;
  if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT);
    EnergyComponentPrinter::printEnergyComponent(
        std::make_pair(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT,
                       getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT)));
  }
  else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT);
    EnergyComponentPrinter::printEnergyComponent(
        std::make_pair(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT,
                       getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT)));
  }
  else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::HF_ENERGY)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::HF_ENERGY);
    EnergyComponentPrinter::printEnergyComponent(
        std::make_pair(ENERGY_CONTRIBUTIONS::HF_ENERGY, getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::HF_ENERGY)));
  }
  else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT);
    EnergyComponentPrinter::printEnergyComponent(
        std::make_pair(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT,
                       getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT)));
  }
  else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT);
    EnergyComponentPrinter::printEnergyComponent(
        std::make_pair(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT,
                       getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT)));
  }
  else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);
    EnergyComponentPrinter::printEnergyComponent(std::make_pair(
        ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY, getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY)));
  }
  if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_KS_DFT_ENERGY)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_KS_DFT_ENERGY);
    EnergyComponentPrinter::printEnergyComponent(
        std::make_pair(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_KS_DFT_ENERGY,
                       getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_KS_DFT_ENERGY)));
  }
  if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_HF_ENERGY)) {
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_HF_ENERGY);
    EnergyComponentPrinter::printEnergyComponent(
        std::make_pair(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_HF_ENERGY,
                       getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SOLV_SCALED_EMBEDDED_HF_ENERGY)));
  }
  std::cout << std::string(100, '-') << std::endl;
}

void EnergyComponentController::printChildrenIfNotStored(const ENERGY_CONTRIBUTIONS energyType) const {
  for (const auto& child : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType)) {
    if (_data.find(child) == _data.end()) {
      EnergyComponentPrinter::printEnergyComponent(std::make_pair(child, getEnergyComponent(child)));
      printChildrenIfNotStored(child);
    }
  }
}

void EnergyComponentController::toFile(std::string fBaseName, std::string id) {
  std::ofstream ofs;
  ofs.open(fBaseName, std::ofstream::out | std::ofstream::trunc);
  ofs << "System_ID " << id << std::endl;
  ofs << "ENUM_INDEX  "
      << "ENERGY" << std::endl;
  for (auto component : _data) {
    std::string contr = std::to_string(int(component.first));
    std::string val = to_string_with_precision(component.second, 10);
    ofs << std::setw(10) << contr << "  " << val << std::endl;
  }
  ofs.close();
}

void EnergyComponentController::fromFile(std::string fBaseName, std::string id) {
  std::string line;
  std::string contr;
  std::string val;
  std::ifstream input(fBaseName);
  (void)id;
  while (getline(input, line)) {
    std::istringstream iss(line);
    contr = "";
    val = "";
    iss >> contr;
    if (contr == "System_ID")
      continue;
    if (contr == "ENUM_INDEX")
      continue;
    if (contr.empty())
      continue;
    if (contr == "#")
      continue;
    iss >> val;
    _data.insert(std::make_pair(static_cast<ENERGY_CONTRIBUTIONS>(std::stoi(contr)), std::stod(val)));
  }
}

} /* namespace Serenity */
