/**
 * @file   EnergyComponentController.cpp
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 * 
 * @date   5. Januar 2015, 14:33
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
/* Include Class Header*/
#include "energies/EnergyComponentController.h"
/* Include Serenity Internal Headers */
#include "energies/EnergyComponentPrinter.h"
#include "energies/EnergyContributions.h"
#include "io/HDF5.h"
/* Include Std and External Headers */
#include <cassert>
#include <iostream>
#include <string>
#include <vector>


namespace Serenity {

EnergyComponentController::EnergyComponentController() :
    _naddXCSub(nullptr),
    _naddKinSub(nullptr){}

void EnergyComponentController::addComponent(
      const std::pair<ENERGY_CONTRIBUTIONS,double>& newData, bool isUnique) {
  // Check if data already exists
  if (_data.find(newData.first) == _data.end()) {
    if (isUnique) {
      // Make sure the energy contribution is not implicitly there already in a different contibution.
      // 1. Check parents
      for (const auto& parent : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_PARENTS_MAP.at(newData.first)) {
        if (_data.find(parent) != _data.end())
          throw SerenityError("EnergComponentyController: Added energy component is not unique.");
      }
      // 2. Check children; will only be a partial contribution, but is still unexpected.
      for (const auto& child : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(newData.first)) {
        if (_data.find(child) != _data.end())
          throw SerenityError("EnergComponentyController: Added energy component is not unique.");
      }
    } else {
      // If the data is not unique it is expected to be only a correction. Thus it is assumed
      // that an energy value has to be present already. This restriction may be released if needed.
      throw SerenityError("EnergComponentyController: non-unique component not found.");
    }
    //assert(EnergyComponentPrinter::isAncestorOf(newData.first, _managedEnergyType));
    _data.insert(newData);
  } else if (!isUnique) {
    // The new data is just a correction. Okay.
    _data[newData.first] += newData.second;
  } else {
    // Data exists already but should be unique. This is not allowed.
    throw std::exception();
  }
}

void EnergyComponentController::replaceComponent(
    const std::pair<ENERGY_CONTRIBUTIONS, double>& newData) {
  // Will fail if an energy component is not yet present.
  _data.at(newData.first) = newData.second;
}

void EnergyComponentController::addOrReplaceComponent(
    const std::pair<ENERGY_CONTRIBUTIONS, double>& newData) {
  if (_data.find(newData.first) != _data.end()){
    _data.at(newData.first) = newData.second;
  } else {
    _data.insert(newData);
  }
}

void EnergyComponentController::addOrReplaceComponent(
    const ENERGY_CONTRIBUTIONS contr, const double newData) {
  if (_data.find(contr) != _data.end()){
    _data.at(contr) = newData;
  } else {
    _data.insert(std::make_pair(contr,newData));
  }
}

double EnergyComponentController::getEnergyComponent(const ENERGY_CONTRIBUTIONS energyType) const {
  if(_data.find(energyType) != _data.end()) {
    // Data is there as is
    return _data.at(energyType);
  } else {
    // Check whether it's a compound energy which can be constructed from known data
    return getEnergyComponentFromChildren(energyType);
  }
}

double EnergyComponentController::getTotalEnergy() const{
  double totalE = 0.0;
  for(auto const &energy : _data) {
    totalE += energy.second;
  }
  return totalE;
}

double EnergyComponentController::getEnergyComponentFromChildren(
        const ENERGY_CONTRIBUTIONS energyType) const {
  // False if the present data is insufficient to calculate the requested energy type.
  assert(EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType).size() > 0);
  double energy = 0.0;
  for (const auto& child : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType)) {
    if (_data.find(child) == _data.end()) {
      energy += getEnergyComponentFromChildren(child);
    } else {
      energy += _data.at(child);
    }
  }
  return energy;
}

bool EnergyComponentController::checkEnergyComponentFromChildren(
        const ENERGY_CONTRIBUTIONS energyType) const {
  // False if the present data is insufficient to calculate the requested energy type.
  if(not (EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType).size() > 0)){
    return false;
  }
  bool result = true;
  for (const auto& child : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType)) {
    if (_data.find(child) == _data.end()) {
      result *= checkEnergyComponentFromChildren(child);
    }
  }
  return result;
}

void EnergyComponentController::printAllComponents() const {
  // Just print data as is
  std::cout << "Energy components as stored:" << std::endl;
  for (const auto& field : _data) EnergyComponentPrinter::printEnergyComponent(field);
  std::cout << std::string(100, '-') << std::endl;
  std::cout << "Derived energy components:" << std::endl;
  if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT)){
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT);
    EnergyComponentPrinter::printEnergyComponent(std::make_pair(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT,
                                 getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_WFT_DFT)));
  } else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT)){
      printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT);
      EnergyComponentPrinter::printEnergyComponent(std::make_pair(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT,
                                   getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_WF_DFT)));
  }else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::HF_ENERGY)){
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::HF_ENERGY);
    EnergyComponentPrinter::printEnergyComponent(std::make_pair(ENERGY_CONTRIBUTIONS::HF_ENERGY,
                                 getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::HF_ENERGY)));
  }
  else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT)){
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT);
    EnergyComponentPrinter::printEnergyComponent(std::make_pair(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT,
                                 getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::FDE_SUPERSYSTEM_ENERGY_DFT_DFT)));
  } else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT)){
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT);
    EnergyComponentPrinter::printEnergyComponent(std::make_pair(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT,
                                 getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::PBE_SUPERSYSTEM_ENERGY_DFT_DFT)));
  } else if (checkEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY)){
    printChildrenIfNotStored(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY);
    EnergyComponentPrinter::printEnergyComponent(std::make_pair(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY,
                                 getEnergyComponentFromChildren(ENERGY_CONTRIBUTIONS::KS_DFT_ENERGY)));
  }
  std::cout << std::string(100, '-') << std::endl;
}

void EnergyComponentController::printChildrenIfNotStored(const ENERGY_CONTRIBUTIONS energyType) const {
  for (const auto& child : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP.at(energyType)) {
    if (_data.find(child) == _data.end()) {
      EnergyComponentPrinter::printEnergyComponent(
           std::make_pair(child, getEnergyComponent(child)));
      printChildrenIfNotStored(child);
    }
  }
}

void EnergyComponentController::printDerived(const ENERGY_CONTRIBUTIONS energyType) const {
  for (auto const &energy : EnergyComponentPrinter::ENERGY_CONTRIBUTIONS_CHILDREN_MAP){
    for (auto const child : energy.second){
      if (child == energyType){
        if (checkEnergyComponentFromChildren(energy.first)){
          EnergyComponentPrinter::printEnergyComponent(std::make_pair(energy.first,getEnergyComponentFromChildren(energy.first) ));
        }
      }
    }
  }
}

void EnergyComponentController::toHDF5(std::string fBaseName, std::string id){
  std::vector<std::string> writeData;
  for(auto component : _data){
    std::string contr=std::to_string(int(component.first));
    std::string val=std::to_string(component.second);
    writeData.push_back(contr);
    writeData.push_back(val);
  }
  std::string name = fBaseName+".h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save_std_vector(file, "Energies", writeData);
  HDF5::save_scalar_attribute(file,"ID",id);
  file.close();
}

void EnergyComponentController::fromHDF5(std::string fBaseName, std::string id){
  std::vector<std::string> readData;
  HDF5::Filepath name(fBaseName+".h5");
  HDF5::H5File file(name.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  HDF5::dataset_exists(file,"Energies");
  HDF5::attribute_exists(file,"ID");
  HDF5::check_attribute(file,"ID",id);
  HDF5::load_std_vector(file,"Energies",readData);
  file.close();
  for (unsigned int i=0; i<readData.size(); i+=2){
    _data.insert(std::make_pair(static_cast<ENERGY_CONTRIBUTIONS>(std::stoi(readData[i])),std::stod(readData[i+1])));
  }
}

} /* namespace Serenity */
