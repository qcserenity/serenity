/**
 * @file   Settings.cpp
 *
 * @date   Mar 4, 2020
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

/* Include Class Header*/
#include "settings/Settings.h"

namespace Serenity {
void Settings::printSettings() {
  std::string field;
  std::string value;
  // Create folder if it does not exist
  if (not makePath(path))
    throw SerenityError((std::string) "Failed to create directory: " + path + "\n" +
                        "The directory does not already exist. A file with the same name may be\n" +
                        "preventing the directory creation.");
  std::ofstream ofs;
  ofs.open((*this).path + (*this).name + ".settings", std::ofstream::out | std::ofstream::trunc);
  ofs << "#=======================================================" << std::endl;
  ofs << "# NOTE:" << std::endl;
  ofs << "# This file contains the entire list of settings" << std::endl;
  ofs << "#  that were stored in one instance of the Settings" << std::endl;
  ofs << "#  class." << std::endl;
  ofs << "# The settings given here have not necessarily been" << std::endl;
  ofs << "#  used, e.g. a functional might not have been needed." << std::endl;
  ofs << "#  in a HF calculations." << std::endl;
  ofs << "#=======================================================" << std::endl;
  ofs << "+system" << std::endl;
  print_visitor visitor(field, value, ofs);
  visit_each((*this), visitor);
  ofs << "+dft" << std::endl;
  visit_each((*this).dft, visitor);
  ofs << "-dft" << std::endl;
  ofs << "+customfunc" << std::endl;
  visit_each((*this).customFunc, visitor);
  ofs << "-customfunc" << std::endl;
  ofs << "+scf" << std::endl;
  visit_each((*this).scf, visitor);
  ofs << "-scf" << std::endl;
  ofs << "+basis" << std::endl;
  visit_each((*this).basis, visitor);
  ofs << "-basis" << std::endl;
  ofs << "+grid" << std::endl;
  visit_each((*this).grid, visitor);
  ofs << "-grid" << std::endl;
  ofs << "+efield" << std::endl;
  visit_each((*this).efield, visitor);
  ofs << "-efield" << std::endl;
  ofs << "+pcm" << std::endl;
  visit_each((*this).pcm, visitor);
  ofs << "-pcm" << std::endl;
  ofs << "+extcharges" << std::endl;
  visit_each((*this).extCharges, visitor);
  ofs << "-extcharges" << std::endl;
  ofs << "-system" << std::endl;
  ofs.close();
}

/**
 * @brief Constructor using text input.
 */
Settings::Settings(std::ifstream& input) : Settings() {
  std::string line;
  std::string word;
  while (getline(input, line)) {
    std::istringstream iss(line);
    word = "";
    iss >> word;
    // check for comments
    if (word[0] == '#')
      continue;
    if (word.empty())
      continue;
    // blocks
    if (word[0] == '+') {
      std::string blockname = word.erase(0, 1);
      for (auto& c : blockname)
        c = std::toupper(c);
      if (blockname == "SYSTEM")
        continue;
      while (getline(input, line)) {
        std::istringstream iss2(line);
        word = "";
        iss2 >> word;
        if (word[0] == '-')
          break;
        // check for comments
        if (word[0] == '#')
          continue;
        if (word.empty())
          continue;
        std::string name = word;
        if (!(iss2 >> word)) {
          throw SerenityError("ERROR: Value missing for keyword '" + name + "'.");
        }
        std::string value = word;
        // cover possible vector inputs
        while (iss2 >> word) {
          value += " " + word;
        }
        // get rid of curly braces
        if (value.front() == '{' && value.back() == '}') {
          value = value.substr(1, value.size() - 2);
        }

        this->set(blockname, name, value);
      }
      continue;
    }
    // end of system block
    if (!word.compare("-system")) {
      /*
       * Capitalize basis label strings
       */
      for (auto& c : this->basis.label)
        c = std::toupper(c);
      for (auto& c : this->basis.auxCLabel)
        c = std::toupper(c);
      for (auto& c : this->basis.auxJLabel)
        c = std::toupper(c);
      // Set path string to end on "/"
      if (path.substr(path.length() - 1) != "/")
        path = path + "/";
      if (spin != 0)
        scfMode = Options::SCF_MODES::UNRESTRICTED;
      return;
    }
    // unblocked options
    std::string name = word;
    if (!(iss >> word)) {
      throw SerenityError("ERROR: Value missing for keyword '" + name + "'.");
    }
    std::string value = word;
    this->set("", name, value);
  }
  /*
   * Cannot be reached
   */
  throw SerenityError("Error while parsing settings file.");
}

void Settings::set(std::string blockname, std::string name, std::string value) {
  /*
   * Definitions
   */
  bool check = false;

  // TODO this templated lambda function version
  //       of the visitor should work with C++14
  //
  //    auto visitor = [&check,&name,&value](auto f){
  //      if(name.compare(f.name()) == 0 ){
  //        std::cout << f.name() << " "<< typeid(f.get()).name() << std::endl;
  //        check = true;
  //      }
  //    };
  /*
   * Parsing
   * (Add new blocks here aswell)
   */

  set_visitor visitor(name, value, check);
  if (!blockname.compare("")) {
    visit_each(*this, visitor);
  }
  else if (!blockname.compare("DFT")) {
    visit_each((*this).dft, visitor);
  }
  else if (!blockname.compare("CUSTOMFUNC")) {
    visit_each(this->customFunc, visitor);
  }
  else if (!blockname.compare("SCF")) {
    visit_each((*this).scf, visitor);
  }
  else if (!blockname.compare("BASIS")) {
    visit_each((*this).basis, visitor);
  }
  else if (!blockname.compare("GRID")) {
    visit_each((*this).grid, visitor);
  }
  else if (!blockname.compare("EFIELD")) {
    visit_each((*this).efield, visitor);
  }
  else if (!blockname.compare("EXTCHARGES")) {
    visit_each((*this).extCharges, visitor);
  }
  else if (!blockname.compare("PCM")) {
    visit_each(this->pcm, visitor);
  }
  else {
    throw SerenityError("ERROR: No block '" + blockname + "' known.");
  }
  if (!check) {
    throw SerenityError("ERROR: No keyword '" + name + "' known in this block.");
  }
}
} /* namespace Serenity */
