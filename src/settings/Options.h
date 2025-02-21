/**
 * @file   Options.h
 *
 * @date   Jan 23, 2017
 * @author David Schnieders
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

#ifndef OPTIONS_H_
#define OPTIONS_H_

/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
#include "misc/WarningTracker.h"
/* Include Std and External Headers */
#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>
#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace Serenity {
namespace Options {

/**
 * @brief Resolve the input string to the given field.
 *        This effectively parses the input to the various options.
 * @param value The input string.
 * @param field The field to be parsed to.
 *
 * Note that this function has to be implemented for every possible option data type.
 * Also note that the resolve function is intended to work in both directions (from string to T when parsing the input,
 * from T to string e.g. when writing the settings to file).
 */
// Default resolve function.
template<class T>
void resolve(std::string& value, T& field) {
  (void)value;
  (void)field;
  throw SerenityError("ERROR: No mapping to string available for this type (Options::resolve)");
}
// Resolve functions for primitives, std::vectors<...> etc.
template<>
inline void resolve<int>(std::string& value, int& field) {
  if (value.empty()) {
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }
  else {
    try {
      field = std::stoi(value);
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into an integer.");
    }
  }
}
template<>
inline void resolve<unsigned int>(std::string& value, unsigned int& field) {
  if (value.empty()) {
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }
  else {
    try {
      if (std::stoi(value) < 0) { // negative number in input could lead to very high numbers
        WarningTracker::printWarning((std::string) "Converted '" + value + "' to its absolute value in the input.", true);
      }
      field = abs(std::stoi(value));
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into an unsigned integer.");
    }
  }
}
template<>
inline void resolve<unsigned long int>(std::string& value, unsigned long int& field) {
  if (value.empty()) {
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }
  else {
    try {
      if (std::stol(value) < 0) { // negative number in input could lead to very high numbers
        WarningTracker::printWarning((std::string) "Converted '" + value + "' to its absolute value in the input.", true);
      }
      field = abs(std::stol(value));
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into an unsigned long integer.");
    }
  }
}
template<>
inline void resolve<double>(std::string& value, double& field) {
  if (value.empty()) {
    std::ostringstream strs;
    strs << field;
    value = strs.str();
  }
  else {
    try {
      field = std::stod(value);
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a double.");
    }
  }
}
template<>
inline void resolve<bool>(std::string& value, bool& field) {
  if (value.empty()) {
    if (field) {
      value = "true";
    }
    else {
      value = "false";
    }
  }
  else {
    std::string copy = value;
    for (auto& c : copy)
      c = std::toupper(c);
    if (!copy.compare("TRUE") or !copy.compare("1")) {
      field = true;
    }
    else if (!copy.compare("FALSE") or !copy.compare("0")) {
      field = false;
    }
    else {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a boolean expression.");
    }
  }
}
template<>
inline void resolve<std::string>(std::string& value, std::string& field) {
  if (value.empty()) {
    value = field;
  }
  else {
    field = value;
  }
}
template<>
inline void resolve<std::vector<unsigned int>>(std::string& value, std::vector<unsigned int>& field) {
  if (value.empty()) {
    if (field.size()) {
      value = "{ ";
      for (unsigned int val : field) {
        std::string varAsString = boost::lexical_cast<std::string>(val);
        value += (varAsString + " ");
      }
      value += "}";
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?+-";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        field.push_back(std::stoi(word));
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of unsigned integers.");
    }
  }
}
template<>
inline void resolve<std::vector<std::vector<unsigned int>>>(std::string& value, std::vector<std::vector<unsigned int>>& field) {
  if (value.empty()) {
    for (unsigned outer = 0; outer < field.size(); ++outer) {
      for (unsigned inner = 0; inner < field[outer].size(); ++inner) {
        std::string varAsString = boost::lexical_cast<std::string>(field[outer][inner]);
        value += (varAsString + " ");
      }
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?+-";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    boost::replace_all(value, ";", " ; ");
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      field.push_back({});
      unsigned item = 0;
      while (iss >> word) {
        if (word == ";") {
          field.push_back({});
          item++;
        }
        else {
          field[item].push_back(std::stoi(word));
        }
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of vector of unsigned integers.");
    }
  }
}
template<>
inline void resolve<std::vector<std::vector<double>>>(std::string& value, std::vector<std::vector<double>>& field) {
  if (value.empty()) {
    for (unsigned outer = 0; outer < field.size(); ++outer) {
      for (unsigned inner = 0; inner < field[outer].size(); ++inner) {
        std::string varAsString = boost::lexical_cast<std::string>(field[outer][inner]);
        value += (varAsString + " ");
      }
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?+";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    boost::replace_all(value, ";", " ; ");
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      field.push_back({});
      unsigned item = 0;
      while (iss >> word) {
        if (word == ";") {
          field.push_back({});
          item++;
        }
        else {
          field[item].push_back(std::stod(word));
        }
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of vector of doubles.");
    }
  }
}
template<>
inline void resolve<std::vector<unsigned long int>>(std::string& value, std::vector<unsigned long int>& field) {
  if (value.empty()) {
    for (unsigned long int val : field) {
      std::string varAsString = boost::lexical_cast<std::string>(val);
      value += (varAsString + " ");
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?+-";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        field.push_back(std::stoul(word));
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of unsigned long integers.");
    }
  }
}
template<>
inline void resolve<std::vector<int>>(std::string& value, std::vector<int>& field) {
  if (value.empty()) {
    for (int val : field) {
      std::string varAsString = boost::lexical_cast<std::string>(val);
      value += (varAsString + " ");
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?+";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        field.push_back(std::stoi(word));
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of integers.");
    }
  }
}
template<>
inline void resolve<std::vector<bool>>(std::string& value, std::vector<bool>& field) {
  if (value.empty()) {
    for (bool val : field) {
      if (val) {
        value += "true ";
      }
      else {
        value += "false ";
      }
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?+-";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    field.clear();
    std::istringstream iss(value);
    std::string word;
    while (iss >> word) {
      std::string copy = word;
      for (auto& c : copy)
        c = std::toupper(c);
      if (!copy.compare("TRUE") or !copy.compare("1")) {
        field.push_back(true);
      }
      else if (!copy.compare("FALSE") or !copy.compare("0")) {
        field.push_back(false);
      }
      else {
        throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of booleans.");
      }
    }
  }
}
template<>
inline void resolve<std::vector<double>>(std::string& value, std::vector<double>& field) {
  if (value.empty()) {
    if (field.size()) {
      value = "{ ";
      for (double val : field) {
        std::string varAsString = boost::lexical_cast<std::string>(val);
        value += (varAsString + " ");
      }
      value += "}";
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        field.push_back(std::stod(word));
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of doubles.");
    }
  }
}
template<>
inline void resolve<std::vector<std::string>>(std::string& value, std::vector<std::string>& field) {
  if (value.empty()) {
    for (std::string val : field) {
      value += (val + " ");
    }
  }
  else {
    std::string bad_symbols = ",?\\'\"&*()^%$#@!{}[]|<>?+-";
    for (auto& c : bad_symbols) {
      if (value.find(c) != std::string::npos)
        throw SerenityError("ERROR: List inputs require spaces as delimiters.");
    }
    try {
      field.clear();
      std::istringstream iss(value);
      std::string word;
      while (iss >> word) {
        field.push_back(word);
      }
    }
    catch (...) {
      throw SerenityError("ERROR: Could not convert '" + value + "' into a vector of strings.");
    }
  }
}

/**
 * @brief A helper function converting a key string to a field of a map, or, conversely, find the first key that matches
 * the given field.
 * @param m The map containing pairs of std::string and members of class T.
 * @param key The key used to lookup map entries. If it is empty, fill it so that its map entry matches the given field.
 * @param field The class T object found as the key's entry. In reverse mode, this is the search value.
 */
template<class T>
void check(const std::map<std::string, T> m, std::string& key, T& field) {
  if (key.empty()) {
    for (auto const& pair : m) {
      if (pair.second == field) {
        key = pair.first;
        return;
      }
    }
    throw SerenityError("ERROR: Option map missing.");
  }
  else {
    try {
      for (auto& c : key)
        c = std::toupper(c);
      field = m.at(key);
    }
    catch (...) {
      std::string error = "ERROR: '" + key + "' is not a valid option.\n";
      error += "Valid options are:\n";
      for (const auto& myPair : m) {
        error += myPair.first + "\n";
      }
      throw SerenityError(error);
    }
  }
}
} /* namespace Options */

} /* namespace Serenity */

#endif /* OPTIONS_H_ */
