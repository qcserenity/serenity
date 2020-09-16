/**
 * @file   FormattedOutput.h
 *
 * @date   Apr 14, 2014
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
#ifndef FORMATTEDOUTPUT_H_
#define FORMATTEDOUTPUT_H_
/* Include Serenity Internal Headers */
#include "io/IOOptions.h"
/* Include Std and External Headers */
#include <sstream>
#include <string>

namespace Serenity {
/**
 * General output precision switch for (floating point) numbers
 */
const unsigned int NUMBER_IN_STRINGSTREAM_PRECISION = 16;
/**
 * @brief Sets the output options for cout
 * TODO add other options, e.g. scientific notation vs. sth. else
 */
void setOutputOptions(unsigned int digits);

/**
 * @brief Prints the header of the program
 */
void printProgramHead();

/**
 * @brief Prints the info on start of the program
 */
void printRunStartInfo();

/**
 * @brief Prints the info on end of the program
 */
void printRunEndInfo();

/**
 * @brief Prints text into section box
 * @param text
 */
void printSectionTitle(std::string text);

/**
 * @brief Prints text into subsection box
 * @param text
 */
void printSubSectionTitle(std::string text);

/**
 * @brief Prints text into small caption, i.e. begin a sub-chapter of the output
 * @param text
 */
void printSmallCaption(std::string text);

/**
 * @brief Prints text into a larger caption, i.e. begin a new chapter of the output
 * @param text
 */
void printBigCaption(std::string text);

/**
 * @brief prints the text
 * @param text
 */
void print(std::string text);

/**
 *
 * @param text Same as printSmallCaption, just without colon.
 */
void printTableHead(std::string text);

/**
 * @brief centers text assuming line length of 78 chars
 * @param s
 */
std::string centerHeadlines(std::string s);

/**
 * @brief centers text assuming line length of 80 chars
 * @param s
 */
std::string center(std::string s);

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(std::string lhs, int rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << rhs;
  return lhs + convert.str();
}

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(int lhs, std::string rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << lhs;
  return convert.str() + rhs;
}

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(std::string lhs, unsigned int rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << rhs;
  return lhs + convert.str();
}

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(unsigned int lhs, std::string rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << lhs;
  return convert.str() + rhs;
}

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(std::string lhs, long unsigned int rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << rhs;
  return lhs + convert.str();
}

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(long unsigned int lhs, std::string rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << lhs;
  return convert.str() + rhs;
}

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(std::string lhs, double rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << rhs;
  return lhs + convert.str();
}

/**
 * @param lhs, rhs are added
 * @returns a string consisting of lhs and rhs
 */
inline std::string operator+(double lhs, std::string rhs) {
  std::ostringstream convert;
  convert.precision(NUMBER_IN_STRINGSTREAM_PRECISION);
  convert << lhs;
  return convert.str() + rhs;
}

} // namespace Serenity

#endif /* FORMATTEDOUTPUT_H_ */
