/**
 * @file FormattedOutputStream.h
 *
 * @date Mar 29, 2019
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

#ifndef IO_FORMATTEDOUTPUTSTREAM_H_
#define IO_FORMATTEDOUTPUTSTREAM_H_
/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

namespace Serenity {

///@brief The global print level. Every task manipulates this setting
///       in its parseGeneralSettings() routine.
extern Options::GLOBAL_PRINT_LEVELS GLOBAL_PRINT_LEVEL;

/**
 * @class IndentStreamBuf FormattedOutputStream.h
 * @brief A custom output stream buffer that only parses to the stream if
 *        the print level is higher or equal to the given template argument.
 */
template<Options::GLOBAL_PRINT_LEVELS PrintLevel>
class IndentStreamBuf : public std::stringbuf {
 public:
  /**
   * @brief Constructor.
   * @param str The output stream into which is printed.
   * @param indent The indentation set before each line if flush is called
   *               for the stream (always happens if std::endl is called).
   */
  IndentStreamBuf(std::ostream& str, std::string indent) : _output(str), _indent(indent) {
  }
  /**
   * @brief Destructor. Puts out charcters left in the buffer.
   */
  ~IndentStreamBuf() {
    if (pbase() != pptr()) {
      putOutput();
    }
  }
  /**
   * @brief Syncs stream and buffer.
   */
  virtual int sync() {
    putOutput();
    return 0;
  }
  /**
   * @brief Put out buffer if the print level is high enough.
   */
  void putOutput() {
    if (GLOBAL_PRINT_LEVEL >= PrintLevel) {
      _output << _indent << str();
    }
    str("");
    _output.flush();
  }

 private:
  std::ostream& _output;
  std::string _indent;
};

/**
 * @class FOut FormattedOutputStream.h
 * @brief Custom output stream that filters the given output according to
 *       a given print level.
 */
template<Options::GLOBAL_PRINT_LEVELS PrintLevel>
class FOut : public std::ostream {
  IndentStreamBuf<PrintLevel> buffer;

 public:
  /**
   * @brief Constructor.
   * @param str The unfiltered output stream.
   * @param indent Indentation of each line in the stream.
   */
  FOut(std::ostream& str, std::string indent) : std::ostream(&buffer), buffer(str, indent) {
  }
  /**
   * @brief Filtered << operator.
   * @param fout Left hand side stream.
   * @param out Right hand side.
   */
  template<typename U>
  friend FOut<PrintLevel>& operator<<(FOut<PrintLevel>& fout, const U& out) {
    if (GLOBAL_PRINT_LEVEL >= PrintLevel)
      static_cast<std::ostream&>(fout) << out;
    return fout;
  }
  /**
   * @brief Overwrites the put function of std::ostream. Filters according to
   *        the print level given.
   * @param c Character which should be put to the stream.
   */
  virtual FOut<PrintLevel>& put(char c) {
    if (GLOBAL_PRINT_LEVEL >= PrintLevel)
      static_cast<std::ostream&>(*this).put(c);
    return *this;
  }
};

/**
 * @class FPrintf FormattedOutputStream.h
 * @brief Custom printf function that filters the given output according to a given printlevel.
 */
template<Options::GLOBAL_PRINT_LEVELS PrintLevel>
class FPrintf {
 public:
  /**
   * @brief Constructor
   * @param indent Indentation of each line in the stream.
   */
  FPrintf(std::string indent) : _indent(indent) {
  }

  /**
   * @brief printf overwrites the standard printf function. Filters according to
   *        the print level given
   * @param format
   * @param ... variadic arguments, contains the arbitrary number of arguments that can be passed to the printf function
   */
  FPrintf<PrintLevel>& printf(const char* format, ...) {
    va_list args;
    va_start(args, format);
    if (GLOBAL_PRINT_LEVEL >= PrintLevel) {
      std::cout << _indent;
      vprintf(format, args);
    }
    va_end(args);
    return *this;
  }

 private:
  std::string _indent;
};

/**
 *  Static output stream for each print level. These can be used for filtered
 *  and indented output.
 */
namespace OutputControl {
static FOut<Options::GLOBAL_PRINT_LEVELS::MINIMUM> mOut(std::cout, "");
static FOut<Options::GLOBAL_PRINT_LEVELS::NORMAL> nOut(std::cout, "");
static FOut<Options::GLOBAL_PRINT_LEVELS::VERBOSE> vOut(std::cout, "==V==  ");
static FOut<Options::GLOBAL_PRINT_LEVELS::VERBOSE> vnOut(std::cout, "");
static FOut<Options::GLOBAL_PRINT_LEVELS::DEBUGGING> dOut(std::cout, "==D==  ");
/**
 * Static objects holding the printf() function for each print level.
 * These can be used for filtered and formatted output.
 * Usage: e.g. OutputControl::n.printf("Your output here!", yourArguments)
 */
static FPrintf<Options::GLOBAL_PRINT_LEVELS::MINIMUM> m("");
static FPrintf<Options::GLOBAL_PRINT_LEVELS::NORMAL> n("");
static FPrintf<Options::GLOBAL_PRINT_LEVELS::VERBOSE> v("==V==  ");
static FPrintf<Options::GLOBAL_PRINT_LEVELS::VERBOSE> vn("");
static FPrintf<Options::GLOBAL_PRINT_LEVELS::DEBUGGING> d("==D==  ");
} /* namespace OutputControl */

} /* namespace Serenity */

#endif /* IO_FORMATTEDOUTPUTSTREAM_H_ */
