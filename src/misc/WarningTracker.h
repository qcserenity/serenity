/**
 * @file WarningTracker.h
 *
 * @date Mar 9, 2017
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

#ifndef MISC_WARNINGTRACKER_H_
#define MISC_WARNINGTRACKER_H_

/* Include Std and External Headers */
#include <memory>
#include <string>

namespace Serenity {

/**
 * @class WarningTracker WarningTracker.h
 * @brief Singleton to handle the warning which may occur during a run.
 *        Warnings are printed to a file and optionally to the output.
 */
class WarningTracker {
  /*
   * Private Constructor -> Singleton
   */
 private:
  WarningTracker() = default;

 public:
  /**
   * @brief Destructor.
   */
  virtual ~WarningTracker() = default;
  /**
   * @brief Getter for the instance of this singleton.
   * @return The instance.
   */
  static WarningTracker& getInstance();
  /**
   * @brief Increase the underlying task number.
   */
  static void increment() {
    _taskNmbr++;
  }
  /**
   * @brief Print the warning to file and to the output.
   * @param text The warning.
   * @param printtoOutput If true:  Print the warning to screen as well.
   *                      If false: Print the warning only to the file.
   */
  static void printWarning(std::string text, bool printtoOutput);

 private:
  // Current task number.
  static unsigned int _taskNmbr;
  // Current warning number.
  static unsigned int _nWarnings;
  // The instance of this singleton.
  static std::unique_ptr<WarningTracker> _instance;
};

} /* namespace Serenity */

#endif /* MISC_WARNINGTRACKER_H_ */
