/**
 * @file   Timing.h
 *
 * @date   Nov 21, 2013
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
#ifndef TIMING_H_
#define TIMING_H_
/* Include Std and External Headers */
#include <ctime>
#include <map>
#include <string>

namespace Serenity {

/**
 * @brief Starts a new timer
 * @param label A label which will be used in the output and also for identification of the timer.
 */
void takeTime(std::string label);
/**
 * @brief Stops a timer and prints the timing.
 * @param printLevel more important timings should get a smaller number here, so that the output
 *                   granularity can be switched.
 * @param label The label identifying the timer. A timer with exactly this label must have been
 *              started (takeTime()) with exactly this label.
 */
void timeTaken(unsigned int printLevel, std::string label);
/**
 * Contains all timers, i.e. their starting time connected with an identification string.
 */
static std::map<const std::string, timespec> STARTING_TIMES;
/**
 * Properly prints a time difference (between 'start' and 'end') without stopping a timer.
 *
 * @param time  A time difference in nanoseconds.
 * @param label the starting time of the timer with this label is used as the 'start' time.
 */
void printTime(long int time, std::string label);
/**
 * @class Timings Timing.h
 * @brief Class to store and accumulate timing over the course of a run.
 */
class Timings {
 public:
  /**
   * @brief Starts a new timer
   * @param label A label which will be used in the output and also for identification of the timer.
   */
  static void takeTime(std::string label);
  /**
   * @brief Stops a timer and prints the timing.
   * @param printLevel more important timings should get a smaller number here, so that the output
   *                   granularity can be switched.
   * @param label The label identifying the timer. A timer with exactly this label must have been
   *              started (takeTime()) with exactly this label.
   */
  static void timeTaken(std::string label);
  /**
   * @brief Prints all stored timings in a timings block.
   *
   * The timings are cleared afterwards.
   */
  static void printTimes();
  /**
   * @brief Clears all stored timings.
   */
  static void clearTimes();

 private:
  static std::map<std::string, timespec> _timings;
  static std::map<std::string, timespec> _tmpTimings;
};

} /* namespace Serenity */
#endif /* TIMING_H_ */
