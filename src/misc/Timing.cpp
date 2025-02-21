/**
 * @file   Timing.cpp
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
/* Include Class Header*/
#include "misc/Timing.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
/* Include Std and External Headers */
#include <stdio.h>
#include <iostream>
#include <utility>

namespace Serenity {

void printTime(const long int time, const std::string label) {
  double timeTakenInSeconds = double(time) * 0.000000001;
  if (timeTakenInSeconds > 60.0) {
    int timeTakenInMinutes = (int)(timeTakenInSeconds / 60);
    print((std::string) "Time taken for " + label + ": " + timeTakenInMinutes + ":" +
          (timeTakenInSeconds - (60 * timeTakenInMinutes)) + " min.");
  }
  else {
    print((std::string) "Time taken for " + label + ": " + timeTakenInSeconds + " s.");
  }
  print("");
}

void takeTime(const std::string label) {
  if (STARTING_TIMES.find(label) == STARTING_TIMES.end()) {
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    STARTING_TIMES.emplace(label, ts);
  }
  else {
    std::cout << "Small error: You try to add a timer with a label which is already present." << std::endl;
    std::cout << "Your label is " << label << "." << std::endl;
  }
}

void timeTaken(const unsigned int printLevel, std::string label) {
  if (STARTING_TIMES.find(label) != STARTING_TIMES.end()) {
    timespec startingTime = STARTING_TIMES[label];
    timespec endTime;
    clock_gettime(CLOCK_REALTIME, &endTime);
    if (printLevel <= iOOptions.timingsPrintLevel) {
      printTime((endTime.tv_sec - startingTime.tv_sec) * 1000000000 + endTime.tv_nsec - startingTime.tv_nsec, label);
    }
    STARTING_TIMES.erase(label);
  }
  else {
    std::cout << "Small error: you want to print a timing of something of which there is no starting time stored."
              << std::endl;
    std::cout << "Maybe you already printed its timing (which stops the timer)? The label is " << label << "." << std::endl;
  }
}

void Timings::takeTime(const std::string label) {
  if (STARTING_TIMES.find(label) == STARTING_TIMES.end()) {
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    _tmpTimings.emplace(label, ts);
  }
  else {
    std::cout << "Small error: You try to add a timer with a label which is already present." << std::endl;
    std::cout << "Your label is " << label << "." << std::endl;
  }
}

void Timings::timeTaken(std::string label) {
  if (_tmpTimings.find(label) != _tmpTimings.end()) {
    timespec starttime = _tmpTimings[label];
    timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    time.tv_sec = time.tv_sec - starttime.tv_sec;
    time.tv_nsec = time.tv_nsec - starttime.tv_nsec;
    if (_timings.find(label) != _timings.end()) {
      _timings[label].tv_sec = _timings[label].tv_sec + time.tv_sec;
      _timings[label].tv_nsec = _timings[label].tv_nsec + time.tv_nsec;
    }
    else {
      _timings[label] = time;
    }
    _tmpTimings.erase(label);
  }
  else {
    std::cout << "Small error: you want to print a timing of something of which there is no starting time stored."
              << std::endl;
    std::cout << "Maybe you already printed its timing (which stops the timer)? The label is " << label << "." << std::endl;
  }
}

void Timings::printTimes() {
  for (auto it = _timings.begin(); it != _timings.end(); it++) {
    long int time = (it->second.tv_sec) * 1000000000 + it->second.tv_nsec;
    double timeTakenInSeconds = double(time) * 0.000000001;
    if (timeTakenInSeconds > 60.0) {
      int timeTakenInMinutes = (int)(timeTakenInSeconds / 60);
      printf("%4s %40s %6i:%06.3f min\n", "", it->first.c_str(), timeTakenInMinutes,
             (timeTakenInSeconds - (60 * timeTakenInMinutes)));
    }
    else {
      printf("%4s %40s %6i:%06.3f min\n", "", it->first.c_str(), 0, timeTakenInSeconds);
    }
  }
  clearTimes();
}

void Timings::clearTimes() {
  _tmpTimings.clear();
  _timings.clear();
}

std::map<std::string, timespec> Timings::_timings;
std::map<std::string, timespec> Timings::_tmpTimings;

} /* namespace Serenity */
