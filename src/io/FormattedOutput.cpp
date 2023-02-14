/**
 * @file   FormattedOutput.cpp
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
/* Include Class Header*/
#include "io/FormattedOutput.h"
/* Include Serenity Internal Headers */
#include "memory/MemoryManager.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <omp.h>
#include <iomanip>
#include <iostream>
#include <vector>

namespace Serenity {

void setOutputOptions(const unsigned int digits) {
  std::cout.precision(digits);
  /**
   * TODO maybe add some cases here
   */
  std::cout << std::scientific;
}

void printProgramHead() {
  std::cout << std::endl;
  std::cout << "#==============================================================================#" << std::endl;
  std::cout << "#" << std::setw(78) << " "
            << "#" << std::endl;
  std::cout << "#" << std::setw(78) << "   SSSSS                                                    #########         #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "   SS    EEEE RRR  EEEE N   N II TTTTTT YY  YY          ####################  #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "   SS    E    R  R E    NN  N II   TT    YYYY        ######################## #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "   SSSSS EEEE RRR  EEEE N N N II   TT     YY      #########             ##### #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "      SS E    R R  E    N  NN II   TT     YY   #######       OOOOOOOO         #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "      SS EEEE R  R EEEE N   N II   TT     YY ########      OOOOOOOOOOOO       #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "   SSSSS                                  ########        OOOOOOOOOOOOOO      #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "                                     ###########   OOO    OOOOOOOOOOOOOO      #"
            << std::endl;
  std::cout << "#" << std::setw(78) << " #########                   ###############      OOOOO   OOOOOOOOOOOOOO      #"
            << std::endl;
  std::cout << "#" << std::setw(78) << " ######################################      OO   OOOOO    OOOOOOOOOOOO       #"
            << std::endl;
  std::cout << "#" << std::setw(78) << "       #########################         O   OO    OOO       OOOOOOOO         #"
            << std::endl;
  std::cout << "#" << std::setw(78) << " "
            << "#" << std::endl;
  std::cout << "#" << std::setw(78) << " "
            << "#" << std::endl;
  std::cout << "#" << std::setw(78) << std::left << centerHeadlines("Serenity") << "#" << std::endl;
  std::cout << "#" << std::setw(78) << " "
            << "#" << std::endl;
  std::cout << "#" << std::setw(78) << std::left << centerHeadlines("A quantum chemistry code") << "#" << std::endl;
  std::cout << "#" << std::setw(78) << std::left << centerHeadlines("developed in the group of Johannes Neugebauer")
            << "#" << std::endl;
  std::cout << "#" << std::setw(78) << std::left << centerHeadlines("at the WWU Münster.") << " #" << std::endl;
  std::cout << "#" << std::setw(78) << " "
            << "#" << std::endl;
  std::cout << "#==============================================================================#" << std::endl;
  std::cout << std::endl;
}

void printRunStartInfo() {
  // time the program run
  takeTime("the entire run");

  // get time
  time_t now = time(0);
  struct tm tstruct;
  char dateAndTime[50];
  tstruct = *localtime(&now);
  strftime(dateAndTime, sizeof(dateAndTime), "%Y-%m-%d %X", &tstruct);

  /**
   * FIXME works only for Linux
   */
  printSmallCaption("Program Info");
  std::cout << "    Version           :   1.5.1" << std::endl;

#ifndef GIT_BRANCH
#define GIT_BRANCH "UNKNOWN"
#endif
  std::cout << "    Git Branch        :   " << GIT_BRANCH << std::endl << std::endl;

  printSmallCaption("Program started");
  std::cout << "    Time              :   " << dateAndTime << std::endl;
  std::cout << "    By                :   " << getenv("USER") << std::endl;
  std::string hostName = (getenv("HOSTNAME") != NULL) ? getenv("HOSTNAME") : "HOSTNAME UNKNOWN";
  std::cout << "    On                :   " << hostName << std::endl << std::endl;

#ifdef _OPENMP
  std::cout << "    Threads           :   " << omp_get_max_threads() << std::endl;
#endif
  size_t mem = MemoryManager::getInstance()->getAvailableSystemMemory() / (1024 * 1024);
  std::cout << "    Memory Limit      :   " << mem << " MB" << std::endl;
}

void printRunEndInfo() {
  // get time
  time_t now = time(0);
  struct tm tstruct;
  char dateAndTime[50];
  tstruct = *localtime(&now);
  strftime(dateAndTime, sizeof(dateAndTime), "%Y-%m-%d %X", &tstruct);

  std::cout << std::endl;
  printSmallCaption("Final Timings");
  Timings::printTimes();
  std::cout << std::endl;
  printSmallCaption("Program ended");
  std::cout << "    "
            << "Time:  " << dateAndTime << std::endl;
  {
    std::string hostName;
    if (getenv("HOSTNAME") != NULL) {
      hostName = getenv("HOSTNAME");
    }
    else {
      hostName = "HOSTNAME UNKNOWN";
    }
    std::cout << "    "
              << "On:    " << hostName << std::endl;
  }
  timeTaken(0, "the entire run");
  std::cout << std::endl;
}

void printSectionTitle(const std::string text) {
  std::cout << std::endl;
  std::cout << "o------------------------------------------------------------------------------o" << std::endl;
  std::cout << "|" << centerHeadlines(text) << "|" << std::endl;
  std::cout << "o------------------------------------------------------------------------------o" << std::endl;
  std::cout << std::endl;
}

void printSubSectionTitle(const std::string text) {
  std::cout << std::endl;
  std::cout << center("------------------------------------------------------------") << std::endl;
  std::cout << center(text) << std::endl;
  std::cout << center("------------------------------------------------------------") << std::endl;
  std::cout << std::endl;
}

void printSmallCaption(const std::string text) {
  std::cout << "  " << text << ":" << std::endl;
  std::cout << " ";
  for (unsigned int i = 0; i < text.length() + 3; ++i) {
    std::cout << "-";
  }
  std::cout << std::endl;
}

void printTableHead(const std::string text) {
  std::cout << "  " << text << std::endl;
  std::cout << " ";
  for (unsigned int i = 0; i < text.length() + 3; ++i) {
    std::cout << "-";
  }
  std::cout << std::endl;
}

void printBigCaption(const std::string text) {
  std::cout << " o";
  for (unsigned int i = 0; i < text.length() + 3; ++i) {
    std::cout << "-";
  }
  std::cout << "o" << std::endl;
  std::cout << " | " << text << ": |" << std::endl;
  std::cout << " o";
  for (unsigned int i = 0; i < text.length() + 3; ++i) {
    std::cout << "-";
  }
  std::cout << "o" << std::endl;
  std::cout << std::endl;
}

void print(const std::string text) {
  std::cout << "    " << text << std::endl;
}

std::string centerHeadlines(const std::string s) {
  std::string outpt;
  int l = s.length();
  int pos = (int)((78 - l) / 2);
  for (int i = 0; i < pos; ++i)
    outpt = outpt + " ";
  outpt = outpt + s;
  for (int i = 0; i < pos + (78 - l) % 2; ++i)
    outpt = outpt + " ";
  return outpt;
}

std::string center(const std::string s) {
  std::string outpt;
  int l = s.length();
  int pos = (int)((80 - l) / 2);
  for (int i = 0; i < pos; ++i)
    outpt = outpt + " ";
  outpt = outpt + s;
  for (int i = 0; i < pos + (80 - l) % 2; ++i)
    outpt = outpt + " ";
  return outpt;
}

} /* namespace Serenity */
