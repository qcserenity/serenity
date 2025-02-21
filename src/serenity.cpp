/**
 * @file serenity.cpp
 * @version 1.6.2
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

/* Include Serenity Internal Headers */
#include "input/Input.h"
#include "io/FormattedOutput.h"
#include "io/IOOptions.h"
#include "misc/SerenityError.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace Serenity;
/**
 * @brief The main function of the program.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
  // remove buffering on output
  setvbuf(stdout, NULL, _IONBF, BUFSIZ);
  setOutputOptions(12);
  printProgramHead();
  printRunStartInfo();

  if (argc < 2) {
    throw SerenityError("Too few arguments were given. Run Serenity with the input file as an argument.");
  }

  std::string inputFileName = argv[1];

  /*
   * Define task vector to work with
   */
  std::vector<std::unique_ptr<Task>> tasks;
  /*
   * Define a map for all the systems sorted by name
   */
  std::map<std::string, std::shared_ptr<SystemController>> systems;
  std::ifstream input(inputFileName);
  std::string line;
  while (getline(input, line)) {
    std::string word;
    std::istringstream iss(line);
    iss >> word;
    std::string upper = word;
    for (auto& c : upper)
      c = std::toupper(c);

    if (word.empty()) {
      continue;
    }
    else if (word[0] == '#') {
      continue;
    }
    else if (!upper.compare("+SYSTEM")) {
      Settings settings(input);
      systems[settings.name] = std::make_shared<SystemController>(settings);
    }
    else if (!upper.compare("+TASK")) {
      iss >> word;
      tasks.push_back(Input::parseTask(systems, word, input));
    }
    else {
      throw SerenityError("ERROR: Unknown text in input: '" + word + "'.");
    }
  }

  auto& warnings = WarningTracker::getInstance();

  /*
   * Loop through all the tasks that need to be done
   */
  printSectionTitle("Running Tasks");
  for (unsigned int tasknumber = 0; tasknumber < tasks.size(); ++tasknumber) {
    warnings.increment();
    printBigCaption((std::string) "Task No. " + (tasknumber + 1));
    takeTime((std::string) "task " + (tasknumber + 1));
    tasks[tasknumber]->parseGeneralSettings();
    tasks[tasknumber]->run();
    timeTaken(1, (std::string) "task " + (tasknumber + 1));
    // Task performed => can be deleted
    tasks[tasknumber].reset();
  }
  /*
   * Cleanup
   */
  printSectionTitle("Cleanup");
  tasks.clear();

  print("All clean and shiny");

  printRunEndInfo();
  return 0;
}
