/**
 * @file   Task.h
 *
 * @date   Mar 7, 2014
 * @author Thomas Dresselhaus
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
#ifndef TASK_H_
#define TASK_H_
/* Include Serenity Internal Headers */
#include "settings/Reflection.h"
/* Include Std and External Headers */
#include <string>

namespace Serenity {
/**
 * @class Task Task.h
 * @brief Interface for jobs that can actually be performed.
 *
 * A Task can be set up in any possible way and is, if part of the task list, triggered during the
 * program run.
 */
class Task {
public:
  Task() = default;
  virtual ~Task() = default;
  /**
   * This method invokes the actual task. Everything may happen here.
   */
  virtual void run() = 0;

  /**
   * @brief Parses the settings to the task. Can be over defined by a new implementation
   *        in the respective task.
   * @param c The settings of any kind.
   * @param v The visitor.
   * @param blockname The name of the input block.
   */
  template<class C, class Visitor>
  void visit(C & c, Visitor v, std::string blockname) {
    (void) blockname;
    visit_each(c, v);
  }

};


} /* namespace Serenity */

#endif /* TASK_H_ */
