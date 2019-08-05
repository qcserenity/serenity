/**
 * @file   MyClass.h
 *
 * @date   Dec 31, 2099
 * @author John Doe
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
#ifndef MYCLASS_H_
#define MYCLASS_H_

namespace Serenity {
using namespace std;

/**
 * @class MyClass MyClass.h
 *
 * @brief A short description of the class.
 *
 * A more extended and maybe very detailed description of this classes capabilities, why it exists,
 * how to use it, how NOT to use it and so on.
 */
class MyClass {
public:
  /**
   *  @param myConstructorArgument Here is explained what I am there for.
   */
  MyClass(MyType myConstructorArgument);
  virtual ~MyClass();

  /**
   * @brief This is a nice method.
   *
   * More information about this method.
   *
   * @param   myMethodArgument is to be used in this way
   * @returns some useful thing
   */
  MyReturnType myMethod(MyType myMethodArgument);

private:
   MyReturnType _myPrivateMethod();

   MyType _myPrivateMemberVariable;
};

} /* namespace QCpack */

#endif /* MYCLASS_H_ */
