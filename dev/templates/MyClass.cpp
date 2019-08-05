/**
 * @file   MyClass.cpp
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
#include "MyClass.h"

namespace Serenity {
using namespace std;

MyClass::MyClass(MyType myConstructorArgument) :
    _myPrivateMemberVariable(myConstructorArgument) {
  // Only if MyType is a GarbageCollectedClass
  _myPrivateMemberVariable->signIn();
}

MyClass::~MyClass() {
  // Only if MyType is a GarbageCollectedClass
  _myPrivateMemberVariable->signOut();
}

MyReturnType MyClass::myMethod(MyType myMethodArgument) {
  // ...
}

MyReturnType MyClass::_myPrivateMethod() {
  // do private things
}

} /* namespace QCpack */

