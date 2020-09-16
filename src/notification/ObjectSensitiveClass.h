/**
 * @file ObjectSensitiveClass.h
 *
 * @date Dec 29, 2014
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
#ifndef OBJECTSENSITIVECLASS_H
#define OBJECTSENSITIVECLASS_H
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward Declarations */
/**
 * @class ObjectSensitiveClass ObjectSensitiveClass.h
 * @brief Template for classes which need to react on changes in a certain external object.
 *
 * Also read the information about NotifyingClass and the corresponding entry in the wiki
 * (Concepts: Notification).\n
 * If a class is dependent on data of another class (T), it can simply inherit from
 * ObjectSensitiveClass<T> to clarify this dependency. There should be an object (notifyingObject)
 * of class NotifyingClass<T> (with the same T of course), which will notify an instance of
 * ObjectSensitiveClass<T> whenever the data (of type T) changes. For this to be possible, you have
 * to make sure that your instance of ObjectSensitiveClass<T> actually signs in to an instance of
 * NotifyingClass<T> by making a call to its method like:
 * notifyingObject.addSensitiveObject(this->_self);
 */
template<class T>
class ObjectSensitiveClass {
 public:
  /**
   * the _self member variable is initialized with an empty deleter, i.e. the shared_ptr _self
   * does not destroy the object when its internal counter runs to 0 (i.e. basically when the
   * shared_ptr itself is destroyed). It must be this way, because _self is only destroyed when
   * the object itself is destroyed and would cause an error if the shared_ptr would try and
   * destruct the object AGAIN.
   */
  ObjectSensitiveClass() : _self(this, [](ObjectSensitiveClass<T>*) {}) {
  }
  virtual ~ObjectSensitiveClass() = default;

  /**
   * @brief This method defines how the sensitive class reacts to a change in an object of type T.
   *
   * The content of this method shall be something like: Mark certain data as invalid or delete them.
   * New data shall, however, not be calculated through this call (if possible). This would
   * completely break lazy evaluation and could lead to a waste of calculation time. E.g. an object
   * to which an object of this class is sensitive may change several times before the object of
   * this class is actually used again.
   */
  virtual void notify() = 0;

 protected:
  /**
   * In case the ObjectSensitiveClass knows about the object that manages the data it is sensitive
   * to, it can directly tell the managing object of its existence.
   * I.e. in your constructor of the ObjectSensitiveClass call
   *    notifyingObject->addSensitiveObject(_self),
   * if you have a notifyingObject at hand.
   *
   * CAUTION: Please do NOT use this member variable for anything else.
   */
  const std::shared_ptr<ObjectSensitiveClass<T>> _self;
};

} /* namespace Serenity */
#endif /* OBJECTSENSITIVECLASS_H */
