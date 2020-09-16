/**
 * @file NotifyingClass.h
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
#ifndef NOTIFYINGCLASS_H
#define NOTIFYINGCLASS_H
/* Include Serenity Internal Headers */
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>

namespace Serenity {
/**
 * @class NotifyingClass NotifyingClass.h
 * @brief Counterpart to ObjectSensitiveClass, template for classes which notify those.
 *
 * Also read the information about ObjectSensitiveClass and the corresponding entry in the wiki
 * (Concepts: Notification).\n
 * Whenever you need to inform other objects about changes in data your class is managing, it
 * should inherit from NotifyingClass<T>, where T is the type of the managed data object (e.g.
 * DensityOnGrid). After your an object of your class changed something in the managed data, it
 * has to call this->notifyObjects() for the notification system to work. Make sure that the
 * ObjectSensitiveClass<T> objects are actually signed in to this (i.e. they were given as an
 * argument in a call to this->addSensitiveObject()).\n
 *
 * Example:\n
 * We will work with some data class Foo. Now we will construct a Controller for it. This
 * controller will also inform other objects (if signed in) about changes that may happen.\n
 *
 * class FooController : public NotifyingClass<Foo> {
 * public:
 *   FooController() : _foo(nullptr), _upToDate(false) {}
 *   void makeInvalid() {
 *     // Something happens which invalidates _foo
 *     notifyObjects();
 *   }
 *   const Foo& getFoo() {
 *     if (!_upToDate) updateFoo();
 *     return *_foo;
 *   }
 * private:
 *   void updateFoo() {
 *     // Create or update the controlled foo object
 *     // ...
 *     _upToDate = true;
 *   }
 *   std::unique_ptr<Foo> _foo;
 *   void _upToDate;
 * };
 * \n
 * This construction will make sure that when calling getFoo() the result is always upToDate.
 * Any object that is signed in to the FooController (see ObjectSensitiveClass) will be
 * informed about a change. But only when a new/updated Foo is actually requested, the probably
 * expensive updateFoo() function is called.
 */
template<class T>
class NotifyingClass {
 public:
  NotifyingClass() = default;
  /**
   * @brief Constructor which already takes a list of objects to be informed about changes.
   */
  NotifyingClass(std::vector<std::weak_ptr<ObjectSensitiveClass<T>>> sensitiveObjects)
    : _sensitiveObjects(sensitiveObjects) {
  }
  virtual ~NotifyingClass() = default;
  /**
   * @brief Adds a sensitive object, which will be informed about changes
   * @param newSensitiveObject which is added to the _sensitiveObjects list
   *
   * Take a look at the class example
   */
  void addSensitiveObject(std::weak_ptr<ObjectSensitiveClass<T>> newSensitiveObject) const {
    _sensitiveObjects.push_back(newSensitiveObject);
  }

 protected:
  /**
   * @brief Informs all attached objects about changes in this object.
   *
   * This method must be called whenever a (critical) change happens inside your class.
   */
  void notifyObjects() const {
    for (const auto& object : _sensitiveObjects)
      if (!object.expired())
        object.lock()->notify();
  }

 private:
  /**
   * @brief The list of objects which will be informed about (critical) changes.
   *
   * This class should not limit the proper usage of the derived classes by adding/removing a
   * dependency information -> mutable
   */
  mutable std::vector<std::weak_ptr<ObjectSensitiveClass<T>>> _sensitiveObjects;
};

} /* namespace Serenity */
#endif /* NOTIFYINGCLASS_H */
