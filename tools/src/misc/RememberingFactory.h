/**
 * @file   RememberingFactory.h
 *
 * @date   May 2, 2014
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
#ifndef REMEMBERINGFACTORY_H_
#define REMEMBERINGFACTORY_H_
/* Include Std and External Headers */
#include <cassert>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <tuple>

namespace Serenity {
/**
 * @class RememberingFactory RememberingFactory.h
 * @brief Abstract implementation of the concept of a caching factory.
 * @tparam Identifier a type determining all specifications needed to uniquely identify how an
 *                    object of type ProductT is to be produced.
 * @tparam ProductT the type of the class for which instances are to be produced by an
 *                  implementation of this factory.
 *
 * A RememberingFactory keeps track of all existing instances of the ProductT class.
 * How to use: Best take a look at other actual implementations. But here are the steps to go:
 * 1. Think about the Identifier template argument(s); it must contain all types to uniquely produce
 *    an object.
 * 2. Write a 'produce' method:
 *    1. Specify whatever arguments you think are needed (they should be system-dependent and
 *       configuration-independent).
 *    2. Call getOrProduce with an argument list which uses the Identifier types and return the result.
 * 3. Write the produceNew([Identifier... id]) method. Use ONLY the id arguments to produce the
 *    actual object.
 */
template<class ProductT, class... Identifier>
class RememberingFactory {
 public:
  RememberingFactory() = default;

  virtual ~RememberingFactory() = default;

 protected:
  /**
   * @param id The full specification of the object to be produced.
   * @returns a freshly produced instance
   *
   * TODO this should actually be static for each template specialization. Don't know whether
   * that's possible.
   * Must be overridden in each implementation / template specialization.
   */
  virtual std::unique_ptr<ProductT> produceNew(Identifier... id) = 0;
  /**
   * @param id The unique identification for the ProductT instance which is to be produced / looked
   *           up.
   * @returns an old instance if one with an Identifier equal to id has already been produced, and
   *          a new instance of the ProductT otherwise.
   *
   * This method should be called inside the 'produce' method of actual implementations. The
   * production of new objects itself must be pulled out into produceNew(Identifier id).
   */
  std::shared_ptr<ProductT> getOrProduce(Identifier... id) {
    std::tuple<Identifier...> idTuple(id...);
    std::lock_guard<std::mutex> lock(_lock);
    if (!(_producedInstances.find(idTuple) != _producedInstances.end())) {
      /*
       * Lets explain the madness happening in the line below:
       *
       * First of all we want to create a shared pointer.
       *
       * That shared pointer is generated from a the unique_ptr of the produceNew() function, but
       *  additionally a special deleter is added. The unique_ptr cannot be moved (because we
       *  override the deleter), but is instead released of the raw pointer.
       *
       * The default deleter of the shared pointer is then replaced by a custom one.
       *  The reason being: we want the shared_ptr to delete itself from the map of
       *  existing (remembered) classes. This leads to an automatic cleanUp of the map.
       *
       * The new deleter is a lambda function, that looks like the default one but calls
       *  the cleanUp() function after the usual delete.
       *
       * - JU
       *
       * Edit: TD
       */
      std::shared_ptr<ProductT> spt(produceNew(id...).release(), [](ProductT* p) {
        delete p;
        cleanUp();
      });
      _producedInstances[idTuple] = spt;
      return spt;
    }
    else if (_producedInstances[idTuple].expired()) {
      /*
       * Although the code below is perfectly valid, we should never end up here, because an
       * expired entry in the map should be deleted instantly due to the construction above.
       *
       * This assertion was actually observed to fail in parallel runs in the following situation:
       * An object is owned by one thread only and runs out of scope. Thus, the object will be
       * destroyed. At the same time another thread asks a RememberingFactory to produce the very
       * same object. It then may happen that the entry is still in the list, but the pointer
       * is already expired.
       *
       * If this occurs, you should change the code to avoid the above situation in the first
       * place. It can't be efficient to at the same time construct and destroy the very same
       * object.
       */
      assert(false);
      /*
       * See above
       */
      std::shared_ptr<ProductT> spt(produceNew(id...).release(), [](ProductT* p) {
        delete p;
        cleanUp();
      });
      _producedInstances[idTuple] = spt;
      return spt;
    }
    else {
      return std::shared_ptr<ProductT>(_producedInstances[idTuple]);
    }
  }

 private:
  /**
   * @brief Erase the entry from the _producedInstances map of which the product has just been destroyed.
   *
   * Upon destruction of an object earlier produced by the RememberingFactory we directly erase the
   * corresponding entry in the map of produced objects. For technical reasons we have to look
   * through the whole map, although only a single entry is erased.
   */
  static void cleanUp() {
    std::lock_guard<std::mutex> lock(_lock);
    auto itr = _producedInstances.begin();
#ifndef NDEBUG
    bool foundEntry = false;
#endif
    while (itr != _producedInstances.end()) {
      if (itr->second.expired()) {
        itr = _producedInstances.erase(itr);
#ifndef NDEBUG
        // Make sure there is only one single entry to be removed from the map
        assert(!foundEntry);
        foundEntry = true;
#else
        // If the check is not used (correct behaviour is assumed) we can just as well leave now.
        break;
#endif
      }
      else {
        ++itr;
      }
    }
#ifndef NDEBUG
    assert(foundEntry);
#endif
  }
  /**
   * A map holding all instances produced by the RememberingFactory. The key is the Identifier
   * template argument, because (as the name says) the equal products are identified by the
   * equality of the identifier.
   */
  static std::map<std::tuple<Identifier...>, std::weak_ptr<ProductT>> _producedInstances;
  /** Guards the _producedInstances map */
  static std::mutex _lock;
};

template<class ProductT, class... Identifier>
std::map<std::tuple<Identifier...>, std::weak_ptr<ProductT>> RememberingFactory<ProductT, Identifier...>::_producedInstances = {};

template<class ProductT, class... Identifier>
std::mutex RememberingFactory<ProductT, Identifier...>::_lock;

} /* namespace Serenity */

#endif /* REMEMBERINGFACTORY_H_ */
