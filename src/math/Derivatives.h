/**
 * @file Derivatives.h
 *
 * @date May 16, 2015
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
#ifndef DERIVATIVES_H
#define DERIVATIVES_H
/* Include Std and External Headers */
#include <memory>
#include <utility>

namespace Serenity {
/**
 * @class Gradient Derivatives.h
 * @brief A generic class for gradients, i.e. the first spatial derivative into each direction, of objects.
 *
 * Contains x, y and z components.\n
 * Hint: construct with makeGradient() or makeGradientPtr().\n
 * Hint: Loop over x, y and z with for(auto& component : myGradient) { ... }.\n
 * Can be copied and moved.
 */
template<class T>
struct Gradient {
  T x;
  T y;
  T z;
  class iterator {
   public:
    inline iterator(T* ptr) : ptr(ptr) {
    }
    inline iterator operator++() {
      iterator i(ptr);
      ++ptr;
      return i;
    }
    inline bool operator!=(const iterator& other) {
      return ptr != other.ptr;
    }
    //    inline bool operator==(const iterator& other) { return ptr == other.ptr; }
    inline T& operator*() const {
      return *ptr;
    }

   private:
    T* ptr;
  };
  class const_iterator {
   public:
    inline const_iterator(const T* ptr) : ptr(ptr) {
    }
    inline const_iterator operator++() {
      const_iterator i(ptr);
      ++ptr;
      return i;
    }
    inline bool operator!=(const const_iterator& other) {
      return ptr != other.ptr;
    }
    //    inline bool operator==(const const_iterator& other) { return ptr == other.ptr; }
    inline const T& operator*() const {
      return *ptr;
    }

   private:
    const T* ptr;
  };
  inline const const_iterator begin() const {
    return const_iterator(&x);
  }
  inline iterator begin() {
    return iterator(&x);
  }
  inline const const_iterator end() const {
    return const_iterator(&z + 1);
  }
  inline iterator end() {
    return iterator(&z + 1);
  }

  inline void operator+=(const Gradient<T>& rhs) {
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;
  }
  inline void operator-=(const Gradient<T>& rhs) {
    this->x -= rhs.x;
    this->y -= rhs.y;
    this->z -= rhs.z;
  }
  inline void operator*=(const double& rhs) {
    this->x *= rhs;
    this->y *= rhs;
    this->z *= rhs;
  }
};

/**
 * @brief Function to create a Gradient of Ts by forwarding the arguments to a constructor of T.
 *
 * I.e.: x, y and z component are all created in the same way.\n
 * Analogous to std::make_shared. Data is constructed on the heap.
 * TODO here and in the following: it might be dangerous to allow r-value arguments, because they
 * may be invalidated in the first call of T()
 */
template<class T, class... Args>
inline std::unique_ptr<Gradient<T>> makeGradientPtr(Args&&... args) {
  return std::unique_ptr<Gradient<T>>(new Gradient<T>{T(args...), T(args...), T(args...)});
  // Perfect forwarding does not work, obviously because the args are used more than once.
  //    std::forward<Args>(args)..., std::forward<Args>(args)..., std::forward<Args>(args)... });
}
/**
 * @brief Function to create a Gradient of Ts by forwarding the arguments to a constructor of T.
 *
 * I.e.: x, y and z component are all created in the same way.\n
 * Analogous to std::make_shared. Data is constructed on the stack.
 */
template<class T, class... Args>
inline Gradient<T> makeGradient(Args&&... args) {
  return Gradient<T>{T(args...), T(args...), T(args...)};

  // Perfect forwarding does not work, obviously because the args are used more than once.
  //    std::forward<Args>(args)..., std::forward<Args>(args)..., std::forward<Args>(args)... });
}

/**
 * @class Hessian Derivatives.h
 * @brief A generic class for the second spatial derivatives (mixed and non-mixed).
 *
 * Hint: construct with makeHessian() or makeHessianPtr().\n
 * Hint: Loop over xx, xy, xz, yy, yz and zz with for(auto& component : myHessian) { ... }.\n
 * Can be copied and moved.
 */
template<class T>
struct Hessian {
  T xx;
  T xy;
  T xz;
  T yy;
  T yz;
  T zz;
  class iterator {
   public:
    inline iterator(T* ptr) : ptr(ptr) {
    }
    inline iterator operator++() {
      iterator i(ptr);
      ++ptr;
      return i;
    }
    inline bool operator!=(const iterator& other) {
      return ptr != other.ptr;
    }
    //    inline bool operator==(const iterator& other) { return ptr == other.ptr; }
    inline T& operator*() const {
      return *ptr;
    }

   private:
    T* ptr;
  };
  class const_iterator {
   public:
    inline const_iterator(const T* ptr) : ptr(ptr) {
    }
    inline const_iterator operator++() {
      const_iterator i(ptr);
      ++ptr;
      return i;
    }
    inline bool operator!=(const const_iterator& other) {
      return ptr != other.ptr;
    }
    //    inline bool operator==(const const_iterator& other) { return ptr == other.ptr; }
    inline const T& operator*() const {
      return *ptr;
    }

   private:
    const T* ptr;
  };
  inline const const_iterator begin() const {
    return const_iterator(&xx);
  }
  inline iterator begin() {
    return iterator(&xx);
  }
  inline const const_iterator end() const {
    return const_iterator(&zz + 1);
  }
  inline iterator end() {
    return iterator(&zz + 1);
  }

  inline void operator+=(const Hessian<T>& rhs) {
    this->xx += rhs.xx;
    this->xy += rhs.xy;
    this->xz += rhs.xz;
    this->yy += rhs.yy;
    this->yz += rhs.yz;
    this->zz += rhs.zz;
  }
  inline void operator-=(const Hessian<T>& rhs) {
    this->xx -= rhs.xx;
    this->xy -= rhs.xy;
    this->xz -= rhs.xz;
    this->yy -= rhs.yy;
    this->yz -= rhs.yz;
    this->zz -= rhs.zz;
  }
};

/**
 * @brief Function to create a Hessian of Ts by forwarding the arguments to a constructor of T.
 *
 * I.e.: xx, xy, xz, yy, yz and zz component are all created in the same way.\n
 * Analogous to std::make_shared. Data is constructed on the heap.
 */
template<class T, class... Args>
inline std::unique_ptr<Hessian<T>> makeHessianPtr(Args&&... args) {
  return std::unique_ptr<Hessian<T>>(new Hessian<T>{T(args...), T(args...), T(args...), T(args...), T(args...), T(args...)});
  // Perfect forwarding does not work, obviously because the args are used more than once.
  //    std::forward<Args>(args)..., std::forward<Args>(args)..., std::forward<Args>(args)...,
  //    std::forward<Args>(args)..., std::forward<Args>(args)..., std::forward<Args>(args)...});
}
/**
 * @brief Function to create a Hessian of Ts by forwarding the arguments to a constructor of T.
 *
 * I.e.: xx, xy, xz, yy, yz and zz component are all created in the same way.\n
 * Analogous to std::make_shared. Data is constructed on the stack.
 */
template<class T, class... Args>
inline Hessian<T> makeHessian(Args&&... args) {
  return Hessian<T>{T(args...), T(args...), T(args...), T(args...), T(args...), T(args...)};
  // Perfect forwarding does not work, obviously because the args are used more than once.
  //    std::forward<Args>(args)..., std::forward<Args>(args)..., std::forward<Args>(args)...,
  //    std::forward<Args>(args)..., std::forward<Args>(args)..., std::forward<Args>(args)...});
}

} /* namespace Serenity */

#endif /* DERIVATIVES_H */
