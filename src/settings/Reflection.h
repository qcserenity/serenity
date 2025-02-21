/**
 * @file Reflection.h
 *
 * @date Jan 16, 2017
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

#ifndef REFLECTION_H_
#define REFLECTION_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"
/* Include Std and External Headers */
#include <boost/bind/bind.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/preprocessor.hpp>
#include <boost/type_traits/add_const.hpp>
#include <fstream>
#include <iostream>
#include <typeinfo>

/*
 * In this program, reflection is the ability of the code to examine its classes, interfaces, fields
 *   and methods at runtime without knowing the names of them at compile time. Reflection is not a native C++ ability,
 *   since the compiler is resolving this information into machine code. Therefore, the names of the variables are
 *   unknown by the class at runtime.
 *
 *   The following implementation is largely based on the
 *   following thread on stackoverflow.com:
 *   http://stackoverflow.com/questions/41453/how-can-i-add-reflection-to-a-c-application
 */
namespace Serenity {
namespace Reflection {

using namespace boost::placeholders;

#define REM(...) __VA_ARGS__
#define EAT(...)

// Retrieve the type
#define TYPEOF(x) DETAIL_TYPEOF(DETAIL_TYPEOF_PROBE x, )
#define DETAIL_TYPEOF(...) DETAIL_TYPEOF_HEAD(__VA_ARGS__)
#define DETAIL_TYPEOF_HEAD(x, ...) REM x
#define DETAIL_TYPEOF_PROBE(...) (__VA_ARGS__),
// Strip off the type
#define STRIP(x) EAT x
// Show the type without parenthesis
#define PAIR(x) REM x

// A helper metafunction for adding const to a type
template<class M, class T>
struct make_const {
  typedef T type;
};

template<class M, class T>
struct make_const<const M, T> {
  typedef typename boost::add_const<T>::type type;
};
#define REFLECTABLE(...)                                           \
  static const int fields_n = BOOST_PP_VARIADIC_SIZE(__VA_ARGS__); \
  friend struct reflector;                                         \
  template<int N, class Self>                                      \
  struct field_data {};                                            \
  BOOST_PP_SEQ_FOR_EACH_I(REFLECT_EACH, data, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define REFLECT_EACH(r, data, i, x)                                                   \
  PAIR(x);                                                                            \
  template<class Self>                                                                \
  struct field_data<i, Self> {                                                        \
    Self& self;                                                                       \
    field_data(Self& self) : self(self) {}                                            \
                                                                                      \
    typename make_const<Self, TYPEOF(x)>::type& get() { return self.STRIP(x); }       \
    typename boost::add_const<TYPEOF(x)>::type& get() const { return self.STRIP(x); } \
    const char* name() const { return BOOST_PP_STRINGIZE(STRIP(x)); }                 \
  };

/*
 * What this does is generate a constant fields_n that is number of reflectable fields in the class.
 *   Then it specializes the field_data for each field. It also friends the reflector class, this is
 *   so it can access the fields even when they are private.
 *   (http://stackoverflow.com/questions/41453/how-can-i-add-reflection-to-a-c-application)
 */

struct reflector {
  // Get field_data at index N
  template<int N, class T>
  static typename T::template field_data<N, T> get_field_data(T& x) {
    return typename T::template field_data<N, T>(x);
  }

  // Get the number of fields
  template<class T>
  struct fields {
    static const int n = T::fields_n;
  };
};

/*
 * Now to iterate over the fields we use the visitor pattern. We create an MPL range from 0 to the
 *   number of fields, and access the field data at that index. Then it passes the field data on to
 *   the user-provided visitor.
 *   (http://stackoverflow.com/questions/41453/how-can-i-add-reflection-to-a-c-application)
 */

struct field_visitor {
  template<class C, class Visitor, class I>
  void operator()(C& c, Visitor v, I) {
    v(reflector::get_field_data<I::value>(c));
  }
};

template<class C, class Visitor>
void visit_each(C& c, Visitor v) {
  typedef boost::mpl::range_c<int, 0, reflector::fields<C>::n> range;
  boost::mpl::for_each<range>(boost::bind<void>(field_visitor(), boost::ref(c), v, _1));
}

/**
 * @brief The set_visitor function checks if the variable and
 *          its value exist and sets the corresponding field.
 * @param n name of the field
 * @param v value of the field
 * @param c check if the field exists
 */
struct set_visitor {
  set_visitor(std::string n, std::string v, bool& c) : name(n), value(v), check(c){};
  std::string name;
  std::string value;
  bool& check;
  template<class FieldData>
  void operator()(FieldData f) {
    // f.name() gives the variable name of the field - e.g. gridType
    std::string fname = f.name();
    for (auto& c : fname)
      c = std::toupper(c);
    for (auto& c : name)
      c = std::toupper(c);
    if (name.compare(fname) == 0) {
      // f.get() gives a reference to the field - with a data type of e.g. Options::GRID_TYPES
      // here, value contains the correct data as a string, but we want to store it with the correct type in the field
      // (f.get())
      Options::resolve(value, f.get());
      check = true;
    }
  }
};

/**
 * @brief The print_visitor function prints the field and its value
 *          to a given stream
 * @param f the field
 * @param v value of the field
 * @param ofs the stream where the field and its value are printed to
 */
struct print_visitor {
  print_visitor(std::string& f, std::string& v, std::ofstream& ofs) : field(f), value(v), ofs(ofs){};
  std::string& field;
  std::string& value;
  std::ofstream& ofs;
  template<class FieldData>
  void operator()(FieldData f) {
    field = f.name();
    for (auto& c : field)
      c = std::toupper(c);
    try {
      // here, value is empty, but f.get() has the data, possibly as some enum class. resolve converts it to a string
      // and fills value
      Options::resolve(value, f.get());
      if (!value.empty()) {
        ofs << field << " " << value << std::endl;
      }
      value.clear();
    }
    catch (...) {
      value.clear();
    }
  }
};

} // namespace Reflection
} // namespace Serenity

#endif /* REFLECTION_H_ */
