/**
 * @file Looper_python.cpp
 *
 * @date Apr 26, 2018
 * @author: Jan Unsleber
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
#include "integrals/looper/CoulombInteractionIntLooper.h"
#include "integrals/looper/ExchangeInteractionIntLooper.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/looper/TwoElecThreeCenterIntLooper.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

template<class Looper>
void loop4c(Looper looper, const std::function<void(const unsigned int&, const unsigned int&, const unsigned int&,
                                                    const unsigned int&, const Eigen::VectorXd&, const unsigned int&)>& func) {
  // TODO
  // These lines are the way that parallel processing using the GIL
  //  should be working with pybin11 according to a few manuals an issue threads
  //  However, this does not seem to be the case.
  //  Obviously the generateion and usage of the integrals in parallel also in the
  //  Python interpreter is desirable. - JU
  //
  //  auto guarded_func = [&func](const unsigned int& i,
  //                          const unsigned int& j,
  //                          const unsigned int& k,
  //                          const unsigned int& l,
  //                          const Eigen::VectorXd& v,
  //                          const unsigned int& t) -> void {
  //    py::gil_scoped_acquire acquire;
  //    func(i,j,k,l,v,t);
  //  };
  //  looper.loop(guarded_func);

  auto n = omp_get_num_threads();
  omp_set_num_threads(1);
  looper.loop(func);
  omp_set_num_threads(n);
}
void loop3c(TwoElecThreeCenterIntLooper looper,
            const std::function<void(const unsigned int&, const unsigned int&, const unsigned int&,
                                     const Eigen::VectorXd&, const unsigned int&)>& func) {
  auto n = omp_get_num_threads();
  omp_set_num_threads(1);
  looper.loop(func);
  omp_set_num_threads(n);
}

void export_Looper(py::module& spy) {
  /* =====================
   *   4c-2e-Supersystem
   * =====================*/
  py::class_<TwoElecFourCenterIntLooper>(spy, "TwoElecFourCenterIntLooper",
                                         R"(A looper for all 2e-4c integrals of a given basis, uses symmetry.
   Due to problems with the GIL this class does not run in parallel,
   instead it automatically runs on a single core, no further user input needed.
   In the future this issue should, however, be resolved.)")
      .def(py::init<LIBINT_OPERATOR, const unsigned int, std::shared_ptr<BasisController>, double, double>(),
           R"(The constructor.
   Args:
     op      (INT_OPERATORS): The operator used in the integrals.
     deriv    (unsigned int): The number of geometrical derivatives to compute.
     basis           (Basis): The AO basis to loop over.
     prescreenTrhld (double): The threshold used to determine insignificant integrals using the Schwarz estimates.
     mu             (double): Parameter for range sepraration. To be used with erf_coulomb operator. [Optional]
   )",
           py::arg("op"), py::arg("deriv"), py::arg("basisController"), py::arg("prescreeningThreshold"), py::arg("mu") = 0)
      //  .def("loop",&loop4c,py::call_guard<py::gil_scoped_release>());
      .def("loop", &loop4c<TwoElecFourCenterIntLooper>,
           R"(The function performing the loop.
   This function expects a function as an argument.
   The function in the argument should return nothing (void)
   and capture output via optional parameters.
   An example running a sum over all computed integrals would look like this:
   > a = [0]
   > def eval(i,j,k,l,ints,tid,a=a):
   >   a[tid]+=ints[0]
   >   return 
   > looper.loop(eval)                                 
   The argument function needs to take the arguments:
   i,j,k,l            (unsigned int): The indeces of the 4c AO integrals.
   ints (Eigen:VectorXd/numpy.array): The integrals and derivatives. 
   tid                (unsigned int): The thread ID of the current thread.)");

  /* =============================
   *   4c-2e-Exchange Interaction
   * =============================*/
  py::class_<ExchangeInteractionIntLooper>(spy, "ExchangeInteractionIntLooper",
                                           R"(A looper for all 2e-4c exchange-like ordered integrals (i_A j_B|k_A l_B) of two given bases (A,B).
   Due to problems with the GIL this class does not run in parallel,
   instead it automatically runs on a single core, no further user input needed.
   In the future this issue should, however, be resolved.)")
      .def(py::init<LIBINT_OPERATOR, const unsigned int, std::shared_ptr<BasisController>, std::shared_ptr<BasisController>, double>(),
           R"(The constructor.
   Args:
     op      (INT_OPERATORS): The operator used in the integrals.
     deriv    (unsigned int): The number of geometrical derivatives to compute.
     basisOne        (Basis): The first AO basis to loop over.
     basisTwo        (Basis): The second AO basis to loop over.
     prescreenTrhld (double): The threshold used to determine insignificant integrals using the Schwarz estimates.
   )")
      //  .def("loop",&loop4c,py::call_guard<py::gil_scoped_release>());
      .def("loop", &loop4c<ExchangeInteractionIntLooper>,
           R"(The function performing the loop.
   This function expects a function as an argument.
   The function in the argument should return nothing (void)
   and capture output via optional parameters.
   An example running a sum over all computed integrals would look like this:
   > a = [0]
   > def eval(i,a,j,b,ints,tid,a=a):
   >   a[tid]+=ints[0]
   >   return 
   > looper.loop(eval)                                 
   The argument function needs to take the arguments:
   i,a,j,b            (unsigned int): The indeces of the 4c AO integrals, 
                                      i,j run over all AOs in the first basis,
                                      a,b run over all AOs in the second basis.
   ints (Eigen:VectorXd/numpy.array): The integrals and derivatives. 
   tid                (unsigned int): The thread ID of the current thread.)");

  /* =============================
   *   4c-2e-Coulomb Interaction
   * =============================*/
  py::class_<CoulombInteractionIntLooper>(spy, "CoulombInteractionIntLooper",
                                          R"(A looper for all 2e-4c Coulomb-like ordered integrals (i_A j_A|k_B l_B) of two given bases (A,B), 
   this looper uses the underlying symmetry.
   Due to problems with the GIL this class does not run in parallel,
   instead it automatically runs on a single core, no further user input needed.
   In the future this issue should, however, be resolved.)")
      .def(py::init<LIBINT_OPERATOR, const unsigned int, std::shared_ptr<BasisController>, std::shared_ptr<BasisController>, double>(),
           R"(The constructor.
   Args:
     op      (INT_OPERATORS): The operator used in the integrals.
     deriv    (unsigned int): The number of geometrical derivatives to compute.
     basisOne        (Basis): The first AO basis to loop over.
     basisTwo        (Basis): The second AO basis to loop over.
     prescreenTrhld (double): The threshold used to determine insignificant integrals using the Schwarz estimates.
   )")
      //  .def("loop",&loop4c,py::call_guard<py::gil_scoped_release>());
      .def("loop", &loop4c<CoulombInteractionIntLooper>,
           R"(The function performing the loop.
   This function expects a function as an argument.
   The function in the argument should return nothing (void)
   and capture output via optional parameters.
   An example running a sum over all computed integrals would look like this:
   > a = [0]
   > def eval(i,j,a,b,ints,tid,a=a):
   >   a[tid]+=ints[0]
   >   return 
   > looper.loop(eval)                                 
   The argument function needs to take the arguments:
   i,j,a,b            (unsigned int): The indeces of the 4c AO integrals, 
                                      i,j run over all AOs in the first basis,
                                      a,b run over all AOs in the second basis.
   ints (Eigen:VectorXd/numpy.array): The integrals and derivatives. 
   tid                (unsigned int): The thread ID of the current thread.)");

  /* =============================
   *   3c-2e-Supersystem
   * =============================*/
  py::class_<TwoElecThreeCenterIntLooper>(spy, "TwoElecThreeCenterIntLooper",
                                          R"(A looper for all 2e-4c Coulomb-like ordered integrals (i_A j_A|k_B) of two given bases (A,B), 
   this looper uses the underlying symmetry.
   Due to problems with the GIL this class does not run in parallel,
   instead it automatically runs on a single core, no further user input needed.
   In the future this issue should, however, be resolved.)")
      .def(py::init<LIBINT_OPERATOR, const unsigned int, std::shared_ptr<BasisController>,
                    std::shared_ptr<BasisController>, double, std::pair<unsigned int, unsigned int>>(),
           R"(The constructor.
   Args:
     op      (INT_OPERATORS): The operator used in the integrals.
     deriv    (unsigned int): The number of geometrical derivatives to compute.
     basisOne        (Basis): The first AO basis to loop over.
     basisTwo        (Basis): The second AO basis to loop over (aux. basis).
     prescreenTrhld (double): The threshold used to determine insignificant integrals using the Schwarz estimates.
     kRange    (pair/touple): The range for the shells in the second basis.
                              The first shell index is 0, the second index should end at the number of shells at maximum.
                              [Optional, default: all shells/ basis functions]
   )",
           py::arg("op"), py::arg("deriv"), py::arg("basisOne"), py::arg("basisTwo"), py::arg("prescreeningThreshold"),
           py::arg("kRange") = std::pair<unsigned int, unsigned int>(0, 0))
      //  .def("loop",&loop4c,py::call_guard<py::gil_scoped_release>());
      .def("loop", &loop3c,
           R"(The function performing the loop.
   This function expects a function as an argument.
   The function in the argument should return nothing (void)
   and capture output via optional parameters.
   An example running a sum over all computed integrals would look like this:
   > a = [0]
   > def eval(i,j,a,ints,tid,a=a):
   >   a[tid]+=ints[0]
   >   return 
   > looper.loop(eval)                                 
   The argument function needs to take the arguments:
   i,j,a              (unsigned int): The indeces of the 4c AO integrals, 
                                      i,j run over all AOs in the first basis,
                                      a runs over all AOs in the second basis.
   ints (Eigen:VectorXd/numpy.array): The integrals and derivatives. 
   tid                (unsigned int): The thread ID of the current thread.)");
}
