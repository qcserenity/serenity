/**
 * @file   CoreCoreRepulsionDerivative.h
 *
 * @date   Nov 21, 2014
 * @author k_klah01
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
#ifndef MATHS_DERIVATIVES_CORECOREREPULSIONDERIVATIVE_H_
#define MATHS_DERIVATIVES_CORECOREREPULSIONDERIVATIVE_H_
/* Include Serenity Internal Headers */
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <array>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class Atom;
class Point;
/**
 * @class CoreCoreRepulsionDerivative CoreCoreRepulsionDerivative.h
 *
 * @brief calculates the partial derivatives of the Coulomb repulsion between two nuclei w.r.t.
 *        changes in the nuclear coordinates.
 */
class CoreCoreRepulsionDerivative {
  /* Purely static class. Must not be instantiated. */
  CoreCoreRepulsionDerivative() = delete;

 public:
  /**
   * @param   atoms
   * @returns First derivative of the total coulomb repulsion energy between the atoms w.r.t.
   *           changes in nuclear coordinates.
   */
  static Matrix<double> calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atoms);
  /**
   * @brief Calculate the derivative for the Coulomb interaction between these charges.
   * @param charges The charges.
   * @return First derivative of the total coulomb repulsion energy between the charges w.r.t.
   *           changes in charge coordinates.
   */
  static Matrix<double> calculateDerivative(const std::vector<std::pair<double, Point>>& charges);

  /**
   * @param   atoms active
   * @param   atoms environment
   * @returns First derivative of the total coulomb repulsion energy between the atoms in the environment and
   *           active region w.r.t. changes in nuclear coordinates.
   */
  static Matrix<double> calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atomsAct,
                                            const std::vector<std::shared_ptr<Atom>>& atomsEnv);

  /**
   * @param   chargesAct active
   * @param   chargesEnv environment
   * @returns First derivative of the total coulomb repulsion energy between the atoms in the environment and
   *           active region w.r.t. changes in nuclear coordinates.
   */
  static Matrix<double> calculateDerivative(const std::vector<std::pair<double, Point>>& chargesAct,
                                            const std::vector<std::pair<double, Point>>& chargesEnv);

  /**
   * @param   atomsAct active
   * @param   chargesEnv environment
   * @returns First derivative of the total coulomb repulsion energy between the atoms in the environment and
   *           active region w.r.t. changes in nuclear coordinates.
   */
  static Matrix<double> calculateDerivative(const std::vector<std::shared_ptr<Atom>>& atomsAct,
                                            const std::vector<std::pair<double, Point>>& chargesEnv);

  /**
   * @param   chargesAct active
   * @param   atomsEnv environment
   * @returns First derivative of the total coulomb repulsion energy between the atoms in the environment and
   *           active region w.r.t. changes in nuclear coordinates.
   */
  static Matrix<double> calculateDerivative(const std::vector<std::pair<double, Point>>& chargesAct,
                                            const std::vector<std::shared_ptr<Atom>>& atomsEnv);

 private:
  static std::vector<std::pair<double, Point>> convertAtomsToCharges(const std::vector<std::shared_ptr<Atom>>& atoms);
};

} /* namespace Serenity */
#endif /* MATHS_DERIVATIVES_CORECOREREPULSIONDERIVATIVE_H_ */
