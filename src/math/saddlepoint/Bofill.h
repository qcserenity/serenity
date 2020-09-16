/**
 * @file Bofill.h
 *
 * @date Jul 13, 2017
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

#ifndef MATH_SADDLEPOINT_BOFILL_H_
#define MATH_SADDLEPOINT_BOFILL_H_

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/**
 * @class  Bofill Bofill.h
 *
 * Reference:
 * Phys. Chem. Chem. Phys., 2002, 4, 11–15
 * (Note that this is not the original Bofill reference but contains
 *  a good write up of the needed steps.)
 * The original paper by Bofill is the following:
 * J. Comput. Chem., 1994, 15, 1
 */
class Bofill {
 public:
  /**
   * @brief Constructor.
   * @param params The starting parameters.
   * @param trustRadius The trust radius for steps.
   * @param searchDirection
   */
  Bofill(const Eigen::VectorXd& params, double trustRadius, std::unique_ptr<Eigen::VectorXd> searchDirection = nullptr);
  /**
   * @brief Default destructor.
   */
  virtual ~Bofill() = default;

  /**
   * @brief The main optimization routine.
   * @param updateFunction The function used to get new gradients for a given set of parameters.
   */
  void optimize(std::function<void(const Eigen::VectorXd&, double&, Eigen::VectorXd&, bool)> updateFunction);

 private:
  Eigen::VectorXd _x;
  std::unique_ptr<Eigen::VectorXd> _ev;
  const double _trustradius;
};

} /* namespace Serenity */

#endif /* MATH_SADDLEPOINT_BOFILL_H_ */
