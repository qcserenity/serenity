/**
 * @file QST.h
 *
 * @date Apr 24, 2017
 * @author Marabel Riesmeier, cleaning and merge Jan Unsleber
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
#ifndef QST_H_
#define QST_H_

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/**
 * @class QST QST.h
 * @brief A Quadratic Synchronous Transit Algorithm.
 *
 * Localization of a transition state (i.e. saddle point on potential energy surface)
 * when starting point and product are known.
 *
 * References:\n
 * Halgren, T. A.; Lipscomb, W. N.
 * Chem. Phys. Lett. 1976, 49, 225–231,\n
 *               and\n
 * Bell, S.; Crighton, J.
 * J. Chem. Phys. 1984, 80, 2464–2475.
 */
class QST {
 public:
  /**
   * @brief Constructor.
   *
   * @param min1 The first minimum.
   * @param min2 The second minimum.
   * @param preopt If true, only runs a pre-optimization.
   *               Meaning [LST(max)->LST(o-min)->]QST(max)->QST(o-min)->QST(max).
   * @param guess Optional argument, takes a starting guess, skipping the initial LST cycle.
   */
  QST(const Eigen::VectorXd& min1, const Eigen::VectorXd& min2, bool preopt, std::unique_ptr<Eigen::VectorXd> guess = nullptr);
  /**
   * @brief Default destructor.
   */
  virtual ~QST() = default;

  /**
   * @brief The optimize function gets an update-function as variable.
   *        The update function updates the energies and forces at a given point.
   *        The update function has to be created by the user.
   * @param updateFunction lambda-function that updates energies and forces at a given position.
   *                       The update function needs a given point which is an Eigen::VectorXd
   *                       and updates the energy (double) and the force (Eigen::VectorXd) at this point.
   *                       The function may be given a boolean, which decides if gradient-information
   *                       is printed out or not.
   *
   * A short example:
   * @code
   * auto updateFunction=[](const Eigen::VectorXd& parameters,
   *                        double& value,
   *                        Eigen::VectorXd& gradients)->void{
   *   ...
   * };
   * Eigen::VectorXd axis = ...;
   * Eigen::vectorXd midpoint = ...;
   * double spread = 0.01;
   *
   * DimerMethod dimer(midPoint, spread, axis);
   * dimer.optimize(updateFunction);
   * @endcode
   */
  void optimize(std::function<void(const Eigen::VectorXd&, double&, Eigen::VectorXd&, bool)> updateFunction);

  /**
   * @brief Getter for the current tangent to the QST path.
   *
   * The tangent to the QST path should be similar to the largest imaginary frequency
   * at the given point.
   *
   * @return The tangent to the QST path.
   */
  Eigen::VectorXd getTangent() {
    return _tangent *= 1 / _tangent.norm();
  }

 private:
  // First minimum.
  Eigen::VectorXd _min1;
  // Second minimum.
  Eigen::VectorXd _min2;
  // Flag for pre-optimization.
  const bool _preopt;
  // Guess vector.
  std::unique_ptr<Eigen::VectorXd> _guess;
  // The transition state.
  Eigen::VectorXd _ts;
  // The gradients.
  Eigen::VectorXd _gradient;
  // The tangent along the largest imaginary frequency.
  Eigen::VectorXd _tangent;
  // The current energy.
  double _energy;
};

} /* namespace Serenity */
#endif /* QST_H_ */
