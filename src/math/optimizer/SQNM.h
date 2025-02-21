/**
 * @file SQNM.h
 *
 * @date Oct 09, 2024
 * @author Kasibek Zumataev, reworked by Thorben Wiegmann
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
#ifndef MATHS_OPTIMIZER_SQNM_H_
#define MATHS_OPTIMIZER_SQNM_H_

/* Include Serenity Internal Headers */
#include "math/optimizer/Optimizer.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

/**
 * @class SQNM SQNM.h
 * @brief Implementation of the stabilized Quasi-Newton method (SQNM) according to
 * [1] J. Chem. Phys. 142, 034112 (2015)
 */
class SQNM : public Optimizer {
 public:
  /**
   * @brief Constructor.
   * @param parameters The parameters to be optimized.
   * @param historyLength The maximum number of previous steps to be considered in the algorithm.
   * @param epsilon The threshold to determine significant eigenvalues of the displacement overlap.
   * @param alpha The initial step length to be updated during the optimization.
   * @param energyThreshold The energy threshold to determine whether an optimization step is accepted.
   * @param trustRadius The maximum step length - step will be scaled to this length if the trust radius is exceeded.
   */
  SQNM(Eigen::VectorXd& parameters, unsigned int historyLength, double epsilon, double alpha, double energyThreshold,
       double trustRadius);
  /**
   * @brief Default destructor.
   */
  virtual ~SQNM() = default;
  /**
   * @brief See documentation of base class Optimizer.
   *
   * The update function parameter 'hessian' is ignored.
   */
  void optimize(std::function<bool(const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients,
                                   std::shared_ptr<Eigen::MatrixXd> hessian, bool print)>
                    updateFunction,
                std::shared_ptr<unsigned int> nRejected = nullptr) override final;
  /**
   * @brief Calculates differences of entries in a vector from a list of vectors and returns a list of vectors.
   * @param vec List of vectors from which the differences will be calculated.
   */
  std::vector<Eigen::VectorXd> getDifferences(std::vector<Eigen::VectorXd> vec);
  /**
   * @brief Normalizes a list of vectors and returns a matrix holding one normalized vector per column.
   * @param vec List of vectors to be normalized.
   */
  Eigen::MatrixXd calcNorms(std::vector<Eigen::VectorXd> vec);
  /**
   * @brief Checks if the trust radius is exceeded and scales the step accordingly if the trust radius is exceeded.
   * @param step Step to be checked and possibly scaled.
   */
  Eigen::VectorXd scaleStep(Eigen::VectorXd step);

 private:
  double _historyLength;
  double _epsilon;
  double _alpha;
  double _alphaStart;
  double _eThresh;
  double _trustRadius;
  // list of coordinates of previous steps - maximum length is _historyLength
  std::vector<Eigen::VectorXd> _coordList;
  // list of gradients of previous steps - maximum length is _historyLength
  std::vector<Eigen::VectorXd> _gradientList;
};

} // namespace Serenity

#endif /* MATHS_OPTIMIZER_SQNM_H_ */