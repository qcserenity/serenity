/**
 * @file SingleSubstitution.h
 *
 * @date May 23, 2019
 * @author Moritz Bensberg
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

#ifndef DATA_SINGLESUBSTITUTION_H_
#define DATA_SINGLESUBSTITUTION_H_
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dense matrices
#include <Eigen/SparseCore> //Sparse matrices
#include <memory>           //smrt_ptr
#include <vector>           //std::vector

namespace Serenity {

/* Forward Declaration */
class OrbitalPair;

/**
 * @class SingleSubstitution SingleSubstitution.h
 * @brief A class that is concerned with the single amplitudes etc. in local-coupled cluster.
 */
class SingleSubstitution {
 public:
  /**
   * @brief Constructor.
   * @param diagonalPair The underlying diagonal pair.
   * @param singlesPNOScaling The scaling for the PNO threshold for this single.
   */
  SingleSubstitution(std::shared_ptr<OrbitalPair> diagonalPair, double singlesPNOScaling);
  /**
   * @brief Destructor.
   */
  ~SingleSubstitution() = default;
  /**
   * @brief The orbital index.
   */
  const unsigned int i;
  ///@brief The singles amplitudes.
  Eigen::VectorXd t_i;
  ///@brief The transformation to the linear independent quasi canonical PAO/PNO basis.
  Eigen::MatrixXd toPAODomain;
  ///@brief The residuals.
  Eigen::VectorXd residual;
  ///@brief The i, PNO fock matrix block.
  Eigen::VectorXd f_ai;
  ///@brief The dressed i,PNO fock matrix block.
  Eigen::VectorXd tf_ai;
  ///@brief The dressed a,b fock matrix block.
  Eigen::MatrixXd tf_ab;
  ///@brief The a,b fock matrix block.
  Eigen::VectorXd f_ab;
  ///@brief eps-F_ii
  Eigen::VectorXd epsMinusF;
  ///@brief The orbital pairs which contain i.
  ///       A loop over this vector is equivalent to a loop over all j.
  std::vector<std::weak_ptr<OrbitalPair>> orbitalPairs;
  ///@brief Flag for singles corresponding to core orbitals.
  bool coreLikeOrbital = false;
  /**
   * @brief Getter for the underlying diagonal pair.
   * @return The diagonal ii-pair.
   */
  inline std::shared_ptr<OrbitalPair> getDiagonalPair() {
    return diagPair.lock();
  }
  /**
   * @brief Getter for the PNO-Threshold.
   * @return The PNO threshold.
   */
  double getPNOThreshold();

 private:
  ///@brief The "dummy" diagonal pair which is used for the domains/PAOs.
  const std::weak_ptr<OrbitalPair> diagPair;
  ///@brief The PNO threshold.
  double _pnoThreshold;
};

} /* namespace Serenity */

#endif /* DATA_SINGLESUBSTITUTION_H_ */
