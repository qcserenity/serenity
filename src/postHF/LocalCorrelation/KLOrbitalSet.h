/**
 * @file KLOrbitalSet.h
 *
 * @author Moritz Bensberg
 * @date Dec 10, 2019
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

#ifndef POSTHF_LOCALCORRELATION_KLORBITALSET_H_
#define POSTHF_LOCALCORRELATION_KLORBITALSET_H_

/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h" //Overlap matrices.
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices
#include <memory>      //smrt_ptr

namespace Serenity {

/* Forward Declarations */
class CouplingOrbitalSet;
class OrbitalPair;

/**
 * @class KLOrbitalSet KLOrbitalSet.h
 * @brief A class that organizes a set of integrals etc. associated with the coupling
 *        between a pair ij and a pair kl.
 *        The class handles partially as a container and has some public variables that will
 *        be initialized during the run.
 */
class KLOrbitalSet {
 public:
  /**
   * @brief Constructor.
   * @param ijPair The ij pair.
   * @param klPair The kl pair.
   */
  KLOrbitalSet(std::shared_ptr<OrbitalPair> ijPair, std::shared_ptr<OrbitalPair> klPair);
  /**
   * @brief Getter for the overlap matrix [ij]x[kl]
   * @return
   */
  const Eigen::MatrixXd& getS_ij_kl();
  /**
   * @brief Getter for the overlap matrix [ij]x[k]
   * @return
   */
  const Eigen::MatrixXd& getS_ij_k();
  /**
   * @brief Getter for the overlap matrix [ij]x[l]
   * @return
   */
  const Eigen::MatrixXd& getS_ij_l();
  /**
   * @brief The integral (ik|jl).
   */
  double ik_jl; // 1x1 matrix
  /**
   * @brief The integral (lk|jk).
   */
  double il_jk; // 1x1 matrix
  /**
   * @brief The dressed integral tilde (ik|jl).
   *        See the pdf-document associated with DLPNO-CCSD for more information.
   */
  double tilde_ik_jl;
  /**
   * @brief The dressed integral tilde (il|jk).
   *        See the pdf-document associated with DLPNO-CCSD for more information.
   */
  double tilde_il_jk;
  /**
   * @brief The integrals (ki|la)
   */
  Eigen::VectorXd ki_la; // a in [j]
  /**
   * @brief The integrals (kj|la)
   */
  Eigen::VectorXd kj_la; // a in [i]
  /**
   * @brief The integrals (li|ka)
   */
  Eigen::VectorXd li_ka; // a in [j]
  /**
   * @brief The integrals (lj|ka)
   */
  Eigen::VectorXd lj_ka; // a in [i]
  /**
   * @brief Getter for kl pair.
   * @return The kl pair.
   */
  inline std::shared_ptr<OrbitalPair> getKLPair() {
    return _klPair.lock();
  }
  /**
   * @brief Initialize the overlap matrix controller.
   * @param domainSController
   */
  void setOverlapMatrixController(std::shared_ptr<DomainOverlapMatrixController> domainSController);
  /**
   * @brief Default destructor.
   */
  virtual ~KLOrbitalSet() = default;
  /**
   * @brief Remove integral vectors from memory.
   */
  void cleanUp();

 private:
  // The overlap matrix controller.
  std::weak_ptr<DomainOverlapMatrixController> _domainSController;
  // The ij pair.
  const std::weak_ptr<OrbitalPair> _ijPair;
  // The kl pair.
  const std::weak_ptr<OrbitalPair> _klPair;
  // The overlap matrix [ij]x[kl]
  std::shared_ptr<Eigen::MatrixXd> _s_ij_kl = nullptr;
  // The overlap matrix [ij]x[k]
  std::shared_ptr<Eigen::MatrixXd> _s_ij_k = nullptr;
  // The overlap matrix [ij]x[l]
  std::shared_ptr<Eigen::MatrixXd> _s_ij_l = nullptr;
};

} /* namespace Serenity */

#endif /* POSTHF_LOCALCORRELATION_KLORBITALSET_H_ */
