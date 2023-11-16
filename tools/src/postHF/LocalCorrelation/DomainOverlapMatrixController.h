/**
 * @file DomainOverlapMatrixController.h
 *
 * @author Moritz Bensberg
 * @date Dec 11, 2019
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

#ifndef POSTHF_LOCALCORRELATION_DOMAINOVERLAPMATRIXCONTROLLER_H_
#define POSTHF_LOCALCORRELATION_DOMAINOVERLAPMATRIXCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h" //MatrixInBasis
#include "math/Matrix.h"                 //Matrix for overlap block storing.
/* Include Std and External Headers */
#include <Eigen/Dense> //Dense matrices
#include <memory>      //smrt_ptr

namespace Serenity {

/* Forward Declarations */
class PAOController;
class OrbitalPair;
class SingleSubstitution;

/**
 * @class DomainOverlapMatrixController DomainOverlapMatrixController.h
 * @brief A class that calculates and holds the overlap matrices between virtual domains for local correlation
 *        approaches. This class will try to calculate all overlap matrices on construction! Lazy evaluation is
 *        not recommended since this would make a rather complicated lock-framework necessary in order to avoid
 *        calculating the same matrix with two different threads at the same time.
 */
class DomainOverlapMatrixController {
 public:
  /**
   * @brief Constructor.
   * @param paoController The PAO controller.
   * @param closeOrbitalPairs The orbital pairs.
   * @param singles The singles.
   * @param closeOrbitalPairIndices The indices of the orbital pairs.
   * @param nOcc The number of occupied orbitals in the system.
   */
  DomainOverlapMatrixController(std::shared_ptr<PAOController> paoController,
                                std::vector<std::shared_ptr<OrbitalPair>> closeOrbitalPairs,
                                std::vector<std::shared_ptr<SingleSubstitution>> singles,
                                const Eigen::MatrixXi closeOrbitalPairIndices, unsigned int nOcc);

  /**
   * @brief Getter for the overlap matrix [ij]x[kl].
   * @param ijPair The ij pair.
   * @param klPair The kl pair.
   * @return The overlap matrix.
   */
  const std::shared_ptr<Eigen::MatrixXd> getS(std::shared_ptr<OrbitalPair> ijPair, std::shared_ptr<OrbitalPair> klPair);
  /**
   * @brief Getter for the overlap matrix [ij]x[k].
   * @param ijPair The ij pair.
   * @param kSingle The k single.
   * @return The overlap matrix.
   */
  const std::shared_ptr<Eigen::MatrixXd> getS(std::shared_ptr<OrbitalPair> ijPair, std::shared_ptr<SingleSubstitution> kSingle);
  /**
   * @brief Getter for the overlap matrix [ij]x[k].
   * @param ijPair The ij pair.
   * @param kSingle The k single.
   * @return The overlap matrix.
   */
  const std::shared_ptr<Eigen::MatrixXd> getS(const OrbitalPair& ijPair, std::shared_ptr<SingleSubstitution> kSingle);
  /**
   * @brief Getter for the overlap matrix [k]x[l].
   * @param kSingle The k single.
   * @param lSingle The l single.
   * @return The overlap matrix.
   */
  const std::shared_ptr<Eigen::MatrixXd> getS(std::shared_ptr<SingleSubstitution> kSingle,
                                              std::shared_ptr<SingleSubstitution> lSingle);
  /**
   * @brief A debug function that allows to manually overwrite
   *        all returned overlap matrix with the given identity matrix.
   * @param identity The debug identity.
   */
  void setIdentity(std::shared_ptr<Eigen::MatrixXd> identity);
  /**
   * @brief Default destructor.
   */
  virtual ~DomainOverlapMatrixController() = default;
  /**
   * @brief Getter for the total size of the overlap matrix in bit.
   * @return The size of the overlap matrix.
   */
  double getOverlapMatrixSize();

 private:
  // The PAO controller.
  std::shared_ptr<PAOController> _paoController;
  // The orbital paris.
  std::vector<std::shared_ptr<OrbitalPair>> _closeOrbitalPairs;
  // The singles.
  std::vector<std::shared_ptr<SingleSubstitution>> _singles;
  // The orbital pair indices.
  const Eigen::MatrixXi _closeOrbitalPairIndices;
  // The singles indices.
  Eigen::VectorXi _singlesIndices;
  // Optional debug identity matrix.
  std::shared_ptr<Eigen::MatrixXd> _debugIdentity = nullptr;
  // The overlap matrices between pairs.
  std::unique_ptr<Matrix<std::shared_ptr<Eigen::MatrixXd>>> _s_ij_kl;
  // The overlap matrices between pairs and singles.
  std::unique_ptr<Matrix<std::shared_ptr<Eigen::MatrixXd>>> _s_ij_k;
  // The overlap matrices between singles.
  std::unique_ptr<Matrix<std::shared_ptr<Eigen::MatrixXd>>> _s_k_l;
};

} /* namespace Serenity */

#endif /* POSTHF_LOCALCORRELATION_DOMAINOVERLAPMATRIXCONTROLLER_H_ */
