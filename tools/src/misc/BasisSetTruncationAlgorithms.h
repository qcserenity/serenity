/**
 * @file BasisSetTruncationAlgorithms.h
 *
 * @date Oct 30, 2017
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

#ifndef MISC_BASISSETTRUNCATIONALGORITHMS_H_
#define MISC_BASISSETTRUNCATIONALGORITHMS_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/* Forward declarations */
class SystemController;
namespace Options {
enum class BASIS_SET_TRUNCATION_ALGORITHMS;
}

/**
 * @class BasisSetTruncationAlgorithms BasisSetTruncationAlgorithms.h
 * @brief A class handling the different truncations algorithms used in projection based embedding.
 *
 * Only basis functions which are centered on dummy atoms of the active system can be truncated by these algorithms. All
 * basis functions which are centered on the non-dummy atoms of the active system are considered to be essential and are
 * not truncated.
 * Make sure that your system has dummy atoms. Otherwise calling methods of this class is pointless.
 *
 * Available truncation algorithms:\n
 * \f$ D^I_{ij}\f$ Density matrix of system I, entry \f$ij\f$.\n
 * \f$ S^I_{ij}\f$ Overlap matrix of system I entry \f$ij\f$.\n
 * \f$\chi_j^B \f$ A basis function centered in the environment.\n\n
 *
 * -NET_POPULATION:\n
 *   Keeps all basis functions with a Mulliken net population of more than the threshold <netPopThreshold>:\n
 *   According to: J. Chem. Phys. 143, 024105 (2015)\n
 *   \f$ q_j = D^A_{jj}S^A_{jj} \f$\n
 *   NOTE: Gives the best results of all truncation algorithms!\n\n
 *
 * -PRIMITIVE_NET_POPULATION:\n
 *   Keeps the <truncationFactor> most important basis functions of the environment. Truncation criterion for basis
 * function \f$ \chi_j^B \f$:\n \f$ q_j = D^A_{jj}S^A_{jj} \f$\n\n
 *
 */
class BasisSetTruncationAlgorithms {
 public:
  /**
   * @brief Default destructor.
   */
  virtual ~BasisSetTruncationAlgorithms() = default;

  /**
   * @brief Constructor.
   * @param activeSystem The system controller for the active system.
   */
  BasisSetTruncationAlgorithms(std::shared_ptr<SystemController> activeSystem);

  /**
   * @brief Selects the truncation scheme and calls the corresponding function.
   * @param truncationAlgorithm The employed truncation algorithm.
   * @param truncationFactor The factor used in primitive truncations.
   * @param dMatActiveSystem The total density matrix of the active system.
   * @param netPopThreshold The threshold used in the "net population" truncations.
   */
  void truncateBasis(Options::BASIS_SET_TRUNCATION_ALGORITHMS truncationAlgorithm, double truncationFactor = 0.0,
                     std::shared_ptr<DensityMatrix<Options::SCF_MODES::RESTRICTED>> dMatActiveSystem = nullptr,
                     double netPopThreshold = 0.0);

 private:
  /// @brief The active system
  std::shared_ptr<SystemController> _activeSystem;

  /* Helper functions */
  /** @brief Builds the truncated basis from a list of important shells. Active AOs always included.
   *
   * @param importantShells A vector with 1 if the shell is important and 0 if the shell can be disregarded
   */
  inline void buildTruncatedBasis(Eigen::VectorXi importantShells);
  /**
   * @brief "One Parameter" basis set truncation scheme proposed by Bennie et al.
   *        According to: J. Chem. Phys. 143, 024105 (2015)
   * @param dMatActiveSystem The total density matrix of the active system.
   * @param netPopThreshold The truncation threshold.
   */
  Eigen::VectorXi netPopulationTruncation(const DensityMatrix<Options::SCF_MODES::RESTRICTED> dMatActiveSystem,
                                          double netPopThreshold);

  /**
   * @brief Similar to the "primitive" truncation, employing the criterion of the netPopulationTruncation.
   * @param dMatActiveSystem The total density matrix of the active system.
   * @param truncationFactor The ratio of kept environmental basis functions.
   */
  Eigen::VectorXi primitiveNetPopTruncation(const DensityMatrix<Options::SCF_MODES::RESTRICTED> dMatActiveSystem,
                                            double truncationFactor);
};

} /* namespace Serenity */

#endif /* MISC_BASISSETTRUNCATIONALGORITHMS_H_ */
