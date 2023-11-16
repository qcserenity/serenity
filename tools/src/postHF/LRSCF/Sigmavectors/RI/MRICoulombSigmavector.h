/**
 * @file MRICoulombSigmavector.h
 *
 * @date Nov 01, 2022
 * @author Niklas Niemeyer
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

#ifndef LRSCF_MRICOULOMBSIGMAVECTOR
#define LRSCF_MRICOULOMBSIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"

/* Include Std and External Headers */

namespace Serenity {

/**
 * @class MRICoulombSigmavector MRICoulombSigmavector.h
 * @brief Performs the calculation of the coulomb sigma vectors:
 *
 *          \f$ \tilde{F}_{ij} = \sum_{k l} \tilde{P}_{kl} (ij|kl) \f$
 * The point of this class is to calculate the Coulomb interaction between two subsystems using
 * not the union of their respective aux. bases (as done in the RICoulombSigmavector), but using
 * a series of metric multiplications so that
 *
 * \sigma_{(ia)_I} = (ia|P)(P|Q)^{-1}(Q|R)(R|S)^{-1}(S|jb) b_{jb},
 *
 * ia, P, Q belong to subsystem 1,
 * jb, R, S belong to subsystem 2.
 *
 * Proposed in J. Chem. Theory Comput. 2016, 12, 2, 549–557.
 *
 */
template<Options::SCF_MODES SCFMode>
class MRICoulombSigmavector : public Sigmavector<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   */
  MRICoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b);

  /**
   * @brief Default like constructor to be used if only the AO representation is required.
   * @param lrscf The (dummy) lrscf controller.
   */
  MRICoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf);

  virtual ~MRICoulombSigmavector() = default;

  // Don't use the pseudo Fock formulation.
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int, unsigned int, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>) final {
    return nullptr;
  };

  /**
   * @brief Override the conventional integral direct build.
   */
  void calcSigma() override;
}; // class MRICoulombSigmavector
} // namespace Serenity

#endif /* LRSCF_MRICOULOMBSIGMAVECTOR */
