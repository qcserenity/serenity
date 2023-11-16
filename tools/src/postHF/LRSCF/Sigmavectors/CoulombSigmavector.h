/**
 * @file CoulombSigmavector.h
 *
 * @date Dec 07, 2018
 * @author Niklas Niemeyer, Johannes Toelle
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

#ifndef LRSCF_COULOMBSIGMAVECTOR
#define LRSCF_COULOMBSIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
/* Include Std and External Headers */

namespace Serenity {

/**
 * @class CoulombSigmavector CoulombSigmavector.h
 * @brief Performs the calculation of the coulomb sigma vectors:
 *
 *          \f$ \tilde{F}_{ij} = \sum_{k l} \tilde{P}_{kl} (ij|kl) \f$
 */
template<Options::SCF_MODES SCFMode>
class CoulombSigmavector : public Sigmavector<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   */
  CoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b);

  /**
   * @brief Default like constructor to be used if only the AO representation is required.
   * @param lrscf The (dummy) lrscf controller.
   */
  CoulombSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf);

  virtual ~CoulombSigmavector() = default;

  /** @brief Calculates the Coulomb pseudo-Fock Matrix contribution of the response matrix with respect to a set of guess
   * vectors. The following calculation is performed: \f$ \tilde{F}_{ij} = \sum_{k l} \tilde{P}_{kl} (ij|kl) \f$\n\n
   *
   *          Note: This function was made public in order to have access to the AO-basis representation of the sigma
   * vector.
   *
   *   @param I/J the particular number of the subsystem
   *   @param P_J pseudo Density Matrix, see Sigmavector.h for definition
   */
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int I, unsigned int J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J) final;

}; // class CoulombSigmavector
} // namespace Serenity

#endif /* LRSCF_COULOMBSIGMAVECTOR */
