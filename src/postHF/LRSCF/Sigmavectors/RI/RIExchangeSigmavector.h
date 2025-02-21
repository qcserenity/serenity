/**
 * @file RIExchangeSigmavector.h
 *
 * @date Oct 07, 2021
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

#ifndef LRSCF_RIEXCHANGESIGMAVECTOR
#define LRSCF_RIEXCHANGESIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Sigmavectors/ExchangeSigmavector.h"

/* Include Std and External Headers */

namespace Serenity {

/**
 * @class RIExchangeSigmavector RIExchangeSigmavector.h
 * @brief Performs the calculation of the exchange sigma vectors:
 *
 *          \f$ \tilde{F}_{ij} = \sum_{k l} \tilde{P}_{kl} (il|kj) \f$
 *
 * AR: The RIExchangeSigmavector inherits from the ExchangeSigmavector because both short- and long-range exchange are
 * handled by these classes and can be separately approximated using density fitting, which necessitates both classes to
 * know the same parameters such as the amount of HF exchange or the range-separation parameter. The
 * RIExchangeSigmavector mainly gets the setParameters() member function from the ExchangeSigmavector.
 */
template<Options::SCF_MODES SCFMode>
class RIExchangeSigmavector : public ExchangeSigmavector<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   * @param pm Signs for the exchange evaluation, the sign for the contribution to A and B needs to be changed for
   * non-Hermitian eigenvalue problem: \f$ (A + B)_{ij} = \sum_{kl} [(il|jk) - (ik|jl)] \tilde{P}_{kl} \f$
   * @param densFitK Use density fitting for HF exchange?
   * @param densFitLRK Use density fitting for LR exchange?
   */
  RIExchangeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b,
                        const std::vector<int> pm, bool densFitK, bool densFitLRK);

  /**
   * @brief Default like constructor to be used if only the AO representation is required.
   * @param lrscf The lrscf controller.
   * @param pm Signs for the exchange evaluation. In the case <pm> is set to 0, only the term
   *              \f$ \tilde{F}_{\mu\nu} = \sum_{\kappa\lambda} \tilde{P}_{\kappa\lambda} (\mu\kappa|\nu\lambda) \f$ is
   * calculated. <pm> = {1} results in \f$ \tilde{F}_{\mu\nu} = \sum_{\kappa\lambda} \tilde{P}_{\kappa\lambda}
   * \left((\mu\kappa|\nu\lambda) + (\mu\lambda|\nu\kappa) \right)\f$
   * @param densFitK Use density fitting for HF exchange?
   * @param densFitLRK Use density fitting for LR exchange?
   */
  RIExchangeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, const std::vector<int> pm,
                        bool densFitK, bool densFitLRK);

  /** @brief Calculates the exchange pseudo-Fock Matrix contribution of the response matrix with respect to a set of
   * guess vectors. The following calculation is performed: \f$ \tilde{F}_{ij} = \sum_{k l} \tilde{P}_{kl} (il|kj)
   * \f$\n\n
   *
   *          Note: This function was made public in order to have access to the AO-basis representation of the sigma
   * vector.
   *
   *   @param I/J the particular number of the subsystem
   *   @param P_J pseudo Density Matrix, see Sigmavector.h for definition
   */
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int I, unsigned int J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J) override final;

}; // class RIExchangeSigmavector
} // namespace Serenity

#endif /* LRSCF_RIEXCHANGESIGMAVECTOR */
