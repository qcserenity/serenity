/**
 * @file GrimmeSigmavector.h
 *
 * @date Oct 01, 2021
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

#ifndef LRSCF_GRIMMESIGMAVECTOR
#define LRSCF_GRIMMESIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "settings/ElectronicStructureOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class LRSCFController;

template<Options::SCF_MODES SCFMode>
class SimplifiedTDDFT;

/**
 * @class GrimmeSigmavector GrimmeSigmavector.h
 * @brief Performs the calculation of the Coulomb/exchange sigma vectors based on Grimme's simplified TDA/TDDFT approach.
 *
 * For the TDA approach (A matrix only), see J. Chem. Phys. 138, 244104 (2013).
 * For the TDDFT approach (also the B matrix), see Comp. Theo. Chem. 1040–1041 (2014) 45–53.
 */
template<Options::SCF_MODES SCFMode>
class GrimmeSigmavector : public Sigmavector<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param lrscf The LRSCFController.
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   * @param simplifiedTDDFT The object holding all the integrals.
   */
  GrimmeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b,
                    std::vector<int> pm, std::shared_ptr<SimplifiedTDDFT<SCFMode>> simplifiedTDDFT);

  /**
   * @brief Destructor.
   */
  virtual ~GrimmeSigmavector();

  // Don't use the pseudo Fock formulation.
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int, unsigned int, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>) final {
    return nullptr;
  };

  /**
   * @brief Override the conventional integral direct build.
   */
  void calcSigma() override;

  ///@brief Signs for the exchange evaluation
  const std::vector<int> _pm;

  ///@brief The SimplifiedTDDFT object belonging to this sigma vector.
  std::shared_ptr<SimplifiedTDDFT<SCFMode>> _simplifiedTDDFT;

}; // class GrimmeSigmavector
} // namespace Serenity

#endif /* LRSCF_GRIMMESIGMAVECTOR */
