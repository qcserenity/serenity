/**
 * @file FockSigmavector.h
 *
 * @date Dec 12, 2018
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

#ifndef LRSCF_FOCKSIGMAVECTOR
#define LRSCF_FOCKSIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Sigmavectors/Sigmavector.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/**
 * @class FockSigmavector FockSigmavector.h
 * @brief Performs the calculation of the delta energy sigma vectors
 */
template<Options::SCF_MODES SCFMode>
class FockSigmavector : public Sigmavector<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   */
  FockSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b);
  virtual ~FockSigmavector() = default;

 private:
  // If I != J, this function will return a nullptr since there will be no contributions.
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int I, unsigned int J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J) final;

  // Calculates the sigma vector.
  void addToSigma(unsigned int I, unsigned int iStartI,
                  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> F_IJ) override final;

}; // class FockSigmavector
} // namespace Serenity

#endif /* LRSCF_FOCKSIGMAVECTOR */
