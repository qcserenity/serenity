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
 */
template<Options::SCF_MODES SCFMode>
class GrimmeSigmavector {
 public:
  /**
   * @brief Constructor.
   * @param lrscf The LRSCFController.
   */
  GrimmeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf);

  /**
   * @brief Destructor.
   */
  virtual ~GrimmeSigmavector();

  /**
   * @brief Used to calculate sigmavectors.
   *
   * @param guessVectors The guessvectors.
   * @param pm Signs in (A pm B) for the sigma vector calculations.
   *
   * @return The sigmavectors for guessvectors.
   */
  std::vector<Eigen::MatrixXd> getSigmavectors(std::vector<Eigen::MatrixXd>& guessVectors, std::vector<int> pm);

 private:
  ///@brief The LRSCFController.
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  ///@brief The SimplifiedTDDFT object belonging to this sigma vector.
  std::unique_ptr<SimplifiedTDDFT<SCFMode>> _simplifiedTDDFT;
}; // class GrimmeSigmavector
} // namespace Serenity

#endif /* LRSCF_GRIMMESIGMAVECTOR */
