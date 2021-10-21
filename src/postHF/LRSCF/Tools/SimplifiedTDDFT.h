/**
 * @file SimplifiedTDDFT.h
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

#ifndef LRSCF_SIMPLIFIEDTDDFT
#define LRSCF_SIMPLIFIEDTDDFT

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/ElectronicStructureOptions.h"

/* Include External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {

struct LRSCFTaskSettings;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

template<Options::SCF_MODES SCFMode>
class SimplifiedTDDFT {
 public:
  /**
   * @brief Constructor
   *
   * @param lrscf All active LRSCFController of this response problem.
   *
   * Note: currently only works for one controller. Coupled is planned for the future.
   *
   * References:
   *  - simplified TDA: J. Chem. Phys. 138, 244104 (2013)
   *  - simplified TDDFT: Computational and Theoretical Chemistry 1040–1041 (2014) 45–53
   */
  SimplifiedTDDFT(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf);

  /**
   * @brief Destructor.
   */
  virtual ~SimplifiedTDDFT() = default;

  /**
   * @brief Sets up intermediates used for the sigma vector calculation.
   */
  void setupSimplifiedTDDFT();

  /**
   * @brief Returns the "left" intermediate used to calculate (ij|ab) exchange integrals.
   *
   * This is not a constant reference since the SimplifiedTDDFTSigmavector needs to use it in a way
   * which prevents that. See also other intermediates.
   *
   * @return A reference to the intermediate used to calculate (ij|ab) exchange integrals.
   */
  SpinPolarizedData<SCFMode, Eigen::MatrixXd>& getJij();

  /**
   * @brief Returns the intermediate used to calculate (ia|jb) Coulomb integrals.
   *
   * @return A reference to the intermediate used to calculate (ia|jb) Coulomb integrals.
   */
  SpinPolarizedData<SCFMode, Eigen::MatrixXd>& getJai();

  /**
   * @brief Returns the "right" intermediate used to calculate (ij|ab) exchange integrals.
   *
   * @return A reference to the intermediate used to calculate (ij|ab) exchange integrals.
   */
  SpinPolarizedData<SCFMode, Eigen::MatrixXd>& getJab();

  /**
   * @brief Returns the HF exchange ratio for this simplified TDDFT problem.
   * @return The HF exchange ratio for this simplified TDDFT problem.
   */
  double getHFExchangeRatio();

 private:
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _Jij;
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _Jai;
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _Jab;

  double _hfExchangeRatio = 1.0;
}; /* class SimplifiedTDDFT */
} /* namespace Serenity */

#endif /* LRSCF_SIMPLIFIEDTDDFT */
