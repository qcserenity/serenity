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
/* Include Std and External Headers */
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
   * @param lrscf All active LRSCFController of this response problem.
   *
   * References:
   *  - simplified TDA: J. Chem. Phys. 138, 244104 (2013)
   *  - simplified TDDFT: Computational and Theoretical Chemistry 1040-1041 (2014) 45-53
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
   * @return A reference to the intermediate used to calculate (ij|ab) exchange integrals (transformed with atom metric
   * matrix).
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>& getJij();

  /**
   * @brief Returns the intermediate used to calculate (ia|jb) Coulomb integrals.
   *
   * @return A reference to the intermediate used to calculate (ia|jb) Coulomb integrals (transformed with atom metric
   * matrix).
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>& getJia();

  /**
   * @brief Returns the "right" intermediate used to calculate (ij|ab) exchange integrals.
   *
   * @return A reference to the intermediate used to calculate (ij|ab) exchange integrals.
   */
  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>& getJab();

  /**
   * @brief Get the Gamma K object.
   *
   * @return Gamma K object for (ia|jb) integrals.
   */
  std::vector<std::vector<Eigen::MatrixXd>>& getGammaK();

  /**
   * @brief Get the Gamma J object.
   *
   * @return Gamma J object for (ab|ji) integrals.
   */
  std::vector<std::vector<Eigen::MatrixXd>>& getGammaJ();

  /**
   * @brief Returns the HF exchange ratio for this simplified TDDFT problem.
   * @return The HF exchange ratio for this simplified TDDFT problem.
   */
  std::vector<double> getHFExchangeRatio();

 private:
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jij;
  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jai;
  std::vector<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jab;

  std::vector<std::vector<Eigen::MatrixXd>> _gammaK;
  std::vector<std::vector<Eigen::MatrixXd>> _gammaJ;

  std::vector<double> _hfExchangeRatio;
}; /* class SimplifiedTDDFT */
} /* namespace Serenity */

#endif /* LRSCF_SIMPLIFIEDTDDFT */
