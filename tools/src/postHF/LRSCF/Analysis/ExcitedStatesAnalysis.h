/**
 * @file ExcitedStatesAnalysis.h
 *
 * @date Apr 21, 2021
 * @author Anton Rikus
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
#ifndef EXCITEDSTATESANALYSIS_H_
#define EXCITEDSTATESANALYSIS_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"

/* Include Std and External Headers */
#include <memory>

namespace Serenity {

/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class MatrixInBasis;
template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class  ExcitedStatesAnalysis ExcitedStatesAnalysis.h
 * @brief  Class to calculate transition densities, particle and hole densities after an
 * LRSCFTask was performed.
 */
template<Options::SCF_MODES SCFMode>
class ExcitedStatesAnalysis {
 public:
  /**
   * @brief Constructor.
   * @param excitations        A vector containing all transitions that shall be analyzed
   * @param lrscfController    An LRSCFController which can be used to gain information on the LRSCFTask
   */
  ExcitedStatesAnalysis(std::vector<unsigned int> excitations, std::shared_ptr<LRSCFController<SCFMode>> lrscfController);

  /**
   * @brief Default destructor.
   */
  virtual ~ExcitedStatesAnalysis() = default;

  /**
   * @brief Returns the transition density matrix in AO basis for the specified transition
   * @param n denotes the transition's index in the input vector
   */
  std::shared_ptr<MatrixInBasis<SCFMode>> getTransitionDensityMatrix(unsigned int n);

  /**
   * @brief Returns the hole density matrix in AO basis for the specified transition
   * @param n denotes the transition's index in the input vector
   */
  std::shared_ptr<MatrixInBasis<SCFMode>> getHoleDensityMatrix(unsigned int n);

  /**
   * @brief Returns the particle density matrix in AO basis for the specified transition
   * @param n denotes the transition's index in the input vector
   */
  std::shared_ptr<MatrixInBasis<SCFMode>> getParticleDensityMatrix(unsigned int n);

  /**
   * @brief Calculates the transition density matrix in AO basis for a specific transition and updates the member
   * variable transdensmatrix
   * @param n denotes the transition's index in the input vector
   */
  void calculateTransitionDensityMatrix(unsigned int n);

  /**
   * @brief Calculates the hole density matrix in AO basis for a specific transition and updates the member variable
   * holedensmatrix
   * @param n denotes the transition's index in the input vector
   */
  void calculateHoleDensityMatrix(unsigned int n);

  /**
   * @brief Calculates the particle density matrix in AO basis for a specific transition and updates the member variable
   * particledensmatrix
   * @param n denotes the transition's index in the input vector
   */
  void calculateParticleDensityMatrix(unsigned int n);

 private:
  /// @brief A vector which stores the transition density matrices for all transitions that are investigated.
  std::vector<std::shared_ptr<MatrixInBasis<SCFMode>>> _transdensmatrix;
  /// @brief A vector which stores the hole density matrices for all transitions that are investigated.
  std::vector<std::shared_ptr<MatrixInBasis<SCFMode>>> _holedensmatrix;
  /// @brief A vector which stores the particle density matrices for all transitions that are investigated.
  std::vector<std::shared_ptr<MatrixInBasis<SCFMode>>> _particledensmatrix;
  /// @brief A vector which stores the indices of the transitions that are to be analyzed by this object.
  std::vector<unsigned int> _excitations;
  /// @brief An LRSCFController which holds information on the LRSCFTask
  std::shared_ptr<LRSCFController<SCFMode>> _lrscfcontroller;
};

} /* namespace Serenity */
#endif /* EXCITEDSTATESANALYSIS_H_ */
