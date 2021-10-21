/**
 * @file   MullikenPopulationCalculator.h
 *
 * @date   Mar 11, 2014
 * @author Thomas Dresselhaus
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
#ifndef MULLIKENPOPULATIONCALCULATOR_H_
#define MULLIKENPOPULATIONCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrix.h"
#include "settings/Options.h"

namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class ElectronicStructure;
template<Options::SCF_MODES SCFMode>
class MatrixInBasis;
template<Options::SCF_MODES SCFMode>
class SPMatrix;
template<class SCFMode>
class Matrix;

class SystemController;
/**
 * @class  MullikenPopulationCalculator MullikenPopulationCalculator.h
 * @brief  Performs a population analysis according to Mulliken.
 *
 * I.e. electrons are counted to an atom if they reside in basis functions which are centered
 * on that atom.
 */
template<Options::SCF_MODES SCFMode>
class MullikenPopulationCalculator {
 public:
  /**
   * @param   system of which the atoms, i.e. also the mapping of basis functions to atoms, are
   *                 taken
   * @param   electronicStructure of which the occupied orbitals, i.e. the density matrix are taken
   * @returns a vector holding the Mulliken population for each atom (in the order of the atoms
   *          as in system)
   */
  static SpinPolarizedData<SCFMode, Eigen::VectorXd>
  calculateMullikenPopulations(std::shared_ptr<SystemController> systemController);
  /**
   * @brief see calculateBasisFunctionPopulations()
   *
   * @param densityMatrix
   * @param overlapMatrix
   * @param atomBasisIndices see AtomCenteredBasisController, forms a connection between the basis
   *                         functions (i.e. their indices) and the atoms (i.e. the atom index)
   *                         For each atom there is a first index and an end index
   * @returns the Mulliken population on the atoms defined by the atomBasisIndices
   */
  static SpinPolarizedData<SCFMode, Eigen::VectorXd>
  calculateAtomPopulations(const DensityMatrix<SCFMode>& densityMatrix,
                           const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix,
                           const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices);
  /**
   * @brief calculates the population of each basis function
   *
   * @param densityMatrix P
   * @param overlapMatrix S
   * @returns the population of each basis function, i.e.
   *          \f$ {\rm Vector}_{\mu} = \sum_{\nu} P_{\mu,\nu} \cdot S_{\mu,\nu} \f$
   */
  static SpinPolarizedData<SCFMode, Eigen::VectorXd>
  calculateBasisFunctionPopulations(const DensityMatrix<SCFMode>& densityMatrix,
                                    const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix);
  /**
   * @brief Analog of calculateAtomPopulations() for orbitals (uses calculateOrbitalPopulations())
   * @param coefficients The orbital coefficients.
   * @param overlapMatrix The overlap matrix.
   * @param atomBasisIndices The atom to basis function map.
   */
  static SPMatrix<SCFMode>
  calculateAtomwiseOrbitalPopulations(const SPMatrix<SCFMode>& coefficients,
                                      const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix,
                                      const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices);
  /**
   * @brief Analog of calculateAtomPopulations() for orbitals.
   * @param system The system controller to calculate the populations for.
   */
  static SPMatrix<SCFMode> calculateAtomwiseOrbitalPopulations(std::shared_ptr<SystemController> system);

 private:
  /**
   * @brief analog to calculateBasisFunctionPopulations, but for a molecular orbital
   *
   * instead of a density matrix. I.e. How much an orbital would contribute to the Mulliken
   * population if it would be occupied.
   *
   * @param orbital
   * @param overlapMatrix
   * @returns a vector containing the information how much a basis function is 'populated' by an
   *          orbital, i.e.
   *          \f$ {\rm Vector}_{\mu} = \sum_{\nu} c_{\mu} \cdot c_{\nu} \cdot S_{\mu,\nu} \f$
   */
  static Eigen::VectorXd calculateOrbitalPopulations(const Eigen::VectorXd& orbitalcoeffitients,
                                                     const MatrixInBasis<Options::SCF_MODES::RESTRICTED>& overlapMatrix);
};

} /* namespace Serenity */

#endif /* MULLIKENPOPULATIONCALCULATOR_H_ */
