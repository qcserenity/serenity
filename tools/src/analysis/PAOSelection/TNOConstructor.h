/**
 * @file TNOConstructor.h
 *
 * @author Moritz Bensberg
 * @date Oct 31, 2019
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

#ifndef ANALYSIS_PAOSELECTION_TNOCONSTRUCTOR_H_
#define ANALYSIS_PAOSELECTION_TNOCONSTRUCTOR_H_
/* Include Serenity Internal Headers */
#include "data/matrices/FockMatrix.h" //Fock matrix definition
/* Include Std and External Headers */
#include <memory> // smrt_ptr

namespace Serenity {
/* Forward Declarations */
class OrbitalTriple;
class PAOController;
class SystemController;
class OrbitalPair;
class SingleSubstitution;
/**
 * @class TNOConstructor TNOConstructor.h
 * @brief A constructor for the triple-natural orbitals (TNOs).\n
 *        See J. Chem. Phys. 139, 134101 (2013) for details.\n
 *        Possible environment orbitals are shifted by a level-shift in the
 *        virtual--virtual Fock matrix block.\n
 *        The TNOs are constructed from PAOs, by\n
 *           removing linear dependencies in the PAO set [ijk],\n
 *           projecting the pair density matrices in this basis,\n
 *           diagonalizing the triple-density matrix,\n
 *           truncating by occupation number,\n
 *           diagonalizing the Fock matrix in this truncated basis.
 */
class TNOConstructor {
 public:
  /**
   * @brief Constructor.
   * @param f The Fock matrix.
   * @param paoController The PAO controller.
   * @param paoOrthogonalizationThreshold The threshold for removing linear dependencies.
   * @param tnoThreshold The TNO truncation threshold.
   * @param coreScaling The scaling factor for triples containing core orbitals.
   * @param environmentSystems The environment systems.
   * @param levelShiftParameter The level-shift for the occupied environment orbitals.
   */
  TNOConstructor(std::shared_ptr<const FockMatrix<Options::SCF_MODES::RESTRICTED>> f,
                 std::shared_ptr<PAOController> paoController, double paoOrthogonalizationThreshold, double tnoThreshold,
                 double coreScaling, std::vector<std::shared_ptr<SystemController>> environmentSystems = {},
                 double levelShiftParameter = 1e+6);
  /**
   * @brief Transforms the virtual basis of the triple to the TNO basis
   *        and calculates all overlap matrices which will become necessary
   *        for DLPNO-CCSD(T0).
   * @param triple The triple.
   */
  void transformToTNOBasis(std::shared_ptr<OrbitalTriple> triple);
  /**
   * @brief Default destructor.
   */
  virtual ~TNOConstructor() = default;

 private:
  // Constructs a pair density matrix for the given pair.
  Eigen::MatrixXd constructDensityMatrix(std::shared_ptr<OrbitalPair> pair);
  // Constructs the overlap matrix between pair-PNOs and ijk-external basis.
  Eigen::MatrixXd constructOverlapMatrix(std::shared_ptr<OrbitalPair> pair, const Eigen::MatrixXd& s_p_ijk);
  // Constructs the overlap matrix between single-Pseudo NOs and ijk-external basis.
  Eigen::MatrixXd constructOverlapMatrix(std::shared_ptr<SingleSubstitution> single, const Eigen::MatrixXd& s_p_ijk);
  // The (shifted) Fock matrix.
  Eigen::MatrixXd _f_pao;
  // The PAO-Controller.
  std::shared_ptr<PAOController> _paoController;
  // The orthogonalization threshold for the linear dependencies.
  double _paoOrthogonalizationThreshold;
  // The TNO truncation threshold.
  double _tnoThreshold;
  // The TNO truncation threshold for triples containing core orbitals.
  double _tnoCoreThreshold;
};

} /* namespace Serenity */

#endif /* ANALYSIS_PAOSELECTION_TNOCONSTRUCTOR_H_ */
