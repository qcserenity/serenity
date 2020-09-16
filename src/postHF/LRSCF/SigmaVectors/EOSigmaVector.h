/**
 * @file EOSigmaVector.h
 *
 * @date April 14, 2018
 * @author Michael Boeckers, Johannes Tölle
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

#ifndef LRSCF_EOSIGMAVECTOR
#define LRSCF_EOSIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/SigmaVectors/SigmaVector.h"

namespace Serenity {
/**
 * @class EOSigmaVector EOSigmaVector.h
 * @brief Calculates the External Orthogonality (EO) coupling matrix contribution for exact sTDDFT/sTDA calculations.\n
 *        For more information the reader is referred to J. Toelle, M. Boeckers, J. Neugebauer, J. Chem. Phys 2019.\n
 *        Further publication for Huzinaga and Hoffmann coming.
 */
template<Options::SCF_MODES SCFMode>
class EOSigmaVector : public SigmaVector<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   * @param densityScreeningThreshold A prescreening threshold. Often, the matrix of guess vectors is rather sparse and
   * has contributions\n from a few subsystems only, i.e. the density matrices of pure environment systems will be close
   * to zero.\n If the maximum density matrix element of the density matrix of a specific subsystem is lower than this
   * threshold,\n the calculation of the contribution of that block to the sigma vectors is skipped.
   * @param levelShiftParameter The levelshift for the levelshift projection operator.
   * @param eoPot The embedding mode used for the sTDDFT calculation.
   */
  EOSigmaVector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b,
                const double densityScreeningThreshold, const double levelShiftParameter,
                const Options::KIN_EMBEDDING_MODES eoPot);
  virtual ~EOSigmaVector() = default;

 private:
  // A function to calculate and return Fock-like matrix F_IJ
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int I, unsigned int J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J) final;

  // The levelshift for the levelshift projection operator.
  double _levelShiftParameter;

  // The embedding mode used in the sTDDFT calculation.
  Options::KIN_EMBEDDING_MODES _eoPot;
};

} /* namespace Serenity */

#endif /* LRSCF_EOSIGMAVECTOR */
