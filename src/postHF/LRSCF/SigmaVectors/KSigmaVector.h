/**
 * @file KSigmaVector.h
 *
 * @date Dec 11, 2018
 * @author Michael Boeckers
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */

#ifndef LRSCF_KSIGMAVECTOR
#define LRSCF_KSIGMAVECTOR

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/SigmaVectors/SigmaVector.h"
#include "settings/Options.h"

/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/**
 * @class KSigmaVector KSigmaVector.h
 * @brief Performs the calculation of the exchange sigma vectors.
 */
template<Options::SCF_MODES SCFMode> class KSigmaVector : public SigmaVector<SCFMode> {
public:
    /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of guess vectors\n
   *          for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1], respectively.\n
   *          For TDA-like problems, the guesses for X are stored in b[0].\n
   *          Note that sigma vectors can, albeit not used in the present implementation, also be calculated for more than two sets of test vectors.
   * @param densityScreeningThreshold A prescreening threshold. Often, the matrix of guess vectors is rather sparse and has contributions\n
   *         from a few subsystems only, i.e. the density matrices of pure environment systems will be close to zero.\n
   *         If the maximum density matrix element of the density matrix of a specific subsystem is lower than this threshold,\n
   *         the calculation of the contribution of that block to the sigma vectors is skipped.
   * @param pm Signs for the exchange evaluation, the sign for the contribution to A and B needs to be changed for non-Hermitian eigenvalue problem:
   *        \f$ (A + B)_{ij} = \sum_{kl} [(il|jk) - (ik|jl)] \tilde{P}_{kl} \f$ 
   */
  KSigmaVector(
      std::vector<std::shared_ptr<LRSCFController<SCFMode> > > lrscf,
      std::vector<Eigen::MatrixXd> b,
      const double densityScreeningThreshold,
      const std::vector<int> pm);
  /**
   * @brief Default destructor.
   */
  virtual ~KSigmaVector() = default;

private:
  /**
  * @brief calculates the coulomb pseudo-fock Matrix contribution of the response matrix with respect to a set of guess vectors.
  * @param I/J the particular number of the subsystem.
  * @param P_J pseudo Density Matrix, see SigmaVector.h for definition.
  */
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode> > > > calcF(
      unsigned int I,
      unsigned int J,
      std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode> > > > P_J) final;
  ///@brief Signs for the exchange evaluation
  const std::vector<int> _pm;

};//class KSigmaVector
}//namespace Serenity

#endif /* LRSCF_KSIGMAVECTOR */
