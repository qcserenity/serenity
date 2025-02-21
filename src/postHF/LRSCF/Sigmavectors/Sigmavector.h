/**
 * @file Sigmavector.h
 *
 * @date Dec 06, 2018
 * @author Michael Boeckers, Johannes Toelle
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

#ifndef LRSCF_SIGMAVECTOR
#define LRSCF_SIGMAVECTOR

/* Include Serenity Internal Headers */
#include "data/matrices/MatrixInBasis.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class Sigmavector Sigmavector.h
 * @brief Interface for the calculation of the sigma vectors\n
 *        used in the orbital direct linear response implementation, that is:\n
 *        \f[ \sigma_{ia} = \sum_{jb} M_{ia,jb}b_{jb} \f]
 */
template<Options::SCF_MODES SCFMode>
class Sigmavector {
 public:
  /**
   * @brief Constructor
   * @param lrscf A controller for the lrscf calculation
   * @param b The sets of guess vectors. Generally, the response problem is non-Hermitian and one can have two sets of
   * guess vectors\n for right (X+Y) and left eigenvectors (X-Y). Per convention, these are stored in b[0] and b[1],
   * respectively.\n For TDA-like problems, the guesses for X are stored in b[0].\n Note that sigma vectors can, albeit
   * not used in the present implementation, also be calculated for more than two sets of test vectors.
   */
  Sigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, std::vector<Eigen::MatrixXd> b);

  /**
   * @brief Default like constructor to be used if only the AO representation is required.
   * @param lrscf The (dummy) lrscf controller.
   */
  Sigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf);

  ///@brief Default destructor.
  virtual ~Sigmavector() = default;

  /**
   * @brief Getter of the sigma vector.
   * @return Returns the sigmavector for each set of guess vectors.
   */
  const std::vector<Eigen::MatrixXd>& getSigma() {
    if (!_hasBeenCalculated) {
      this->calcSigma();
    }
    return _sigma;
  }

  std::vector<MatrixInBasis<SCFMode>> getPerturbedFockMatrix();

  /** @brief Function to calculate and return Fock-like matrix F_IJ.
   *  @param I/J the particular number of the subsystem
   *  @param P_J pseudo density matrices that are to be contracted with ERIs. The dimensions of P_J are nSets and within
   * that nGuess. Similarly, the return dimensions are nSets and nGuess (and for the MatrixInBasis obviously nBasis *
   * nBasis).
   *  @return Fock-like matrix F_IJ.
   */
  virtual std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
  calcF(unsigned int I, unsigned int J, std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> P_J) = 0;

 protected:
  /**
   * @brief Calculates and returns pseudo density matrices P_I for guess vectors of subsystem I
   *        \f$ P_{\lambda \sigma} = \sum_{jb} c_{\lambda j} b_{jb} c_{\sigma b} \f$
   */
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> calcP(unsigned int I);

  ///@brief Calculates the actual sigma vector.
  virtual void calcSigma();

  ///@brief Calculate and store sigma vectors for sets of guess vectors from subsystem I from Fock-like matrix.
  virtual void addToSigma(unsigned int I, unsigned int iStartI,
                          std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> F_IJ);

  ///@brief Returns an (nShells x nShells) matrix for the set of density matrices containing maximum entries for
  /// prescreening.
  void setShellWiseMaxDens(unsigned J, std::vector<std::vector<MatrixInBasis<SCFMode>>>& densityMatrices);

  ///@brief Decides if a specific subsystem interaction is to be skipped.
  bool skipThisInteraction(unsigned I, unsigned J);

  ///@brief Maximum entry of all density matrices.
  double _maxDens;

  ///@brief Shell-wise maximum density matrix entries.
  Eigen::MatrixXd _maxDensMat;

  ///@brief LRSCFController for each reference subsystem.
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  ///@brief Number of subsystems.
  const unsigned int _nSub;

  ///@brief The number of guess vectors for each set.
  const unsigned int _nGuess;

  ///@brief The number of guess vector sets.
  const unsigned int _nSet;

  ///@brief Number of threads used for Fock contractions.
  unsigned int _nThreads = omp_get_max_threads();

  ///@brief Number of threads used for Fock contractions before restriction.
  unsigned int _maxThreads = omp_get_max_threads();

  ///@brief Prescreening threshold.
  double _prescreeningThreshold = 1e-13;

  ///@brief Sets of guess vectors for each subsystem.
  std::vector<std::vector<Eigen::MatrixXd>> _b;

  ///@brief Sigma vectors for each set of guess vectors.
  std::vector<Eigen::MatrixXd> _sigma;

  ///@brief True if already calculated.
  bool _hasBeenCalculated;
}; // class Sigmavector
} // namespace Serenity

#endif /* LRSCF_SIGMAVECTOR */