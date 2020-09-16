/**
 * @file EOSigmaVector.cpp
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

#include "EOSigmaVector.h"

/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "grid/GridControllerFactory.h"
#include "integrals/wrappers/Libint.h"

/* ABFockMatrixConstruction */
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundle.h"
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundleFactory.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
EOSigmaVector<SCFMode>::EOSigmaVector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                      std::vector<Eigen::MatrixXd> b, const double densityScreeningThreshold,
                                      const double levelShiftParameter, const Options::KIN_EMBEDDING_MODES eoPot)
  : SigmaVector<SCFMode>(lrscf, b, densityScreeningThreshold), _levelShiftParameter(levelShiftParameter), _eoPot(eoPot) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
EOSigmaVector<SCFMode>::calcF(unsigned int I, unsigned int J,
                              std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> densityMatrices) {
  Timings::takeTime("LRSCF -  Fock-like matrix: EO");

  assert(I != J && "No EO contribution for intra-subsystem blocks.");

  // Set dimensions for Fock like matrices
  // Final dimensions are the dimensions of subsystem I
  std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> fock(
      new std::vector<std::vector<MatrixInBasis<SCFMode>>>(this->_nSet));

  for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Return zero if not needed
  if (_eoPot == Options::KIN_EMBEDDING_MODES::NADD_FUNC || _eoPot == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION ||
      _eoPot == Options::KIN_EMBEDDING_MODES::NONE ||
      // For readCoupledExcitation:
      this->_lrscf[J]->getBasisController() == this->_lrscf[I]->getBasisController())
    return fock;

  // Calculate overlap Matrix between the different subsystems:
  auto& libint = Libint::getInstance();
  Eigen::MatrixXd SAB = libint.compute1eInts(libint2::Operator::overlap, this->_lrscf[J]->getBasisController(),
                                             this->_lrscf[I]->getBasisController());

  if (_eoPot == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
    // Calculate sigma vectors for each guess vector with pseudo density matrix
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& f = (*fock)[iSet][iGuess];
        auto& p = (*densityMatrices)[iSet][iGuess];
        for_spin(f, p) {
          f_spin = SAB * p_spin * SAB.transpose();
          f_spin *= _levelShiftParameter;
        };
      }
    }
  }
  else if (_eoPot == Options::KIN_EMBEDDING_MODES::HUZINAGA || _eoPot == Options::KIN_EMBEDDING_MODES::HOFFMANN) {
    // Calculates ABPotential
    // ToDo: Additional environment systems are not taken into account for the EO contribution
    std::vector<std::shared_ptr<SystemController>> environmentSystem;
    for (unsigned int iLRSCF = 0; iLRSCF < this->_lrscf.size(); iLRSCF++) {
      // Put systems to the environemnt with a different system controller than I
      if (this->_lrscf[I]->getSys() != this->_lrscf[iLRSCF]->getSys()) {
        // Makes sure that the environment system wasn't already added
        if (environmentSystem.size() > 0) {
          bool push_back = false;
          for (auto envSys : environmentSystem) {
            if (envSys != this->_lrscf[iLRSCF]->getSys())
              push_back = true;
          }
          if (push_back)
            environmentSystem.push_back(this->_lrscf[iLRSCF]->getSys());
        }
        else {
          environmentSystem.push_back(this->_lrscf[iLRSCF]->getSys());
        }
      }
    }
    for (auto& pas : this->_lrscf[I]->getEnvSystems()) {
      environmentSystem.push_back(pas);
    }
    auto abPotential = ABEmbeddedBundleFactory<SCFMode>::produce(
        this->_lrscf[I]->getSys(), this->_lrscf[J]->getBasisController(), this->_lrscf[J]->getSys()->getGeometry(),
        environmentSystem, std::make_shared<EmbeddingSettings>(this->_lrscf[I]->getLRSCFSettings().embedding), true);
    SPMatrix<SCFMode> F_AB = abPotential->getABMatrix();
    for_spin(F_AB) {
      F_AB_spin = F_AB_spin * 0.5;
    };
    // Calculate sigma vectors for each guess vector with pseudo density matrix
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& f = (*fock)[iSet][iGuess];
        auto& p = (*densityMatrices)[iSet][iGuess];
        for_spin(f, p, F_AB) {
          f_spin -= F_AB_spin * p_spin * SAB.transpose();
          f_spin -= SAB * p_spin * F_AB_spin.transpose();
        };
      }
    }
  }

  Timings::timeTaken("LRSCF -  Fock-like matrix: EO");
  return fock;
}

template class EOSigmaVector<Options::SCF_MODES::RESTRICTED>;
template class EOSigmaVector<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
