/**
 * @file EOSigmavector.cpp
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

/* Include Class Header*/
#include "postHF/LRSCF/Sigmavectors/EOSigmavector.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <iomanip>
/* Construction of embedded fock matrix. */
#include "basis/BasisFunctionMapper.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "potentials/ABFockMatrixConstruction/ABBundles/ABEmbeddedBundleFactory.h"
#include "potentials/ECPInteractionPotential.h"
#include "potentials/HuzinagaProjectionPotential.h"
#include "potentials/LevelshiftHybridPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "potentials/ZeroPotential.h"
#include "potentials/bundles/ESIPotentials.h"
#include "potentials/bundles/FDEPotentialBundleFactory.h"
#include "potentials/bundles/FDEPotentials.h"
#include "potentials/bundles/PBEPotentials.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
EOSigmavector<SCFMode>::EOSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                      std::vector<Eigen::MatrixXd> b, const EmbeddingSettings& embeddingSettings)
  : Sigmavector<SCFMode>(lrscf, b),
    _levelShiftParameter(embeddingSettings.levelShiftParameter),
    _eoPot(embeddingSettings.embeddingMode),
    _fermiShift(embeddingSettings.fermiShift) {
}

template<Options::SCF_MODES SCFMode>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>>
EOSigmavector<SCFMode>::calcF(unsigned int I, unsigned int J,
                              std::unique_ptr<std::vector<std::vector<MatrixInBasis<SCFMode>>>> densityMatrices) {
  Timings::takeTime("LRSCF -   Sigmavector:     EO");

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
  if (_eoPot == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION)
    return fock;
  // Only when the two different orbital transitions spaces are different
  if (I != J) {
    // Calculate overlap Matrix between the different subsystems:
    auto& libint = Libint::getInstance();
    Eigen::MatrixXd SAB = libint.compute1eInts(LIBINT_OPERATOR::overlap, this->_lrscf[J]->getBasisController(),
                                               this->_lrscf[I]->getBasisController());
    // FockMatrix
    SPMatrix<SCFMode> F_AB;
    // For list input
    if (_eoPot == Options::KIN_EMBEDDING_MODES::NONE &&
        this->_lrscf[I]->getLRSCFSettings().embedding.embeddingModeList.size() == 0)
      return fock;
    if (this->_lrscf[I]->getLRSCFSettings().embedding.embeddingModeList.size() > 0) {
      if (this->_lrscf[I]->getLRSCFSettings().embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::NADD_FUNC ||
          this->_lrscf[J]->getLRSCFSettings().embedding.embeddingModeList[J] == Options::KIN_EMBEDDING_MODES::NADD_FUNC) {
        return fock;
      }
      else if (this->_lrscf[I]->getLRSCFSettings().embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
        _eoPot = this->_lrscf[I]->getLRSCFSettings().embedding.embeddingModeList[I];
      }
      else if (this->_lrscf[I]->getLRSCFSettings().embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA) {
        throw SerenityError(
            "Huzinaga/Hoffmann kernel contribution for mixed exact/approximate embedding not yet implemented!");
      }
    }
    // Check if two subsystems have the same density
    bool sameDensity = false;
    if (this->_lrscf[I]->getLRSCFSettings().samedensity.size() > 0) {
      if (this->_lrscf[I]->getLRSCFSettings().samedensity[I] == this->_lrscf[I]->getLRSCFSettings().samedensity[J])
        sameDensity = true;
    }
    // Evaluation or reading of FAB if needed
    if (sameDensity || _eoPot == Options::KIN_EMBEDDING_MODES::HUZINAGA || _eoPot == Options::KIN_EMBEDDING_MODES::HOFFMANN ||
        _eoPot == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA) {
      auto activeSystem = this->_lrscf[I]->getSys();
      auto es_active = activeSystem->template getElectronicStructure<SCFMode>();
      // Try to use Fock matrix in electronic structure
      if (es_active->checkFock() && sameDensity) {
        std::cout << "   Fock matrix from a previous embedding/supermolecular calculation used!" << std::endl;
        F_AB = es_active->getFockMatrix();
      }
      else {
        std::vector<std::shared_ptr<SystemController>> environmentSystem;
        if (this->_lrscf[I]->getLRSCFSettings().samedensity.size() > 0) {
          std::vector<std::shared_ptr<SystemController>> supersystems;
          for (unsigned int iLRSCF = 0; iLRSCF < this->_lrscf.size(); iLRSCF++) {
            supersystems.push_back(this->_lrscf[iLRSCF]->getSys());
          }
          for (auto& pas : this->_lrscf[I]->getEnvSystems()) {
            supersystems.push_back(pas);
          }
          unsigned int iRefDens = this->_lrscf[I]->getLRSCFSettings().samedensity[I];
          unsigned int counter = 0;
          std::vector<unsigned int> included;
          for (auto i : supersystems) {
            if (this->_lrscf[I]->getLRSCFSettings().samedensity[counter] != iRefDens) {
              if (included.size() == 0) {
                included.push_back(this->_lrscf[I]->getLRSCFSettings().samedensity[counter]);
                environmentSystem.push_back(i);
              }
              else {
                if (std::find(included.begin(), included.end(), this->_lrscf[I]->getLRSCFSettings().samedensity[counter]) ==
                    included.end()) {
                  environmentSystem.push_back(i);
                }
              }
            }
            counter++;
          }
        }
        else {
          for (unsigned int iLRSCF = 0; iLRSCF < this->_lrscf.size(); iLRSCF++) {
            if (iLRSCF != I) {
              environmentSystem.push_back(this->_lrscf[iLRSCF]->getSys());
            }
          }
          for (auto& pas : this->_lrscf[I]->getEnvSystems()) {
            environmentSystem.push_back(pas);
          }
        }
        // Calculates ABPotential
        auto abPotential = ABEmbeddedBundleFactory<SCFMode>::produce(
            this->_lrscf[I]->getSys(), this->_lrscf[J]->getBasisController(), this->_lrscf[J]->getSys()->getGeometry(),
            {environmentSystem}, std::make_shared<EmbeddingSettings>(this->_lrscf[I]->getLRSCFSettings().embedding), true);
        F_AB = abPotential->getABMatrix();
      }
    }
    for (unsigned int iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned int iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        auto& f = (*fock)[iSet][iGuess];
        auto& p = (*densityMatrices)[iSet][iGuess];
        for_spin(f, p, F_AB) {
          if (_eoPot == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
            if (!sameDensity) {
              f_spin += _levelShiftParameter * SAB * p_spin * SAB.transpose();
            }
            else {
              f_spin += SAB * p_spin * F_AB_spin.transpose();
            }
          }
          else if (_eoPot == Options::KIN_EMBEDDING_MODES::HUZINAGA || _eoPot == Options::KIN_EMBEDDING_MODES::HOFFMANN) {
            if (sameDensity) {
              f_spin += SAB * p_spin * F_AB_spin.transpose();
            }
            else {
              f_spin -= F_AB_spin * p_spin * SAB.transpose();
            }
          }
          else if (_eoPot == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA) {
            if (sameDensity) {
              f_spin += SAB * p_spin * (F_AB_spin.transpose() - 2.0 * _fermiShift * SAB.transpose());
            }
            else {
              f_spin -= (F_AB_spin - 2.0 * _fermiShift * SAB) * p_spin * SAB.transpose();
            }
          }
        };
      }
    }
  } /*I!=J*/
  Timings::timeTaken("LRSCF -   Sigmavector:     EO");
  return fock;
}
template class EOSigmavector<Options::SCF_MODES::RESTRICTED>;
template class EOSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */