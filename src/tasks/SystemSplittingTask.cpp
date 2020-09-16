/**
 * @file SystemSplittingTask.cpp
 *
 * @author Moritz Bensberg
 * @date Jan 8, 2020
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
#include "tasks/SystemSplittingTask.h"
/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/SPADEAlgorithm.h"              //SPADE algorithm.
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h" //Population based separation.
#include "basis/AtomCenteredBasisController.h"                        //Population based separation and basis label.
#include "data/ElectronicStructure.h"                                 //Electronic structure split.
#include "data/OrbitalController.h"                                   //Coefficients.
#include "geometry/Geometry.h"                                        //Geometry and atoms.
#include "integrals/OneElectronIntegralController.h"                  //Overlap integrals for Mulliken populations.
#include "io/FormattedOutput.h"                                       //Captions.
#include "io/FormattedOutputStream.h"                                 //Filtered output streams.
#include "misc/SerenityError.h"                                       //Error messages.
#include "misc/SystemSplittingTools.h"                                //System partitioning.
#include "misc/WarningTracker.h"                                      //Warnings.
#include "settings/Settings.h"                                        //Settings.
/* Include Std and External Headers */
#include <iomanip> //setw(...) for ostream.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SystemSplittingTask<SCFMode>::SystemSplittingTask(std::shared_ptr<SystemController> supersystem,
                                                  std::vector<std::shared_ptr<SystemController>> subsystems)
  : _supersystem(supersystem), _subsystems(subsystems) {
}

template<Options::SCF_MODES SCFMode>
std::vector<unsigned int> SystemSplittingTask<SCFMode>::findAtoms(std::shared_ptr<SystemController> subsystem) {
  unsigned int nSupersystemAtoms = _supersystem->getGeometry()->getNAtoms();
  std::vector<unsigned int> subsystemAtomIndices;
  for (unsigned int iAtom = 0; iAtom < subsystem->getGeometry()->getNAtoms(); ++iAtom) {
    auto subAtom = subsystem->getGeometry()->getAtoms()[iAtom];
    if (subAtom->isDummy())
      continue;
    unsigned int indexInSupersystem = SystemSplittingTools<SCFMode>::matchAtom(_supersystem->getGeometry(), subAtom);
    if (indexInSupersystem >= nSupersystemAtoms)
      throw SerenityError((std::string) "ERROR: Subsystem does not constitute a sub part of the supersystem!\n" +
                          "       Subsystem " + subsystem->getSystemName() + " Atom index " + iAtom + " " +
                          subAtom->getAtomType()->getElementSymbol());
    subsystemAtomIndices.push_back(indexInSupersystem);
  }
  return subsystemAtomIndices;
}

template<Options::SCF_MODES SCFMode>
void SystemSplittingTask<SCFMode>::checkInput() {
  // Check charge/spin.
  switch (settings.systemPartitioning) {
    case Options::SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES: {
      int totCharge = _supersystem->getCharge();
      int totSpin = _supersystem->getSpin();
      for (const auto subsystem : _subsystems) {
        totCharge -= subsystem->getCharge();
        totSpin -= subsystem->getSpin();
      } // for subsystem
      if (totCharge != 0 || totSpin != 0)
        throw SerenityError("ERROR: Enforce charges was used but the subsystem charges/spins do not sum up to the "
                            "supersystem charge/spin!");
      printSubSectionTitle("Running Orbital Partitioning with Enforced Charges");
      break;
    }
    case Options::SYSTEM_SPLITTING_ALGORITHM::BEST_MATCH: {
      printSubSectionTitle("Running Orbital Partitioning with Best Match Assignment");
      break;
    }
    case Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD: {
      printSubSectionTitle("Running Orbital Partitioning with Population-based Assignment");
      break;
    }
    case Options::SYSTEM_SPLITTING_ALGORITHM::SPADE: {
      if (_subsystems.size() != 2)
        throw SerenityError("The SPADE algorithm does only work with two subsystems!");
      printSubSectionTitle("Running Orbital Partitioning using the SPADE Algorithm");
      break;
    }
    case Options::SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES: {
      if (_subsystems.size() != 2)
        throw SerenityError("The SPADE algorithm does only work with two subsystems!");
      printSubSectionTitle("Running Orbital Partitioning using the SPADE Algorithm with Enforced Charges");
      break;
    }
  }
  bool warningTriggered = false;
  for (auto subsystem : _subsystems) {
    if (subsystem->getSettings().scfMode != SCFMode) {
      std::string scfModeString = (SCFMode == RESTRICTED) ? "RESTRICTED" : "UNRESTRICTED";
      WarningTracker::printWarning((std::string) "    WARNING: The SCFMode of the subsystem " + subsystem->getSystemName() +
                                       "\n             is different from the SCFMode used in the task!\n" +
                                       "             The system's SCFMode will be changed to " + scfModeString,
                                   true);
      warningTriggered = true;
    }
  }
  if (warningTriggered)
    std::cout << std::endl;
}
template<>
SpinPolarizedData<RESTRICTED, std::string> SystemSplittingTask<RESTRICTED>::getOutputPraefixes() {
  return "";
}

template<>
SpinPolarizedData<UNRESTRICTED, std::string> SystemSplittingTask<UNRESTRICTED>::getOutputPraefixes() {
  SpinPolarizedData<UNRESTRICTED, std::string> toReturn;
  toReturn.alpha = "Alpha";
  toReturn.beta = "Beta";
  return toReturn;
}

template<Options::SCF_MODES SCFMode>
void SystemSplittingTask<SCFMode>::run() {
  checkInput();
  SpinPolarizedData<SCFMode, Eigen::VectorXi> assignment;
  unsigned int nSubsystems = _subsystems.size();
  const auto nOccSuper = _supersystem->getNOccupiedOrbitals<SCFMode>();
  /*
   * 1. Calculate orbital populations.
   * 2. Construct supersystem atom to subsystem atom mapping.
   * 3. Loop orbitals and assign system indices by population.
   * 4. Split the electronic structure and save it all on disk.
   */
  if (settings.systemPartitioning != Options::SYSTEM_SPLITTING_ALGORITHM::SPADE and
      settings.systemPartitioning != Options::SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES) {
    // Populations.
    auto orbitalWisePopulations = MullikenPopulationCalculator<SCFMode>::calculateAtomwiseOrbitalPopulations(
        _supersystem->getActiveOrbitalController<SCFMode>()->getCoefficients(),
        _supersystem->getOneElectronIntegralController()->getOverlapIntegrals(),
        _supersystem->getAtomCenteredBasisController()->getBasisIndices());
    // Atom mapping.
    std::vector<std::vector<unsigned int>> atomIndices;
    std::vector<bool> atomMatched(_supersystem->getGeometry()->getNAtoms(), false);
    for (const auto& subsystem : _subsystems) {
      auto indices = findAtoms(subsystem);
      for (const auto index : indices) {
        if (!atomMatched[index]) {
          atomMatched[index] = true;
        }
        else {
          throw SerenityError((std::string) "ERROR: There were duplicated atoms between the subsystems.\n" +
                              "       Subsystem: " + subsystem->getSystemName() + " Supersystem atom index " + index);
        }
      }
      atomIndices.push_back(indices);
    }
    OutputControl::nOut << "  Partitioning the supersystem > " << _supersystem->getSystemName()
                        << " < into the subsystems" << std::endl;
    for (unsigned int iSub = 0; iSub < _subsystems.size(); ++iSub) {
      auto subsystem = _subsystems[iSub];
      OutputControl::nOut << "    " << subsystem->getSystemName() << std::setw(12) << " Atom-indices: ";
      for (auto atomIndex : atomIndices[iSub])
        OutputControl::nOut << atomIndex << " ";
      OutputControl::nOut << std::endl;
    }

    SPMatrix<SCFMode> subsystemWisePopulations;
    // Try to generate a "best match" assignment of the orbitals based on their
    // predominant localization.
    for_spin(subsystemWisePopulations, orbitalWisePopulations, nOccSuper, assignment) {
      subsystemWisePopulations_spin = Eigen::MatrixXd::Zero(nSubsystems, nOccSuper_spin);
      assignment_spin = Eigen::VectorXi::Constant(nOccSuper_spin, nSubsystems);
      for (unsigned int iOrb = 0; iOrb < nOccSuper_spin; ++iOrb) {
        for (unsigned int iSub = 0; iSub < nSubsystems; ++iSub) {
          for (const auto subAtomIndex : atomIndices[iSub]) {
            subsystemWisePopulations_spin(iSub, iOrb) += std::fabs(orbitalWisePopulations_spin(subAtomIndex, iOrb));
          } // for subAtomIndex
        }   // for iSub
        unsigned int subMax;
        subsystemWisePopulations_spin.col(iOrb).maxCoeff(&subMax);
        assignment_spin(iOrb) = subMax;
        // Assign the orbitals based on a population threshold for the first (0th) subsystem.
        if (settings.systemPartitioning == Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD) {
          if (subsystemWisePopulations_spin(0, iOrb) >= settings.orbitalThreshold)
            assignment_spin(iOrb) = 0;
        } // if settings.systemPartitioning == Options::SYSTEM_SPLITTING_ALGORITHM::POPULATION_THRESHOLD
      }   // for iOrb
    };
    // Assign orbitals to the subsystems until charge and spin is correct.
    // The order of the assignment is based on the population analysis done previously.
    if (settings.systemPartitioning == Options::SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES) {
      for (unsigned int iSub = 0; iSub < nSubsystems; ++iSub) {
        const auto subsystem = _subsystems[iSub];
        const auto nOccSub = subsystem->getNOccupiedOrbitals<SCFMode>();
        for_spin(nOccSub, assignment, subsystemWisePopulations) {
          for (unsigned int i = 0; i < nOccSub_spin; ++i) {
            unsigned int maxOrb;
            subsystemWisePopulations_spin.row(iSub).maxCoeff(&maxOrb);
            assignment_spin(maxOrb) = iSub;
            subsystemWisePopulations_spin(iSub, maxOrb) = 0.0;
          } // for i
        };
      } // for iSub
    }
  } // if
  else {
    SPADEAlgorithm<SCFMode> spade(_supersystem, _subsystems[0]);
    assignment = spade.run();
    if (settings.systemPartitioning == Options::SYSTEM_SPLITTING_ALGORITHM::SPADE_ENFORCE_CHARGES) {
      auto nOccAct = _subsystems[0]->getNOccupiedOrbitals<SCFMode>();
      for_spin(assignment, nOccAct) {
        assignment_spin = Eigen::VectorXi::Constant(assignment_spin.size(), 1);
        assignment_spin.head(nOccAct_spin).setZero();
      };
    }
  }
  for (unsigned int iSub = 0; iSub < nSubsystems; ++iSub) {
    auto subsystem = _subsystems[iSub];
    // Generate the final orbital selection based on the criteria constructed above.
    SpinPolarizedData<SCFMode, std::vector<bool>> orbitalSelection;
    for_spin(nOccSuper, orbitalSelection, assignment) {
      orbitalSelection_spin = std::vector<bool>(nOccSuper_spin, false);
      for (unsigned int iOrb = 0; iOrb < nOccSuper_spin; ++iOrb) {
        if ((int)iSub == assignment_spin(iOrb))
          orbitalSelection_spin[iOrb] = true;
      } // for iOrb
    };
    // Update basis, electronic structure and geometry and save them on disk.
    subsystem->getGeometry()->addAsDummy(*_supersystem->getGeometry());
    subsystem->getGeometry()->deleteIdenticalAtoms();
    subsystem->setBasisController(nullptr);
    subsystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_COULOMB);
    subsystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_CORREL);

    // Select the part of the electronic structure.
    auto subsystemES = SystemSplittingTools<SCFMode>::splitElectronicStructure(_supersystem, orbitalSelection).first;

    if (subsystem->getAtomCenteredBasisController()->getBasisLabel() !=
        _supersystem->getAtomCenteredBasisController()->getBasisLabel()) {
      OutputControl::dOut << "NOTE: The basis label of the subsystem " << subsystem->getSystemName() << std::endl;
      OutputControl::dOut << "      is different from the basis label of the supersystem!" << std::endl;
      OutputControl::dOut << "      The density matrix will be projected into the subsystem basis." << std::endl;
      OutputControl::dOut << "      Orbitals will not be available until the next Fock-matrix diagonalization." << std::endl;
      OutputControl::dOut << "      Note that the projection may be inaccurate." << std::endl;
      const DensityMatrix<SCFMode>& selectedDensity = subsystemES->getDensityMatrix();
      DensityMatrix<SCFMode> projectedDensity =
          SystemSplittingTools<SCFMode>::projectMatrixIntoNewBasis(selectedDensity, subsystem->getBasisController());
      subsystemES = std::make_shared<ElectronicStructure<SCFMode>>(
          subsystem->getBasisController(), subsystem->getGeometry(), subsystemES->getNOccupiedOrbitals());
      subsystemES->getDensityMatrixController()->setDensityMatrix(projectedDensity);
    }
    else {
      // The same basis set is used for the subsystem and the supersystem.
      // The coefficients may need to be resorted, but no projection is necessary!
      subsystemES = SystemSplittingTools<SCFMode>::resortBasisSetOfElectronicStructure(
          subsystemES, subsystem->getBasisController(), subsystem->getOneElectronIntegralController());
    }

    subsystem->setSCFMode(SCFMode);
    subsystem->setElectronicStructure<SCFMode>(subsystemES);
    subsystem->getGeometry()->printToFile(subsystem->getHDF5BaseName(), subsystem->getSettings().identifier);
    subsystemES->toHDF5(subsystem->getHDF5BaseName(), subsystem->getSettings().identifier);
    // Check charges and spin of the system, since it may have changed during the partitioning.
    if (settings.systemPartitioning != Options::SYSTEM_SPLITTING_ALGORITHM::ENFORCE_CHARGES) {
      auto oldNElSub = subsystem->getNElectrons<SCFMode>();
      auto newNElSub = SystemSplittingTools<SCFMode>::getNElectrons(orbitalSelection).first;
      unsigned int nTotOld = 0;
      unsigned int nTotNew = 0;
      for_spin(oldNElSub, newNElSub) {
        nTotOld += oldNElSub_spin;
        nTotNew += newNElSub_spin;
      };
      int subSpin = SystemSplittingTools<SCFMode>::getSpin(newNElSub);
      unsigned int minusNElChanged = nTotOld - nTotNew;
      subsystem->setSpin(subSpin);
      subsystem->setCharge(subsystem->getCharge() + minusNElChanged);
    } // if !settings.enforceCharges
    // Print information about the system selection to the output.
    OutputControl::nOut << "--------------------------------------------------------" << std::endl;
    OutputControl::nOut << "  Partitioning for subsystem " << subsystem->getSystemName() << std::endl;
    OutputControl::nOut << "  Charge: " << subsystem->getCharge() << std::endl;
    OutputControl::nOut << "  Spin:   " << subsystem->getSpin() << std::endl;
    OutputControl::vOut << "Orbital selection:" << std::endl;
    auto praefixes = getOutputPraefixes();
    for_spin(orbitalSelection, praefixes) {
      if (SCFMode == UNRESTRICTED)
        OutputControl::vOut << praefixes_spin << " orbitals:" << std::endl;
      unsigned int printCounter = 0;
      for (unsigned int iOrb = 0; iOrb < orbitalSelection_spin.size(); ++iOrb) {
        if (orbitalSelection_spin[iOrb]) {
          OutputControl::vOut << iOrb << " ";
          ++printCounter;
          if (printCounter % 15 == 0)
            OutputControl::vOut << std::endl;
        } // if orbitalSelection_spin[iOrb]
      }   // for iOrb
      if (printCounter % 15 != 0)
        OutputControl::vOut << std::endl;
      OutputControl::vOut << "Number of occ. orbitals " << printCounter << std::endl;
    };
  } // for iSub
  OutputControl::nOut << "--------------------------------------------------------" << std::endl;
}

template class SystemSplittingTask<Options::SCF_MODES::RESTRICTED>;
template class SystemSplittingTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
