/**
 * @file SystemAdditionTask.cpp
 *
 * @author Moritz Bensberg
 * @date Jan 9, 2020
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
#include "tasks/SystemAdditionTask.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h" //Basis set labels.
#include "basis/BasisController.h"             //Supersystem basis construction.
#include "basis/BasisFunctionMapper.h"         //Coefficient matrix resorting.
#include "data/ElectronicStructure.h"          //Definition of an electronic structure.
#include "data/OrbitalController.h"            //Definition of an orbital controller.
#include "data/matrices/CoefficientMatrix.h"   //Coefficient matrix definition.
#include "geometry/Geometry.h"                 //Build super. geometry.
#include "integrals/wrappers/Libint.h"         //A--B overlap matrix for projection.
#include "io/FormattedOutput.h"                //Captions.
#include "io/FormattedOutputStream.h"          //Filtered output streams.
#include "misc/SerenityError.h"                //Error messages.
#include "misc/SystemSplittingTools.h"         //Searching for atoms.
#include "system/SystemController.h"           //System controller definition.
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
SystemAdditionTask<SCFMode>::SystemAdditionTask(std::shared_ptr<SystemController> supersystem,
                                                std::vector<std::shared_ptr<SystemController>> subsystems)
  : _supersystem(supersystem), _subsystems(subsystems) {
}

template<Options::SCF_MODES SCFMode>
void SystemAdditionTask<SCFMode>::checkGeometry(std::shared_ptr<Geometry> supersystemGeometry) {
  unsigned int nSuperSystemAtoms = _supersystem->getGeometry()->getNAtoms();
  unsigned int sumOverSubsystemAtoms = supersystemGeometry->getNAtoms();
  if (nSuperSystemAtoms != sumOverSubsystemAtoms)
    throw SerenityError("ERROR: The number of supersystem and subsystem atoms do not match!");

  auto subsystemAtoms = supersystemGeometry->getAtoms();
  for (const auto& subAtom : subsystemAtoms) {
    unsigned int supersystemIndex = SystemSplittingTools<SCFMode>::matchAtom(_supersystem->getGeometry(), subAtom);
    if (supersystemIndex > nSuperSystemAtoms)
      throw SerenityError((std::string) "ERROR: The subsystem atoms do not add up to the supersystem atoms");
  } // for subAtom
}

template<Options::SCF_MODES SCFMode>
void SystemAdditionTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, {_supersystem}, _subsystems);
  printSubSectionTitle("Adding Subsystem up to a Supersystem");
  /*
   * 1. Build joined geometry.
   * 2. Build supersystem basis.
   * 3. Build supersystem coefficient matrix by sparse projection of occupied orbitals.
   * 4. Set electronic structure of supersystem.
   */
  OutputControl::nOut << "The supersystem > " << _supersystem->getSystemName()
                      << " < is constructed from the subsystems" << std::endl;
  auto supersystemGeometry = std::make_shared<Geometry>();
  int totalCharge = 0;
  int totalSpin = 0;
  for (unsigned int i = 0; i < _subsystems.size(); ++i) {
    *supersystemGeometry += *_subsystems[i]->getGeometry();
    totalCharge += _subsystems[i]->getCharge();
    totalSpin += _subsystems[i]->getSpin();
    OutputControl::nOut << "    " << _subsystems[i]->getSystemName() << std::endl;
  }
  supersystemGeometry->deleteIdenticalAtoms();
  if (settings.checkSuperGeom)
    checkGeometry(supersystemGeometry);
  if (settings.checkSuperCharge) {
    if (totalCharge != _supersystem->getCharge() || totalSpin != _supersystem->getSpin())
      throw SerenityError("ERROR: The subsystem charges/spins do not add up to the supersystem charge/spin");
  }
  // Delete old geometry.
  while (_supersystem->getGeometry()->getNAtoms() > 0)
    _supersystem->getGeometry()->deleteAtom(0);
  // Add new geometry.
  *_supersystem->getGeometry() += *supersystemGeometry;
  // Reset the basis
  _supersystem->setBasisController(nullptr);
  _supersystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_COULOMB);
  _supersystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_CORREL);
  // Set new charge/spin
  _supersystem->setSpin(totalSpin);
  _supersystem->setCharge(totalCharge);
  // Get the occupied orbitals.
  if (settings.addOccupiedOrbitals) {
    auto supersystemBasisController = _supersystem->getBasisController();
    unsigned int nSuperBasisFunc = supersystemBasisController->getNBasisFunctions();
    auto newCoefficientMatrixPtr =
        std::unique_ptr<CoefficientMatrix<SCFMode>>(new CoefficientMatrix<SCFMode>(supersystemBasisController));
    auto newEigenvaluesPtr = std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(
        new SpinPolarizedData<SCFMode, Eigen::VectorXd>(Eigen::VectorXd::Zero(nSuperBasisFunc)));
    auto newCoreOrbitalPtr =
        std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>(Eigen::VectorXi::Zero(nSuperBasisFunc));
    // Tools for creating sorting or projection matrices, which are needed in to
    // transform the subsystem MOs to the supersystem basis.
    BasisFunctionMapper basisFuncMapperSubToSuper(supersystemBasisController);
    auto& libint = Libint::getInstance();
    Eigen::MatrixXd overlapSuper =
        libint.compute1eInts(LIBINT_OPERATOR::overlap, supersystemBasisController, supersystemBasisController);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(overlapSuper, Eigen::ComputeThinU | Eigen::ComputeThinV);
    SpinPolarizedData<SCFMode, unsigned int> nOccSuper(0);
    SpinPolarizedData<SCFMode, unsigned int> nOccSuper2(0);
    for (auto subsystem : _subsystems) {
      CoefficientMatrix<SCFMode>& newCoefficientMatrix = *newCoefficientMatrixPtr;
      SpinPolarizedData<SCFMode, Eigen::VectorXd>& newEigenvalues = *newEigenvaluesPtr;
      SpinPolarizedData<SCFMode, Eigen::VectorXi>& newCoreOrbitals = *newCoreOrbitalPtr;
      CoefficientMatrix<SCFMode> subsystemCoefficients =
          subsystem->template getActiveOrbitalController<SCFMode>()->getCoefficients();
      SpinPolarizedData<SCFMode, Eigen::VectorXd> subsystemEigenvalues =
          subsystem->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
      SpinPolarizedData<SCFMode, Eigen::VectorXi> subsystemCoreOrbitals =
          subsystem->template getActiveOrbitalController<SCFMode>()->getOrbitalFlags();
      auto nOccSub = subsystem->template getNOccupiedOrbitals<SCFMode>();
      Eigen::MatrixXd projection;
      if (subsystem->getAtomCenteredBasisController()->getBasisLabel() ==
          _supersystem->getAtomCenteredBasisController()->getBasisLabel()) {
        // The subsystem-basis set will be at least a subset of the supersystem-basis set.
        // Resorting will do the job.
        projection = *basisFuncMapperSubToSuper.getSparseProjection(subsystem->getBasisController());
      }
      else {
        // The orbital coefficients need to be projected into the new basis-set.
        // Performs least square fit.
        OutputControl::dOut << "NOTE: The basis-set labels of the supersystem is different from the basis-set label"
                            << std::endl;
        OutputControl::dOut << "      of subsystem " << subsystem->getSystemName() << std::endl;
        OutputControl::dOut << "      The subsystem orbitals will be projected to the supersystem basis-set." << std::endl;
        OutputControl::dOut << "      Note that this projection may be inaccurate." << std::endl;
        Eigen::MatrixXd overlapSubSuper =
            libint.compute1eInts(LIBINT_OPERATOR::overlap, subsystem->getBasisController(), supersystemBasisController);
        projection = svd.solve(overlapSubSuper).transpose();
      }
      // Get the new coefficients and eigenvalues from the subsystem orbitals.
      for_spin(subsystemCoefficients, subsystemEigenvalues, newCoefficientMatrix, newEigenvalues, nOccSuper, nOccSub) {
        newCoefficientMatrix_spin.block(0, nOccSuper_spin, nSuperBasisFunc, nOccSub_spin) =
            projection.transpose() * subsystemCoefficients_spin.leftCols(nOccSub_spin);
        newEigenvalues_spin.segment(nOccSuper_spin, nOccSub_spin) = subsystemEigenvalues_spin.head(nOccSub_spin);
        nOccSuper_spin += nOccSub_spin;
      };
      for_spin(subsystemCoreOrbitals, newCoreOrbitals, nOccSuper2, nOccSub) {
        newCoreOrbitals_spin.segment(nOccSuper2_spin, nOccSub_spin) = subsystemCoreOrbitals_spin.head(nOccSub_spin);
        nOccSuper2_spin += nOccSub_spin;
      };
    } // for iSub
    // Build the final electronic structure.
    auto newOrbitalController =
        std::make_shared<OrbitalController<SCFMode>>(std::move(newCoefficientMatrixPtr), supersystemBasisController,
                                                     std::move(newEigenvaluesPtr), std::move(newCoreOrbitalPtr));
    auto newElectronicStructure = std::make_shared<ElectronicStructure<SCFMode>>(
        newOrbitalController, _supersystem->getOneElectronIntegralController(), nOccSuper);
    _supersystem->setElectronicStructure<SCFMode>(newElectronicStructure);
    newElectronicStructure->toHDF5(_supersystem->getHDF5BaseName(), _supersystem->getSystemIdentifier());
  } // if settings.addOccupiedOrbitals
  _supersystem->getGeometry()->printToFile(_supersystem->getHDF5BaseName(), _supersystem->getSystemIdentifier());
}

template class SystemAdditionTask<Options::SCF_MODES::RESTRICTED>;
template class SystemAdditionTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
