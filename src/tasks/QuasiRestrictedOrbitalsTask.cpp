/**
 * @file QuasiRestrictedOrbitalsTask.cpp
 *
 * @date Mai 10, 2021
 * @author Moritz Bensberg
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
#include "tasks/QuasiRestrictedOrbitalsTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"                //Access density matrix.
#include "data/OrbitalController.h"                  //Update orbitals.
#include "integrals/OneElectronIntegralController.h" //Overlap integrals.
#include "io/FormattedOutputStream.h"                //Filtered output streams.
#include "potentials/bundles/PotentialBundle.h"      //Update energy.
#include "scf/SCFAnalysis.h"                         //Calculate S^2.
#include "system/SystemController.h"                 //Access electronic structure.
#include "tasks/FDETask.h"                           //Energy calculation.
/* Include Std and External Headers */

namespace Serenity {

template<Options::SCF_MODES SCFMode>
QuasiRestrictedOrbitalsTask<SCFMode>::QuasiRestrictedOrbitalsTask(std::shared_ptr<SystemController> activeSystem,
                                                                  std::vector<std::shared_ptr<SystemController>> environmentSystems)
  : _activeSystem(activeSystem), _environmentSystems(environmentSystems) {
}

template<>
void QuasiRestrictedOrbitalsTask<RESTRICTED>::run() {
  printSectionTitle("Quasi-Restricted orbital construction");
  OutputControl::nOut << "  Reference orbitals are already restricted!" << std::endl;
  OutputControl::nOut << "  Nothing to be done here. Continuing." << std::endl;
  return;
}
template<>
Eigen::MatrixXd QuasiRestrictedOrbitalsTask<UNRESTRICTED>::calculateNaturalOrbitals() {
  MatrixInBasis<RESTRICTED> densityMatrix =
      _activeSystem->template getElectronicStructure<UNRESTRICTED>()->getDensityMatrix().total();
  /*
   * Chose an orthogonal basis to express the orbitals in. We will work with the alpha-coefficients.
   */
  const Eigen::MatrixXd alphaCoeff = _activeSystem->getActiveOrbitalController<UNRESTRICTED>()->getCoefficients().alpha;
  Eigen::MatrixXd alphaCoeff_inv =
      alphaCoeff.transpose() * _activeSystem->getOneElectronIntegralController()->getOverlapIntegrals();
  Eigen::MatrixXd tmp = alphaCoeff_inv * densityMatrix * alphaCoeff_inv.transpose();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(tmp);
  const Eigen::VectorXd& eigenvalues = eigenvalueSolver.eigenvalues();
  _occupationNumbers = eigenvalues;
  _nVirtuals = (_occupationNumbers.array() < 0.49).count();
  _nDOMOs = (_occupationNumbers.array() > 1.5).count();
  _nSOMOs = _occupationNumbers.size() - _nVirtuals - _nDOMOs;
  Eigen::MatrixXd eigenvectors = alphaCoeff * eigenvalueSolver.eigenvectors();
  unsigned int nCols = eigenvectors.cols();
  Eigen::MatrixXd naturalOrbitals = Eigen::MatrixXd::Zero(nCols, nCols);
  for (unsigned int iCol = 0; iCol < nCols; ++iCol) {
    const unsigned int index = nCols - 1 - iCol;
    naturalOrbitals.col(index) = eigenvectors.col(iCol);
    _occupationNumbers(index) = eigenvalues(iCol);
  }
  _initialized = true;
  return naturalOrbitals;
}
template<>
unsigned int QuasiRestrictedOrbitalsTask<UNRESTRICTED>::getNDOMOs() {
  if (not _initialized)
    calculateNaturalOrbitals();
  return _nDOMOs;
}
template<>
unsigned int QuasiRestrictedOrbitalsTask<UNRESTRICTED>::getNSOMOs() {
  if (not _initialized)
    calculateNaturalOrbitals();
  return _nSOMOs;
}
template<>
unsigned int QuasiRestrictedOrbitalsTask<UNRESTRICTED>::getNVirtuals() {
  if (not _initialized)
    calculateNaturalOrbitals();
  return _nVirtuals;
}
template<>
Eigen::VectorXd QuasiRestrictedOrbitalsTask<UNRESTRICTED>::getOccupationNumbers() {
  if (not _initialized)
    calculateNaturalOrbitals();
  return _occupationNumbers;
}
template<>
void QuasiRestrictedOrbitalsTask<UNRESTRICTED>::updateEnergy() {
  if (_environmentSystems.size() > 0) {
    FDETask<UNRESTRICTED> fde(_activeSystem, _environmentSystems);
    fde.settings.embedding = settings.embedding;
    fde.settings.skipSCF = true;
    fde.run();
  }
  auto potBundle = _activeSystem->template getElectronicStructure<UNRESTRICTED>()->getPotentialBundle();
  FockMatrix<UNRESTRICTED> fock = potBundle->getFockMatrix(
      _activeSystem->template getElectronicStructure<UNRESTRICTED>()->getDensityMatrix(),
      _activeSystem->template getElectronicStructure<UNRESTRICTED>()->getEnergyComponentController());
}

template<>
void QuasiRestrictedOrbitalsTask<UNRESTRICTED>::run() {
  printSectionTitle("Quasi-Restricted orbital construction");
  /*
   * Diagonalize UHF/UKS-density matrix:
   * A) MOs with occ. number n= 1.0: SOMOs.
   * B) MOs with occ. number n~2.0:  DOMOs (double occ. orbitals).
   * C) MOs with occ. number n~0.0:  virtual.
   * Canonicalization.
   * DOMOs: diag. F-beta.
   * virtual: diag. F-alpha.
   * SOMOs: diag. (F-alpha + F-beta)/2, two orbital energies:
   *        e-alpha = <SOMO|f-alpha|SOMO> and e-beta = <SOMO|f-beta|SOMO>
   */
  Eigen::MatrixXd eigenvectors = this->calculateNaturalOrbitals();

  OutputControl::nOut << std::string(100, '-') << std::endl;
  OutputControl::nOut << "  Number of DOMOs:       " << this->getNDOMOs() << std::endl;
  OutputControl::nOut << "  Number of SOMOs:       " << this->getNSOMOs() << std::endl;
  OutputControl::nOut << "  Number of virtual MOs: " << this->getNVirtuals() << std::endl;
  {
    SCFAnalysis<UNRESTRICTED> scfAnalysis({_activeSystem});
    double oldS2 = scfAnalysis.getS2();
    double S = fabs(0.5 * _activeSystem->getSpin());
    double targetS2 = S * (S + 1);
    OutputControl::nOut << "  Initial <S*S>:  " << oldS2 << "  (should be " << targetS2 << ")" << std::endl;
  }
  OutputControl::dOut << "Occupations" << std::endl;
  for (unsigned int iOrb = 0; iOrb < _occupationNumbers.size(); ++iOrb)
    OutputControl::dOut << _occupationNumbers(iOrb) << std::endl;
  updateToNonCanonicalNaturalOrbitals(eigenvectors);
  this->updateEnergy();
  if (this->settings.canonicalize) {
    canonicalizeOrbitals(_nDOMOs, _nSOMOs, _nVirtuals);
  }
  {
    SCFAnalysis<UNRESTRICTED> scfAnalysis({_activeSystem});
    double newS2 = scfAnalysis.getS2();
    OutputControl::nOut << "  Final   <S*S>:  " << newS2 << std::endl;
  }
  OutputControl::nOut << std::string(100, '-') << std::endl;
  _activeSystem->getElectronicStructure<UNRESTRICTED>()->getEnergyComponentController()->printAllComponents();
  return;
}

template<Options::SCF_MODES SCFMode>
void QuasiRestrictedOrbitalsTask<SCFMode>::updateToNonCanonicalNaturalOrbitals(const Eigen::MatrixXd& naturalOrbitals) {
  CoefficientMatrix<UNRESTRICTED> newCoefficients(_activeSystem->getBasisController());
  newCoefficients.alpha = naturalOrbitals;
  newCoefficients.beta = naturalOrbitals;
  auto orbitalController = _activeSystem->template getActiveOrbitalController<UNRESTRICTED>();
  orbitalController->updateOrbitals(newCoefficients, orbitalController->getEigenvalues());
}

template<Options::SCF_MODES SCFMode>
void QuasiRestrictedOrbitalsTask<SCFMode>::canonicalizeOrbitals(unsigned int nDOMOs, unsigned int nSOMOs, unsigned int nVirt) {
  auto potBundle = _activeSystem->template getElectronicStructure<UNRESTRICTED>()->getPotentialBundle();
  FockMatrix<UNRESTRICTED> fock = potBundle->getFockMatrix(
      _activeSystem->template getElectronicStructure<UNRESTRICTED>()->getDensityMatrix(),
      _activeSystem->template getElectronicStructure<UNRESTRICTED>()->getEnergyComponentController());
  auto orbitalController = _activeSystem->template getActiveOrbitalController<UNRESTRICTED>();
  const auto& coefficients = orbitalController->getCoefficients();

  CoefficientMatrix<UNRESTRICTED> newCoefficients(_activeSystem->getBasisController());
  SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd> eigenvalues(Eigen::VectorXd::Zero(coefficients.cols()));
  // DOMOs
  if (nDOMOs > 0) {
    Eigen::MatrixXd c_domo = coefficients.alpha.leftCols(nDOMOs);
    const Eigen::MatrixXd f_domo = c_domo.transpose() * fock.beta * c_domo;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eDOMO(f_domo);
    c_domo = c_domo * eDOMO.eigenvectors();
    newCoefficients.alpha.leftCols(nDOMOs) = c_domo;
    newCoefficients.beta.leftCols(nDOMOs) = c_domo;
    eigenvalues.alpha.head(nDOMOs) = eDOMO.eigenvalues();
    eigenvalues.beta.head(nDOMOs) = eDOMO.eigenvalues();
  }
  // SOMOs
  if (nSOMOs > 0) {
    Eigen::MatrixXd c_somo = coefficients.alpha.middleCols(nDOMOs, nSOMOs);
    const Eigen::MatrixXd f_somo = c_somo.transpose() * (0.5 * (fock.beta + fock.alpha)) * c_somo;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eSOMO(f_somo);
    c_somo = c_somo * eSOMO.eigenvectors();
    Eigen::VectorXd epsSOMOalpha = (c_somo.transpose() * fock.alpha * c_somo).diagonal();
    Eigen::VectorXd epsSOMObeta = (c_somo.transpose() * fock.beta * c_somo).diagonal();
    eigenvalues.alpha.segment(nDOMOs, nSOMOs) = epsSOMOalpha;
    eigenvalues.beta.segment(nDOMOs, nSOMOs) = epsSOMObeta;
    newCoefficients.alpha.middleCols(nDOMOs, nSOMOs) = c_somo;
    newCoefficients.beta.middleCols(nDOMOs, nSOMOs) = c_somo;
  }
  // Virtuals
  if (nVirt > 0) {
    Eigen::MatrixXd c_virt = coefficients.alpha.rightCols(nVirt);
    const Eigen::MatrixXd f_virt = c_virt.transpose() * fock.alpha * c_virt;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eVirt(f_virt);
    c_virt = c_virt * eVirt.eigenvectors();
    newCoefficients.alpha.rightCols(nVirt) = c_virt;
    newCoefficients.beta.rightCols(nVirt) = c_virt;
    eigenvalues.alpha.tail(nVirt) = eVirt.eigenvalues();
    eigenvalues.beta.tail(nVirt) = eVirt.eigenvalues();
  }

  orbitalController->updateOrbitals(newCoefficients, eigenvalues);
}

template class QuasiRestrictedOrbitalsTask<Options::SCF_MODES::RESTRICTED>;
template class QuasiRestrictedOrbitalsTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
