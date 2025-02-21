/**
 * @file   FiniteFieldTask.cpp
 *
 * @date   Apr 29, 2021
 * @author Patrick Eschenbach, Niklas Niemeyer
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
#include "tasks/FiniteFieldTask.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrix.h"
#include "geometry/Atom.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutput.h"
#include "io/FormattedOutputStream.h"
#include "parameters/Constants.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FDETask.h"
#include "tasks/LRSCFTask.h"
#include "tasks/ScfTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FiniteFieldTask<SCFMode>::FiniteFieldTask(const std::shared_ptr<SystemController>& activeSystem,
                                          const std::vector<std::shared_ptr<SystemController>>& environmentSystems)
  : _activeSystem(activeSystem),
    _environmentSystems(environmentSystems),
    _polarizability(Eigen::Matrix3d::Zero(3, 3)),
    _hyperPolarizability(std::vector<Eigen::Matrix3d>(3, Eigen::Matrix3d::Zero(3, 3))),
    _isotropicPolarizability(0.0) {
  _embeddingScheme = Options::EMBEDDING_SCHEME::ISOLATED;
}

template<Options::SCF_MODES SCFMode>
void FiniteFieldTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _activeSystem);
  // run the task
  printSubSectionTitle("Finite Field Task");
  printBigCaption("Calculation for " + _activeSystem->getSystemName());
  if (_environmentSystems.size() > 0) {
    _embeddingScheme = Options::EMBEDDING_SCHEME::FDE;
  }
  bool rsHybrid = false;
  auto func = _activeSystem->getSettings().customFunc.basicFunctionals.size()
                  ? Functional(_activeSystem->getSettings().customFunc)
                  : resolveFunctional(_activeSystem->getSettings().dft.functional);
  if (func.isRSHybrid() && _activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    rsHybrid = true;
  }
  /* calculate dipolmoments for unperturbed active system */
  Eigen::MatrixXd unperturbedDipMoment = this->perturbedSCF(0, 0, BASE_PROPERTY::DIPOLE_MOMENT, 0.0);

  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  if (rsHybrid) {
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::erf_coulomb, 0, 4);
  }
  Eigen::Matrix<Eigen::MatrixXd, 3, 4> baseProperty;

  /* purely numerical calculation for static frequency */
  if (settings.frequency == 0.0) {
    /* calculate dipole moments for all finite fields */
    for (unsigned iCoord = 0; iCoord < 3; ++iCoord) {
      baseProperty(iCoord, 0) =
          this->perturbedSCF(iCoord, +1.0 * settings.finiteFieldStrength, BASE_PROPERTY::DIPOLE_MOMENT, settings.frequency);
      baseProperty(iCoord, 1) =
          this->perturbedSCF(iCoord, -1.0 * settings.finiteFieldStrength, BASE_PROPERTY::DIPOLE_MOMENT, settings.frequency);
      baseProperty(iCoord, 2) =
          this->perturbedSCF(iCoord, +2.0 * settings.finiteFieldStrength, BASE_PROPERTY::DIPOLE_MOMENT, settings.frequency);
      baseProperty(iCoord, 3) =
          this->perturbedSCF(iCoord, -2.0 * settings.finiteFieldStrength, BASE_PROPERTY::DIPOLE_MOMENT, settings.frequency);
    }
    /* Polarizability */
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        _polarizability(i, j) = (2.0 / 3.0 * (baseProperty(j, 0)(i) - baseProperty(j, 1)(i)) -
                                 1.0 / 12.0 * (baseProperty(j, 2)(i) - baseProperty(j, 3)(i))) /
                                settings.finiteFieldStrength;
      }
    }
    _isotropicPolarizability = (1.0 / 3.0) * _polarizability.trace();
    /* Hyperpolarizability */
    double fieldStrengthSquared = settings.finiteFieldStrength * settings.finiteFieldStrength;
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        double b =
            ((baseProperty(i, 2)(j) + baseProperty(i, 3)(j)) / 3.0 - (baseProperty(i, 0)(j) + baseProperty(i, 1)(j)) / 3.0) /
            (fieldStrengthSquared);
        _hyperPolarizability[i](i, j) = b;
        _hyperPolarizability[i](j, i) = b;
        _hyperPolarizability[j](i, i) = b;
      }
    }
    _hyperPolarizability[0](1, 2) = std::numeric_limits<double>::infinity();
    _hyperPolarizability[0](2, 1) = std::numeric_limits<double>::infinity();
    _hyperPolarizability[1](2, 0) = std::numeric_limits<double>::infinity();
    _hyperPolarizability[1](0, 2) = std::numeric_limits<double>::infinity();
    _hyperPolarizability[2](0, 1) = std::numeric_limits<double>::infinity();
    _hyperPolarizability[2](1, 0) = std::numeric_limits<double>::infinity();
  }
  /* analytical calculation for dynamic frequency via the LRSCFTask */
  else {
    /* calculate polarizabilities for all finite fields */
    _polarizability = this->perturbedSCF(0, 0, BASE_PROPERTY::DYNAMIC_POLARIZABILITY, settings.frequency);
    _isotropicPolarizability = (1.0 / 3.0) * _polarizability.trace();
    baseProperty.resize(3, 2);
    if (settings.hyperPolarizability) {
      for (unsigned iCoord = 0; iCoord < 3; ++iCoord) {
        baseProperty(iCoord, 0) = this->perturbedSCF(iCoord, +1.0 * settings.finiteFieldStrength,
                                                     BASE_PROPERTY::DYNAMIC_POLARIZABILITY, settings.frequency);
        baseProperty(iCoord, 1) = this->perturbedSCF(iCoord, -1.0 * settings.finiteFieldStrength,
                                                     BASE_PROPERTY::DYNAMIC_POLARIZABILITY, settings.frequency);
      }
      for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {
          for (unsigned k = 0; k < 3; ++k) {
            _hyperPolarizability[k](i, j) =
                (baseProperty(k, 0)(i, j) - baseProperty(k, 1)(i, j)) / (2.0 * settings.finiteFieldStrength);
          }
        }
      }
    }
  }
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  if (rsHybrid) {
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::erf_coulomb, 0, 4);
  }

  // Reset system.
  _activeSystem->setElectronicStructure<SCFMode>(nullptr);
  _activeSystem->setBasisController(nullptr);
  _activeSystem->setElectricField({0.0, 0.0, 0.0}, 0.0, true, false);

  printBigCaption("Finite-Field Task Results");

  if (settings.frequency == 0) {
    printf("  Note that these results were obtained by numerical differentiation.\n\n");
  }
  else {
    printf("  Note that the polarizability tensor was obtained analytically via the LRSCFTask.\n\n");
    if (settings.hyperPolarizability) {
      printf("  The hyperpolarizability was obtained by numerical differentiation\n");
      printf("  of the analytical dynamic polarizability with respect to a static electric field.\n\n");
    }
  }

  printf("  %-30s %12.3e \n", "Frequency / a.u.", settings.frequency * EV_TO_HARTREE);
  printf("  %-30s %12.3e \n", "Frequency / eV", settings.frequency);
  printf("  %-30s %12.3f \n\n\n", "Wavelength / nm", HARTREE_TO_NM / (settings.frequency * EV_TO_HARTREE));
  printSmallCaption("Dipole Moment Vector / a.u.");
  printf("%4s %12s %12s %12s \n", "", "x", "y", "z");
  printf("%4s %12.8f %12.8f %12.8f \n\n", "", unperturbedDipMoment(0), unperturbedDipMoment(1), unperturbedDipMoment(2));
  OutputControl::nOut << std::endl;
  OutputControl::nOut << std::endl;

  std::vector<std::string> coord({"x", "y", "z"});
  printSmallCaption("Polarizability Tensor / a.u.");
  printf("%4s %12s %12s %12s \n", "", "x", "y", "z");
  for (unsigned iRow = 0; iRow < _polarizability.rows(); ++iRow)
    printf(" %4s %12.8f %12.8f %12.8f \n", coord[iRow].c_str(), _polarizability.row(iRow)(0),
           _polarizability.row(iRow)(1), _polarizability.row(iRow)(2));
  OutputControl::nOut << std::endl;
  printf("%-33s %12.8f\n", "  Isotropic Polarizability / a.u.:", _isotropicPolarizability);
  OutputControl::nOut << std::endl;
  OutputControl::nOut << std::endl;
  if (settings.hyperPolarizability) {
    printSmallCaption("Hyperpolarizability Tensor / a.u.");
    for (unsigned iRow = 0; iRow < _hyperPolarizability.size(); ++iRow) {
      printf("%2s\n", coord[iRow].c_str());
      printf("%4s %12s %12s %12s \n", "", "x", "y", "z");
      for (unsigned iCoord = 0; iCoord < _hyperPolarizability.size(); ++iCoord) {
        printf(" %4s %12.8f %12.8f %12.8f \n", coord[iCoord].c_str(), _hyperPolarizability[iRow](iCoord, 0),
               _hyperPolarizability[iRow](iCoord, 1), _hyperPolarizability[iRow](iCoord, 2));
      }
      OutputControl::nOut << std::endl;
    }
  }
} /* enf of run function */

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd FiniteFieldTask<SCFMode>::perturbedSCF(unsigned direction, double fStrength, BASE_PROPERTY propertyType,
                                                       double frequency) {
  if (direction > 2) {
    throw SerenityError("You are trying to manipulate a vector element that does not exist!");
  }
  Eigen::MatrixXd property;
  std::vector<double> position = {0.0, 0.0, 0.0};
  position[direction] = 1.0;
  _activeSystem->setElectricField(position, fStrength, true, true);
  switch (propertyType) {
      /* dipole moment */
    case BASE_PROPERTY::DIPOLE_MOMENT: {
      if (_embeddingScheme == Options::EMBEDDING_SCHEME::ISOLATED) {
        ScfTask<SCFMode> scf(_activeSystem);
        scf.run();
      }
      else if (_embeddingScheme == Options::EMBEDDING_SCHEME::FDE) {
        FDETask<SCFMode> fde(_activeSystem, _environmentSystems);
        fde.settings.embedding = settings.embedding;
        fde.run();
      }
      else {
        throw SerenityError("Unknown embedding scheme! Allowed: ISOLATED, FDE");
      }
      auto dipoleLength = _activeSystem->getOneElectronIntegralController()->getDipoleLengths();
      auto densityMatrix = _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrix().total();
      property = Eigen::MatrixXd::Zero(3, 1);
      for (unsigned j = 0; j < 3; j++) {
        /* electronic */
        property(j) = (densityMatrix.transpose() * dipoleLength[j]).trace();
        /* nuclear */
        for (unsigned int A = 0; A < _activeSystem->getNAtoms(); A++) {
          property(j) += _activeSystem->getAtoms()[A]->coords()(j) * _activeSystem->getAtoms()[A]->getNuclearCharge();
        }
      }
      break;
    }
    case BASE_PROPERTY::DYNAMIC_POLARIZABILITY: {
      /* polarizability */
      if (_embeddingScheme == Options::EMBEDDING_SCHEME::ISOLATED) {
        ScfTask<SCFMode> scf(_activeSystem);
        scf.run();
        LRSCFTask<SCFMode> lrscf({_activeSystem});
        lrscf.settings.frequencies = {frequency};
        lrscf.settings.densFitJ = _activeSystem->getSettings().basis.densFitJ;
        lrscf.settings.densFitK = _activeSystem->getSettings().basis.densFitK;
        lrscf.settings.densFitLRK = _activeSystem->getSettings().basis.densFitLRK;
        lrscf.settings.nEigen = 0;
        lrscf.run();
        property = std::get<1>(lrscf.getProperties()[0]);
      }
      else if (_embeddingScheme == Options::EMBEDDING_SCHEME::FDE) {
        FDETask<SCFMode> fde(_activeSystem, _environmentSystems);
        fde.settings.embedding = settings.embedding;
        fde.run();
        LRSCFTask<SCFMode> lrscf({_activeSystem}, _environmentSystems);
        lrscf.settings.embedding = settings.embedding;
        lrscf.settings.frequencies = {frequency};
        lrscf.settings.densFitJ = _activeSystem->getSettings().basis.densFitJ;
        lrscf.settings.densFitK = _activeSystem->getSettings().basis.densFitK;
        lrscf.settings.densFitLRK = _activeSystem->getSettings().basis.densFitLRK;
        lrscf.settings.nEigen = 0;
        lrscf.run();
        property = std::get<1>(lrscf.getProperties()[0]);
      }
      else {
        throw SerenityError("Unknown embedding scheme! Allowed: ISOLATED, FDE");
      }
      break;
    }
    default:
      throw SerenityError("You requested a non-existing property!");
      break;
  }
  return property;
}

template<Options::SCF_MODES SCFMode>
void FiniteFieldTask<SCFMode>::visit(FiniteFieldTaskSettings& c, set_visitor v, std::string blockname) {
  if (!blockname.compare("")) {
    visit_each(c, v);
    return;
  }
  // If reached, the blockname is unknown.
  if (c.embedding.visitAsBlockSettings(v, blockname))
    return;
  throw SerenityError((std::string) "Unknown block in FiniteFieldTaskSettings: " + blockname);
}

template<Options::SCF_MODES SCFMode>
Eigen::Matrix3d FiniteFieldTask<SCFMode>::getPolarizability() {
  return _polarizability;
}

template<Options::SCF_MODES SCFMode>
std::vector<Eigen::Matrix3d> FiniteFieldTask<SCFMode>::getHyperPolarizability() {
  return _hyperPolarizability;
}

template<Options::SCF_MODES SCFMode>
double FiniteFieldTask<SCFMode>::getIsotropicPolarizability() {
  return _isotropicPolarizability;
}

template class FiniteFieldTask<Options::SCF_MODES::RESTRICTED>;
template class FiniteFieldTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
