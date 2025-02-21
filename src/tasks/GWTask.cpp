/**
 * @file   GWTask.cpp
 *
 * @date   24.03.2020
 * @author J. Toelle
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

/* Include Class Header*/
#include "tasks/GWTask.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "postHF/MBPT/GW_Analytic.h"
#include "postHF/MBPT/GW_AnalyticContinuation.h"
#include "postHF/MBPT/GW_ContourDeformation.h"
#include "potentials/ExchangePotential.h"
#include "potentials/FuncPotential.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"
/* Include Std and External Headers */
#include <math.h>
#include <Eigen/Dense>
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GWTask<SCFMode>::GWTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
                        const std::vector<std::shared_ptr<SystemController>>& passiveSystems)
  : _act(activeSystems), _env(passiveSystems) {
}

template<Options::SCF_MODES SCFMode>
void GWTask<SCFMode>::run() {
  this->avoidMixedSCFModes(SCFMode, _act);
  if (_env.size() == 0 && settings.environmentScreening) {
    settings.environmentScreening = false;
  }
  // Create LRSCFController
  LRSCFTaskSettings settingsLRSCF;
  settingsLRSCF.nafThresh = settings.nafThresh;
  settingsLRSCF.frozenCore = settings.frozenCore;
  settingsLRSCF.coreOnly = settings.coreOnly;
  settingsLRSCF.densFitCache = settings.densFitCache;
  auto lrscfContr = std::make_shared<LRSCFController<SCFMode>>(_act[0], settingsLRSCF);

  auto C_old = lrscfContr->getCoefficients();
  auto E_old = lrscfContr->getEigenvalues();

  // Apply frozen core
  if (settings.frozenCore) {
    lrscfContr->applyFrozenCore();
  }
  else if (settings.coreOnly) {
    lrscfContr->applyCoreOnly();
  }

  if (settings.mbpttype == Options::MBPT::GW) {
    printSectionTitle("GW Module");

    if (_act.size() > 1) {
      throw SerenityError("Only one active subsystem supported!");
    }
    auto orb = _act[0]->template getElectronicStructure<SCFMode>()->getMolecularOrbitals();
    // Initializing for GW calculation
    auto coeffs = lrscfContr->getCoefficients();
    auto orbEigenValues = lrscfContr->getEigenvalues();
    const auto nOcc = lrscfContr->getNOccupied();
    // Use truncated virtual orbital space if it was manipulated before
    const auto nVirt = lrscfContr->getNVirtual();
    // Initialize the quantities
    SpinPolarizedData<SCFMode, Eigen::VectorXd> qpEnergy = orbEigenValues;
    SpinPolarizedData<SCFMode, Eigen::VectorXd> correlation = orbEigenValues;
    SpinPolarizedData<SCFMode, Eigen::VectorXd> dsigma_de = orbEigenValues;
    SpinPolarizedData<SCFMode, Eigen::VectorXd> lienarizazion_z = orbEigenValues;
    // States included
    int endOrb = 0;
    int startOrb = 0;
    for_spin(nOcc, nVirt, dsigma_de, lienarizazion_z, correlation) {
      dsigma_de_spin.setZero();
      lienarizazion_z_spin.setZero();
      correlation_spin.setZero();
      if (settings.nVirt > nVirt_spin)
        settings.nVirt = nVirt_spin;
      if (settings.nOcc > nOcc_spin)
        settings.nOcc = nOcc_spin;
      endOrb = nOcc_spin + settings.nVirt;
      startOrb = nOcc_spin - settings.nOcc;
    };
    // atoms of all subsystems
    std::vector<std::shared_ptr<Atom>> superSystemAtoms;
    superSystemAtoms.insert(superSystemAtoms.end(), _act[0]->getAtoms().begin(), _act[0]->getAtoms().end());
    // if environment systems are given. Evaluate the xc contribution on grid including env subsystems
    if (_env.size() > 0) {
      for (auto sys : _env) {
        for (auto atom : sys->getAtoms()) {
          for (auto check : _act[0]->getAtoms()) {
            double dist = distance(*atom, *check);
            if (dist < settings.gridCutOff) {
              superSystemAtoms.push_back(atom);
              break;
            }
          }
        }
      }
    }
    auto superSystemGeometry = std::make_shared<Geometry>(superSystemAtoms);
    superSystemGeometry->deleteIdenticalAtoms();
    // supersystem grid
    std::shared_ptr<GridController> supersystemgrid =
        AtomCenteredGridControllerFactory::produce(superSystemGeometry, _act[0]->getSettings().grid);
    _act[0]->setGridController(supersystemgrid);
    // XC Func
    auto functional = _act[0]->getSettings().customFunc.basicFunctionals.size()
                          ? Functional(_act[0]->getSettings().customFunc)
                          : resolveFunctional(_act[0]->getSettings().dft.functional);
    if (functional.isRSHybrid())
      throw SerenityError(" RS-Hybrid Functional not supported for GW calculations! ");
    std::shared_ptr<FuncPotential<SCFMode>> Vxc(new FuncPotential<SCFMode>(
        _act[0], _act[0]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        _act[0]->getGridController(), functional));
    // Exact exchange potential
    if (functional.getHfExchangeRatio() > 0.0)
      std::cout << " WARNING: Exchange energy contribution is already modified by exact exchange of XC-Functional!\n"
                << std::endl;
    std::shared_ptr<ExchangePotential<SCFMode>> exchange(new ExchangePotential<SCFMode>(
        _act[0], _act[0]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        1.0 - functional.getHfExchangeRatio(), _act[0]->getSettings().basis.integralThreshold,
        _act[0]->getSettings().basis.integralIncrementThresholdStart,
        _act[0]->getSettings().basis.integralIncrementThresholdEnd, _act[0]->getSettings().basis.incrementalSteps));

    auto pot_vxc = Vxc->getMatrix();
    auto pot_x = exchange->getMatrix();
    // Delete the potential
    Vxc = nullptr;
    exchange = nullptr;

    SpinPolarizedData<SCFMode, Eigen::VectorXd> vxc_energies;
    SpinPolarizedData<SCFMode, Eigen::VectorXd> x_energies;
    for_spin(coeffs, pot_vxc, pot_x, vxc_energies, x_energies) {
      vxc_energies_spin.resize(coeffs_spin.cols());
      x_energies_spin.resize(coeffs_spin.cols());
      for (unsigned int i = 0; i < coeffs_spin.cols(); i++) {
        vxc_energies_spin(i) = coeffs_spin.col(i).transpose() * pot_vxc_spin * coeffs_spin.col(i);
        x_energies_spin(i) = coeffs_spin.col(i).transpose() * pot_x_spin * coeffs_spin.col(i);
      }
      // Delete the potential
      pot_vxc_spin.resize(0, 0);
      pot_x_spin.resize(0, 0);
    };
    std::cout << std::fixed << std::setprecision(9);
    if (settings.gwtype == Options::GWALGORITHM::ANALYTIC) {
      printSubSectionTitle("Analytic-GW calculation");
      // Some information output
      std::cout << " Number of states included : " << endOrb - startOrb << std::endl;
      std::cout << " Eta                       : " << settings.eta << std::endl;
      std::cout << " QP-Iterations             : " << settings.qpiterations << std::endl;
      if (settings.qpiterations > 0) {
        std::cout << " Convergence Threshold     : " << settings.ConvergenceThreshold << std::endl;
        if (settings.diis)
          std::cout << " DIIS vectors stored       : " << settings.diisMaxStore << std::endl;
      }
      std::cout << " ------------------------------------------------ \n" << std::endl;
      auto gw_analytic = std::make_shared<GW_Analytic<SCFMode>>(lrscfContr, settings, _env, nullptr, startOrb, endOrb);
      gw_analytic->calculateGWOrbitalenergies(qpEnergy, vxc_energies, x_energies, correlation, dsigma_de, lienarizazion_z);
    }
    if (settings.gwtype == Options::GWALGORITHM::CD || settings.gwtype == Options::GWALGORITHM::AC) {
      if (settings.gwtype == Options::GWALGORITHM::CD)
        printSubSectionTitle("Contour Deformation (CD)-GW calculation");
      if (settings.gwtype == Options::GWALGORITHM::AC)
        printSubSectionTitle("Analytic Continuation (AC)-GW calculation");

      // Some information output
      std::cout << " Number of states included : " << endOrb - startOrb << std::endl;
      std::cout << " Eta                       : " << settings.eta << std::endl;
      std::cout << " QP-Iterations             : " << settings.qpiterations << std::endl;
      if (settings.evGW)
        std::cout << " evGW cycles               : " << settings.evGWcycles << std::endl;
      if (settings.qpiterations > 0 || settings.evGW) {
        std::cout << " Convergence Threshold     : " << settings.ConvergenceThreshold << std::endl;
        if (settings.diis)
          std::cout << " DIIS vectors stored       : " << settings.diisMaxStore << std::endl;
      }
      std::cout << " Integration points        : " << settings.integrationPoints << std::endl;
      if (settings.gwtype == Options::GWALGORITHM::AC) {
        std::cout << " Pade points               : " << settings.padePoints << std::endl;
        std::cout << " Fermi shift               : " << settings.fermiShift << std::endl;
        if (settings.linearized) {
          std::cout << " Derivative shift          : " << settings.derivativeShift << std::endl;
          std::cout << " Imaginary shift           : " << settings.imagShift << std::endl;
        }
      }
      std::cout << " ------------------------------------------------ \n" << std::endl;
      auto riInts = std::make_shared<RIIntegrals<SCFMode>>(lrscfContr, LIBINT_OPERATOR::coulomb, 0.0, true, startOrb,
                                                           endOrb, this->superMolGeo());
      if (settings.gwtype == Options::GWALGORITHM::CD) {
        auto gw_cd = std::make_shared<GW_ContourDeformation<SCFMode>>(lrscfContr, settings, _env, riInts, startOrb, endOrb);
        gw_cd->calculateGWOrbitalenergies(qpEnergy, vxc_energies, x_energies, correlation);
      }
      else if (settings.gwtype == Options::GWALGORITHM::AC) {
        auto gw_ac = std::make_shared<GW_AnalyticContinuation<SCFMode>>(lrscfContr, settings, _env, riInts, startOrb, endOrb);
        gw_ac->calculateGWOrbitalenergies(qpEnergy, vxc_energies, x_energies, correlation, dsigma_de, lienarizazion_z);
      }
    }
    _correlationEnergies = correlation;
    printGWResult(orbEigenValues, qpEnergy, vxc_energies, x_energies, correlation, lienarizazion_z, dsigma_de, startOrb, endOrb);

    if (settings.frozenCore) {
      unsigned nCore = 0;
      for (const auto& atom : _act[0]->getGeometry()->getAtoms()) {
        nCore += atom->getNCoreElectrons();
      }
      nCore /= 2;

      std::cout << "\n  ! Sorting back " << nCore << " core orbitals !" << std::endl;

      coeffs = C_old;
      for_spin(qpEnergy, E_old) {
        Eigen::VectorXd qp = qpEnergy_spin;
        unsigned len = qp.size();
        for (unsigned j = 0; j < len; ++j) {
          qpEnergy_spin(j) = (j < nCore) ? E_old_spin(j) : qp(j - nCore);
        }
      };
    }

    // Set new orbital energies
    orb->setDiskMode(false, _act[0]->getHDF5BaseName(), _act[0]->getSettings().identifier);
    orb->updateOrbitals(coeffs, qpEnergy);
    orb->toHDF5(_act[0]->getHDF5BaseName(), _act[0]->getSettings().identifier);
  }
  // RPA correlation energy calculation
  else if (settings.mbpttype == Options::MBPT::RPA) {
    printSectionTitle("RPA Module");
    // XC Func
    auto functional = _act[0]->getSettings().customFunc.basicFunctionals.size()
                          ? Functional(_act[0]->getSettings().customFunc)
                          : resolveFunctional(_act[0]->getSettings().dft.functional);
    std::shared_ptr<FuncPotential<SCFMode>> Vxc(new FuncPotential<SCFMode>(
        _act[0], _act[0]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        _act[0]->getGridController(), functional));
    // Exact exchange potential
    if (functional.getHfExchangeRatio() > 0.0)
      std::cout << " WARNING: Exchange energy contribution is already modified by exact exchange of XC-Functional!\n"
                << std::endl;
    std::shared_ptr<ExchangePotential<SCFMode>> exchange(new ExchangePotential<SCFMode>(
        _act[0], _act[0]->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
        1.0 - functional.getHfExchangeRatio(), _act[0]->getSettings().basis.integralThreshold,
        _act[0]->getSettings().basis.integralIncrementThresholdStart,
        _act[0]->getSettings().basis.integralIncrementThresholdEnd, _act[0]->getSettings().basis.incrementalSteps));

    auto ener_vxc = Vxc->getEnergy(
        _act[0]->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());
    auto ener_x = exchange->getEnergy(
        _act[0]->template getElectronicStructure<SCFMode>()->getDensityMatrixController()->getDensityMatrix());

    auto riInts =
        std::make_shared<RIIntegrals<SCFMode>>(lrscfContr, LIBINT_OPERATOR::coulomb, 0.0, true, 0, 0, this->superMolGeo());
    auto mbpt = std::make_shared<MBPT<SCFMode>>(MBPT<SCFMode>(lrscfContr, settings, _env, riInts, 0, 0));
    _correlationEnergy = mbpt->calculateRPACorrelationEnergy();
    std::cout << std::fixed << std::setprecision(9);
    std::cout << "\n Integration points        : " << settings.integrationPoints << std::endl;
    std::cout << " ------------------------------------------------ " << std::endl;
    std::cout << " Exchange Energy / Ha: " << ener_x << std::endl;
    std::cout << " Exchange Energy / eV: " << ener_x * HARTREE_TO_EV << std::endl;
    std::cout << " XC-Func. Energy/  Ha: " << ener_vxc << std::endl;
    std::cout << " XC-Func. Energy/  eV: " << ener_vxc * HARTREE_TO_EV << std::endl;
    std::cout << " ------------------------------------------------ " << std::endl;
    std::cout << " RPA Correlation Energy/ Ha: " << _correlationEnergy << std::endl;
    std::cout << " RPA Correlation Energy/ eV: " << _correlationEnergy * HARTREE_TO_EV << std::endl;
    std::cout << " ------------------------------------------------ " << std::endl;
    std::cout << " Total Correction Energy / Ha: " << ener_x + _correlationEnergy - ener_vxc << "\n" << std::endl;
  }
}

template<Options::SCF_MODES SCFMode>
void GWTask<SCFMode>::printGWResult(SpinPolarizedData<SCFMode, Eigen::VectorXd>& orbEigenValues,
                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& vxc_energies,
                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& x_energies,
                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation,
                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& z,
                                    SpinPolarizedData<SCFMode, Eigen::VectorXd>& dsigma_de, int& start, int& end) {
  std::cout << "\n  GW Summary:    " << std::endl;
  std::cout << "  ------------------- " << std::endl;
  // Write output
  printf("%4s %3s %12s %12s %12s %12s %12s %12s %8s %9s \n", "", " # ", " Orb. Ener. / eV ", " QP Ener. / eV ",
         "   Sigma / eV  ", "   Sigma_x / eV   ", "   Sigma_c / eV   ", "   Sigma_vxc / eV  ", "   Z   ", "   dS/dE   ");
  printf("%4s %3s %12s %12s %12s %12s %12s %12s %8s %9s \n", "", "---", "-----------------", "---------------",
         "---------------", "------------------", "------------------", "-------------------", "-------", "-----------");
  int spincounter = 0;
  for_spin(orbEigenValues, qpEnergy, vxc_energies, x_energies, dsigma_de, z, correlation) {
    if (spincounter == 1)
      std::cout << "     "
                   "---------------------------------------------------------------------------------------------------"
                   "---------------------------------"
                << std::endl;
    for (int i = start; i < end; i++) {
      printf("%4s %3d %+14.5f %+14.5f %1s %+14.5f %1s %+14.5f %2s %+14.5f %4s %+14.5f %3s %9.3f %+10.3f \n", "",
             (i + 1), orbEigenValues_spin(i) * HARTREE_TO_EV, qpEnergy_spin(i) * HARTREE_TO_EV, "",
             (x_energies_spin(i) + correlation_spin(i) - vxc_energies_spin(i)) * HARTREE_TO_EV, "",
             x_energies_spin(i) * HARTREE_TO_EV, "", correlation_spin(i) * HARTREE_TO_EV, "",
             vxc_energies_spin(i) * HARTREE_TO_EV, "", z_spin(i), dsigma_de_spin(i));
    }
    spincounter++;
  };
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Geometry> GWTask<SCFMode>::superMolGeo() {
  // Create supersystem geometry
  std::vector<std::shared_ptr<Atom>> superAtoms;
  auto sysAtoms = _act[0]->getAtoms();
  for (auto& atom : sysAtoms) {
    superAtoms.push_back(atom);
  }
  if (settings.subsystemAuxillaryBasisOnly == false && _env.size() > 0) {
    std::cout << " Supersystem auxillary basis is used!" << std::endl;
    for (auto& env : _env) {
      auto envAtoms = env->getAtoms();
      for (auto& atom : envAtoms) {
        superAtoms.push_back(atom);
      }
    }
  }
  auto superGeo = std::make_shared<Geometry>(superAtoms);
  superGeo->deleteIdenticalAtoms();
  // ToDo: Test if properly working; Carefull in case of PbE??
  superGeo->deleteGhostAtoms();
  return superGeo;
}

template class GWTask<Options::SCF_MODES::RESTRICTED>;
template class GWTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
