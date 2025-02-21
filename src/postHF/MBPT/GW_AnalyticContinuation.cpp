/**
 * @file GW_AnalyticContinuation.cpp
 *
 * @date Apr 14, 2020
 * @author Johannes Toelle
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
#include "postHF/MBPT/GW_AnalyticContinuation.h"
/* Include Serenity Internal Headers */
#include "grid/GaussLegendre.h"
#include "math/diis/DIIS.h"
#include "math/linearAlgebra/PadeApproximation.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
namespace Serenity {

template<Options::SCF_MODES SCFMode>
GW_AnalyticContinuation<SCFMode>::GW_AnalyticContinuation(std::shared_ptr<LRSCFController<SCFMode>> lrscf,
                                                          GWTaskSettings settings,
                                                          std::vector<std::shared_ptr<SystemController>> envSystemController,
                                                          std::shared_ptr<RIIntegrals<SCFMode>> riInts, int startOrb, int endOrb)
  : MBPT<SCFMode>(lrscf, settings, envSystemController, riInts, startOrb, endOrb) {
  auto integrator = GaussLegendre(this->_settings.padePoints);
  _nodesPade = Eigen::VectorXcd(integrator.getGridPoints());
  _nodesPade.real().setZero();
  Eigen::VectorXd nodes = integrator.getGridPoints();
  _nodesPade.imag() = (1.0 + nodes.array()) / (1.0 - nodes.array());
}

template<Options::SCF_MODES SCFMode>
void GW_AnalyticContinuation<SCFMode>::calculateGWOrbitalenergies(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& vxc_energies,
                                                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& x_energies,
                                                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation,
                                                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& dsigma_de,
                                                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& z) {
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED && this->_env.size() > 0) {
    throw SerenityError("Unrestricted GW-AC with environmental screening not supported, yet!");
  }
  // Check if calculation if evGW
  unsigned int evGWCycles = 1;
  if (this->_settings.evGW) {
    evGWCycles = this->_settings.evGWcycles;
  }
  // DIIS initialization
  std::shared_ptr<DIIS> diis_qpiterations = nullptr;
  std::shared_ptr<DIIS> diis_evGW = nullptr;
  if (this->_settings.diis) {
    if (this->_settings.qpiterations > 0)
      diis_qpiterations = std::make_shared<DIIS>(this->_settings.diisMaxStore);
    if (this->_settings.evGW)
      diis_evGW = std::make_shared<DIIS>(this->_settings.diisMaxStore);
  };
  // the initial orbital energies
  auto initialOrbEig = this->_orbEig;
  // evGW iterations
  for (unsigned int iGW = 0; iGW < evGWCycles; iGW++) {
    auto oldEvGWQP = qpEnergy;
    // calculate and set new quantities due to new Orbenergies in second interation
    if (iGW > 0) {
      this->setQPEner(qpEnergy);
      auto eia = this->calculateEia(qpEnergy, this->_nOcc, this->_nVirt);
      this->setEia(eia);
    }
    SpinPolarizedData<SCFMode, Eigen::MatrixXd> wnm;
    // fermi level
    if (this->_settings.ltconv == 0) {
      wnm = this->calculateWnmComplex();
    }
    else {
      wnm = this->calculateWnmComplexLT();
    }
    // QP interations
    for (unsigned int iQP = 0; iQP < this->_settings.qpiterations + 1; iQP++) {
      // set old qp energies
      SpinPolarizedData<SCFMode, Eigen::VectorXd> oldQP = qpEnergy;
      // do the actual analytic continuation
      Timings::takeTime("MBPT - Analytic Continuation");
      this->analyticContinuation(qpEnergy, correlation, wnm, dsigma_de, z);
      Timings::timeTaken("MBPT - Analytic Continuation");
      // Set new correlation energy, qpEnergies etc.
      for_spin(correlation, qpEnergy, x_energies, vxc_energies, initialOrbEig, z) {
        if (!this->_settings.linearized) {
          qpEnergy_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) =
              initialOrbEig_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
              x_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
              correlation_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) -
              vxc_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb);
        }
        else {
          qpEnergy_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) =
              initialOrbEig_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
              z_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb)
                  .cwiseProduct(x_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
                                correlation_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) -
                                vxc_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb));
        }
      };
      // Check for convergence and use convergence accelaration
      if (this->_settings.qpiterations > 0) {
        bool converged_qp = this->convergenCheck(qpEnergy, oldQP, diis_qpiterations, false, iQP);
        if (converged_qp)
          break;
      }
    }
    // Shift not inlcuded orbital energies by gap
    if (this->_settings.gap) {
      this->shiftByGap(qpEnergy, initialOrbEig);
    }
    // Check for convergence and use convergence accelaration
    if (this->_settings.evGW) {
      bool converged_evGW = this->convergenCheck(qpEnergy, oldEvGWQP, diis_evGW, true, iGW);
      if (converged_evGW)
        break;
    }
  } /* Loop evGW iterations */
}

template<Options::SCF_MODES SCFMode>
void GW_AnalyticContinuation<SCFMode>::analyticContinuation(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                                            SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation,
                                                            SpinPolarizedData<SCFMode, Eigen::MatrixXd>& wnm,
                                                            SpinPolarizedData<SCFMode, Eigen::VectorXd>& dsigma_de,
                                                            SpinPolarizedData<SCFMode, Eigen::VectorXd>& z) {
  auto& nOcc = this->_nOcc;
  auto& nVirt = this->_nVirt;
  auto& orbEig = this->_orbEig;
  const double overpi = (1.0 / (2.0 * PI));
  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned iState = unsigned(this->_startOrb); iState < unsigned(this->_endOrb); iState++) {
    for_spin(qpEnergy, orbEig, wnm, correlation, nOcc, nVirt, dsigma_de, z) {
      int states = nOcc_spin + nVirt_spin;
      double fermi = (orbEig_spin(nOcc_spin - 1) + orbEig_spin(nOcc_spin)) / 2;
      double factor = (iState < nOcc_spin) ? 1.0 : -1.0;
      fermi += factor * this->_settings.fermiShift * EV_TO_HARTREE;
      Eigen::VectorXcd frequencies = _nodesPade;
      Eigen::VectorXcd funcValues = Eigen::VectorXcd::Zero(this->_settings.padePoints);
      for (unsigned int iPade = 0; iPade < this->_settings.padePoints; iPade++) {
        for (int mState = 0; mState < states; mState++) {
          for (unsigned int iFreq = 0; iFreq < this->_settings.integrationPoints; iFreq++) {
            std::complex<double> k_nm =
                1.0 / (frequencies(iPade) - std::complex<double>(0.0, this->_nodes(iFreq)) + fermi - orbEig_spin(mState));
            k_nm += 1.0 / (frequencies(iPade) + std::complex<double>(0.0, this->_nodes(iFreq)) + fermi - orbEig_spin(mState));
            k_nm *= overpi * this->_weights(iFreq);
            funcValues(iPade) -= k_nm * wnm_spin((iState - this->_startOrb) * states + mState, iFreq);
          }
        }
      }
      auto padeApprox = PadeApproximation(frequencies, funcValues);
      correlation_spin(iState) =
          (padeApprox.padeApproximation(std::complex<double>(qpEnergy_spin(iState) - fermi, this->_settings.eta))).real();
      // Numerical Derivative
      if (this->_settings.linearized) {
        auto positive =
            padeApprox.padeApproximation(std::complex<double>(qpEnergy_spin(iState) + this->_settings.derivativeShift, 0.0));
        auto negative =
            padeApprox.padeApproximation(std::complex<double>(qpEnergy_spin(iState) - this->_settings.derivativeShift, 0.0));
        dsigma_de_spin(iState) = (positive.real() - negative.real()) / (2.0 * this->_settings.derivativeShift);
        // Reevaluate with imaginary shift
        if (std::abs(dsigma_de_spin(iState)) > 1.0) {
          correlation_spin(iState) =
              (padeApprox.padeApproximation(std::complex<double>(qpEnergy_spin(iState), this->_settings.imagShift))).real();
          positive = padeApprox.padeApproximation(
              std::complex<double>(qpEnergy_spin(iState) + this->_settings.derivativeShift, this->_settings.imagShift));
          negative = padeApprox.padeApproximation(
              std::complex<double>(qpEnergy_spin(iState) - this->_settings.derivativeShift, this->_settings.imagShift));
          dsigma_de_spin(iState) = (positive.real() - negative.real()) / (2.0 * this->_settings.derivativeShift);
        }
        double temp = 1.0 / (1.0 - dsigma_de_spin(iState));
        if (temp > 1.0 || temp < 0.0) {
          temp = 1.0;
        }
        z_spin(iState) = temp;
      }
    };
  }
  Eigen::setNbThreads(0);
}

template class GW_AnalyticContinuation<Options::SCF_MODES::RESTRICTED>;
template class GW_AnalyticContinuation<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
