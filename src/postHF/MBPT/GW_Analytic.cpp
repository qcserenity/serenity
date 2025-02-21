/**
 * @file GW_Analytic.cpp
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
#include "postHF/MBPT/GW_Analytic.h"
/* Include Serenity Internal Headers */
#include "math/diis/DIIS.h"
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Sigmavectors/CoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/RICoulombSigmavector.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GW_Analytic<SCFMode>::GW_Analytic(std::shared_ptr<LRSCFController<SCFMode>> lrscf, GWTaskSettings settings,
                                  std::vector<std::shared_ptr<SystemController>> envSystemController,
                                  std::shared_ptr<RIIntegrals<SCFMode>> riInts, int startOrb, int endOrb)
  : MBPT<SCFMode>(lrscf, settings, envSystemController, riInts, startOrb, endOrb) {
}

template<Options::SCF_MODES SCFMode>
void GW_Analytic<SCFMode>::calculateGWOrbitalenergies(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                                      SpinPolarizedData<SCFMode, Eigen::VectorXd>& vxc_energies,
                                                      SpinPolarizedData<SCFMode, Eigen::VectorXd>& x_energies,
                                                      SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation,
                                                      SpinPolarizedData<SCFMode, Eigen::VectorXd>& dsigma_de,
                                                      SpinPolarizedData<SCFMode, Eigen::VectorXd>& z) {
  auto nocc = this->_nOcc;
  auto nvirt = this->_nVirt;
  auto orbEig = this->_orbEig;
  auto lrscf = this->_lrscfController;
  // Coeffs
  auto coeffs = lrscf->getCoefficients();
  // Set loading type for calculation
  Options::LRSCF_TYPE load_type = Options::LRSCF_TYPE::ISOLATED;
  if (this->_env.size() == 0) {
    load_type = Options::LRSCF_TYPE::ISOLATED;
  }
  else if (this->_env.size() >= 1) {
    load_type = Options::LRSCF_TYPE::UNCOUPLED;
  }
  // Add everything to sigma vectors
  auto eigenvectors = lrscf->getExcitationVectors(load_type);
  // obtain X and Y from loaded (X+Y) and (X-Y).
  // X = 0.5 * {(X+Y) + (X-Y)}.
  (*eigenvectors)[0] = 0.5 * ((*eigenvectors)[0] + (*eigenvectors)[1]);
  // Y = -{(X-Y) - X}.
  (*eigenvectors)[1] = -((*eigenvectors)[1] - (*eigenvectors)[0]);

  auto eigenvalues = lrscf->getExcitationEnergies(load_type);
  // Coulomb
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> temp;
  temp.push_back(lrscf);

  std::vector<Eigen::MatrixXd> xpm = {(*eigenvectors)[0] + (*eigenvectors)[1]};
  eigenvectors = nullptr;

  std::shared_ptr<Sigmavector<SCFMode>> coulomb;
  if (lrscf->getSysSettings().basis.densFitCorr == Options::DENS_FITS::RI) {
    coulomb = std::make_shared<RICoulombSigmavector<SCFMode>>(temp, xpm);
  }
  else {
    coulomb = std::make_shared<CoulombSigmavector<SCFMode>>(temp, xpm);
  }
  // Get perturbed Fock Matrix
  auto perturbedFockMatrices = coulomb->getPerturbedFockMatrix();

  // Initialize everything if frequency scan should be performed
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> frequency_dependency;
  Eigen::VectorXd frequencies;
  if (this->_settings.freq.size() == 3) {
    unsigned int nSteps = std::abs(this->_settings.freq[1] - this->_settings.freq[0]) / this->_settings.freq[2];
    for_spin(frequency_dependency) {
      frequency_dependency_spin = Eigen::MatrixXd::Zero(nSteps, this->_endOrb - this->_startOrb);
    };
    frequencies = Eigen::VectorXd::Zero(nSteps);
    for (unsigned int iStep = 0; iStep < nSteps; iStep++) {
      frequencies(iStep) = (this->_settings.freq[0] + iStep * this->_settings.freq[2]) * (1.0 / HARTREE_TO_EV);
    }
  }
  // DIIS initialization
  std::shared_ptr<DIIS> diis_qpiterations = nullptr;
  if (this->_settings.diis) {
    if (this->_settings.qpiterations > 0)
      diis_qpiterations = std::make_shared<DIIS>(this->_settings.diisMaxStore);
  };
  // QP interations
  for (unsigned int iQP = 0; iQP < this->_settings.qpiterations + 1; iQP++) {
    // set old qp energies
    SpinPolarizedData<SCFMode, Eigen::VectorXd> oldQP = qpEnergy;
    for_spin(correlation) {
      correlation_spin.setZero();
    };
    // Calculation of correlation self energy and derivative
    for (unsigned int iState = 0; iState < perturbedFockMatrices.size(); iState++) {
      auto& perturbedFock = perturbedFockMatrices[iState];
      for (int nState = this->_startOrb; nState < this->_endOrb; nState++) {
        // MO-Fock transformation
        SpinPolarizedData<SCFMode, Eigen::VectorXd> moFock;
        for_spin(nocc, nvirt, perturbedFock, coeffs, moFock) {
          // Transform perturbed fock matrix
          moFock_spin =
              coeffs_spin.col(nState).transpose() * perturbedFock_spin * coeffs_spin.leftCols(nocc_spin + nvirt_spin);
          moFock_spin = moFock_spin.cwiseProduct(moFock_spin);
        };
        for_spin(nocc, nvirt, orbEig, moFock, correlation, dsigma_de, frequency_dependency, qpEnergy) {
          double freq = qpEnergy_spin(nState);
          auto d_Vector = this->calculatedVector(freq, orbEig_spin, (*eigenvalues)(iState), nocc_spin, nvirt_spin);
          correlation_spin(nState) += (moFock_spin.cwiseProduct(d_Vector)).sum();
          // if the correlation self-energy should be calculated for mutliple frequencies
          if (this->_settings.freq.size() == 3) {
            for (unsigned int i = 0; i < frequencies.size(); i++) {
              auto d_freq = this->calculatedVector(frequencies(i), orbEig_spin, (*eigenvalues)(iState), nocc_spin, nvirt_spin);
              frequency_dependency_spin(i, nState - this->_startOrb) += (moFock_spin.cwiseProduct(d_freq)).sum();
            }
          }
          if (this->_settings.linearized) {
            Eigen::VectorXd d2_Vector = Eigen::VectorXd::Zero(nocc_spin + nvirt_spin);
            // Occupied
            d2_Vector.segment(0, nocc_spin) =
                pow(this->_settings.eta, 2) -
                (freq - orbEig_spin.segment(0, nocc_spin).array() + (*eigenvalues)(iState)).square();
            d2_Vector.segment(0, nocc_spin) =
                d2_Vector.segment(0, nocc_spin).array() /
                ((freq - orbEig_spin.segment(0, nocc_spin).array() + (*eigenvalues)(iState)).square() +
                 pow(this->_settings.eta, 2))
                    .square();
            // Virtual
            d2_Vector.segment(nocc_spin, nvirt_spin) =
                pow(this->_settings.eta, 2) -
                (freq - orbEig_spin.segment(nocc_spin, nvirt_spin).array() - (*eigenvalues)(iState)).square();
            d2_Vector.segment(nocc_spin, nvirt_spin) =
                d2_Vector.segment(nocc_spin, nvirt_spin).array() /
                ((freq - orbEig_spin.segment(nocc_spin, nvirt_spin).array() - (*eigenvalues)(iState)).square() +
                 pow(this->_settings.eta, 2))
                    .square();
            dsigma_de_spin(nState) += (moFock_spin.cwiseProduct(d2_Vector)).sum();
          }
        };
      }
    }
    // prefactor for results
    double prefactor = 1.0;
    if (SCFMode == Options::SCF_MODES::RESTRICTED)
      prefactor = 2.0;
    for_spin(qpEnergy, dsigma_de, z, vxc_energies, x_energies, correlation, orbEig) {
      correlation_spin = correlation_spin * prefactor;
      dsigma_de_spin = dsigma_de_spin * prefactor;
      z_spin = 1.0 / (1.0 - dsigma_de_spin.array());
      if (this->_settings.linearized) {
        qpEnergy_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) =
            orbEig_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
            z_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb)
                .cwiseProduct(x_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
                              correlation_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) -
                              vxc_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb));
      }
      else {
        qpEnergy_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) =
            orbEig_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
            x_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
            correlation_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) -
            vxc_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb);
      }
    };

    // Check for convergence and use convergence acceleration
    if (this->_settings.qpiterations > 0) {
      bool converged_qp = this->convergenCheck(qpEnergy, oldQP, diis_qpiterations, false, iQP);
      if (converged_qp)
        break;
    }
  } /* Loop QP iterations */

  // Write frequencies to file
  if (this->_settings.freq.size() == 3) {
    unsigned int spincounter = 0;
    for_spin(frequency_dependency) {
      std::string name = "FrequencyDependentSelfEnergy_";
      if (SCFMode == Options::SCF_MODES::RESTRICTED) {
        name += "res.txt";
      }
      else {
        std::string spin = (spincounter == 0) ? "unres_alpha.txt" : "unres_beta.txt";
        name += spin;
      }
      std::ofstream file(name);
      if (file.is_open()) {
        file << frequency_dependency_spin;
      }
      file.close();
      spincounter++;
    };
  }
}

template<Options::SCF_MODES SCFMode>
Eigen::VectorXd GW_Analytic<SCFMode>::calculatedVector(double freq, Eigen::VectorXd orbEigenValues,
                                                       double excitationEnergy, unsigned int nOcc, unsigned int nVirt) {
  Eigen::VectorXd d_Vector = Eigen::VectorXd::Zero(nOcc + nVirt);
  // Occupied
  d_Vector.segment(0, nOcc) = freq - orbEigenValues.segment(0, nOcc).array() + excitationEnergy;
  // Virtual
  d_Vector.segment(nOcc, nVirt) = freq - orbEigenValues.segment(nOcc, nVirt).array() - excitationEnergy;
  d_Vector = d_Vector.array() / (d_Vector.array().square() + pow(this->_settings.eta, 2));
  return d_Vector;
}

template class GW_Analytic<Options::SCF_MODES::RESTRICTED>;
template class GW_Analytic<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
