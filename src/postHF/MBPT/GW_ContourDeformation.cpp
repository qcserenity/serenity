/**
 * @file GW_ContourDeformation.cpp
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
#include "postHF/MBPT/GW_ContourDeformation.h"
/* Include Serenity Internal Headers */
#include "math/diis/DIIS.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
GW_ContourDeformation<SCFMode>::GW_ContourDeformation(std::shared_ptr<LRSCFController<SCFMode>> lrscf, GWTaskSettings settings,
                                                      std::vector<std::shared_ptr<SystemController>> envSystemController,
                                                      std::shared_ptr<RIIntegrals<SCFMode>> riInts, int startOrb, int endOrb)
  : MBPT<SCFMode>(lrscf, settings, envSystemController, riInts, startOrb, endOrb) {
}

template<Options::SCF_MODES SCFMode>
void GW_ContourDeformation<SCFMode>::calculateGWOrbitalenergies(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                                                SpinPolarizedData<SCFMode, Eigen::VectorXd>& vxc_energies,
                                                                SpinPolarizedData<SCFMode, Eigen::VectorXd>& x_energies,
                                                                SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation) {
  if (SCFMode == Options::SCF_MODES::UNRESTRICTED && this->_env.size() > 0) {
    throw SerenityError("Unrestricted GW-CD with environmental screening not supported, yet!");
  }
  // Check if calculation is evGW
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
    // fermi level
    auto fermiLevel = this->calculateFermiLevel();
    SpinPolarizedData<SCFMode, Eigen::MatrixXd> wnm;
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
      // correlation for imaginary axis integration
      Timings::takeTime("MBPT -     Numerical Integr.");
      auto imagCorr = this->calculateImagCorr(wnm, qpEnergy);
      Timings::timeTaken("MBPT -     Numerical Integr.");
      // real correlation contribution
      Timings::takeTime("MBPT -     Countour Residues");
      auto realCorr = this->calculateContourResidues(qpEnergy, fermiLevel);
      Timings::timeTaken("MBPT -     Countour Residues");
      // Set new correlation energy, qpEnergies etc.
      for_spin(realCorr, imagCorr, correlation, qpEnergy, x_energies, vxc_energies, initialOrbEig) {
        correlation_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) = realCorr_spin + imagCorr_spin;
        qpEnergy_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) =
            initialOrbEig_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
            x_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) +
            correlation_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) -
            vxc_energies_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb);
      };
      // Check for convergence and use convergence acceleration
      if (this->_settings.qpiterations > 0) {
        bool converged_qp = this->convergenCheck(qpEnergy, oldQP, diis_qpiterations, false, iQP);
        if (converged_qp)
          break;
      }
    } /* Loop QP iterations */
    // Shift not included orbital energies by gap
    if (this->_settings.gap) {
      this->shiftByGap(qpEnergy, initialOrbEig);
    }
    // Check for convergence and use convergence acceleration
    if (this->_settings.evGW) {
      bool converged_evGW = this->convergenCheck(qpEnergy, oldEvGWQP, diis_evGW, true, iGW);
      if (converged_evGW)
        break;
    }
  } /* Loop evGW iterations */
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd>
GW_ContourDeformation<SCFMode>::calculateContourResidues(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                                                         SpinPolarizedData<SCFMode, double>& fermiLevel) {
  // contour residue correlation contribution along real axis
  SpinPolarizedData<SCFMode, Eigen::VectorXd> realPart;
  for_spin(realPart) {
    realPart_spin = Eigen::VectorXd::Zero(this->_endOrb - this->_startOrb);
  };
  // check if environment screening should be included
  bool environment = false;
  if (this->_env.size() > 0) {
    environment = true;
  }
  Eigen::MatrixXd unit = Eigen::MatrixXd::Identity(this->_nAux, this->_nAux);
  for (int iState = this->_startOrb; iState < this->_endOrb; iState++) {
    // needed integrals
    auto& jia = *(this->_Jia);
    auto& jpq = *(this->_Jpq);
    auto& jia_transformed = *(this->_Jia_transformed);
    auto& jpq_transformed = *(this->_Jpq_transformed);
    // other variables
    auto& eia = this->_eia;
    auto& orbEig = this->_orbEig;
    auto& nOcc = this->_nOcc;
    auto& nVirt = this->_nVirt;
    for_spin(realPart, orbEig, jpq, qpEnergy, fermiLevel, jpq_transformed, nOcc, nVirt) {
      unsigned states = nOcc_spin + nVirt_spin;
      for (unsigned int iRes = 0; iRes < states; iRes++) {
        // residue contribution
        double f = 0.0;
        if (double_smaller(fermiLevel_spin, orbEig_spin(iRes), 1e-6) &&
            double_smaller(orbEig_spin(iRes), qpEnergy_spin(iState), 1e-6)) {
          f = 1.0;
        }
        else if (double_smaller(qpEnergy_spin(iState), orbEig_spin(iRes), 1e-6) &&
                 double_smaller(orbEig_spin(iRes), fermiLevel_spin, 1e-6)) {
          f = -1.0;
        }
        else if (double_smaller(fermiLevel_spin, orbEig_spin(iRes), 1e-6) &&
                 double_equals(orbEig_spin(iRes), qpEnergy_spin(iState), 1e-6)) {
          f = 0.5;
        }
        else if (double_equals(orbEig_spin(iRes), qpEnergy_spin(iState), 1e-6) &&
                 double_smaller(orbEig_spin(iRes), fermiLevel_spin, 1e-6)) {
          f = -0.5;
        }
        if (f == 0.0)
          continue;
        // residue frequency
        auto freq = std::complex<double>(std::abs(orbEig_spin(iRes) - qpEnergy_spin(iState)), this->_settings.eta);
        Timings::takeTime("MBPT -CD Dielectric Residues");
        Eigen::MatrixXd dielectric = unit;
        for_spin(eia, jia, jia_transformed) {
          if (!environment) {
            dielectric.noalias() -= this->calculatePiOmega(freq, jia_spin, eia_spin);
          }
          else {
            auto dielectricScreen = this->calculatePiOmega(freq, jia_transformed_spin, eia_spin);
            dielectric.noalias() -= dielectricScreen * (this->_eValues.asDiagonal());
          }
        };
        Timings::timeTaken("MBPT -CD Dielectric Residues");
        double value = 0.0;
        // other environmental screening
        if (environment) {
          Eigen::VectorXd temp = dielectric.householderQr().solve(
              jpq_transformed_spin.row((iState - this->_startOrb) * states + iRes).transpose());
          value += jpq_transformed_spin.row((iState - this->_startOrb) * states + iRes) * (this->_eValues.asDiagonal()) * temp;
          value -= jpq_spin.row((iState - this->_startOrb) * states + iRes) *
                   jpq_spin.row((iState - this->_startOrb) * states + iRes).transpose();
        }
        else {
          Timings::takeTime("MBPT -   CD LU Decomposition");
          Eigen::VectorXd jpqVec = jpq_spin.row((iState - this->_startOrb) * states + iRes);
          Eigen::VectorXd jpqtransformed = dielectric.lu().solve(jpqVec);
          Timings::timeTaken("MBPT -   CD LU Decomposition");
          value += jpqtransformed.transpose() * jpq_spin.row((iState - this->_startOrb) * states + iRes).transpose();
          value -= jpq_spin.row((iState - this->_startOrb) * states + iRes) *
                   jpq_spin.row((iState - this->_startOrb) * states + iRes).transpose();
        }
        realPart_spin(iState - this->_startOrb) += f * value;
      }
    };
  }
  return realPart;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd>
GW_ContourDeformation<SCFMode>::calculateImagCorr(SpinPolarizedData<SCFMode, Eigen::MatrixXd>& wnm,
                                                  SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy) {
  SpinPolarizedData<SCFMode, Eigen::VectorXd> imagPart;
  for_spin(imagPart) {
    imagPart_spin = Eigen::VectorXd::Zero(this->_endOrb - this->_startOrb);
  };
  auto& nOcc = this->_nOcc;
  auto& nVirt = this->_nVirt;
  auto& orbEig = this->_orbEig;
  for (int iState = this->_startOrb; iState < this->_endOrb; iState++) {
    for_spin(imagPart, qpEnergy, orbEig, wnm, nOcc, nVirt) {
      int states = nOcc_spin + nVirt_spin;
      for (int mState = 0; mState < states; mState++) {
        for (unsigned int iFreq = 0; iFreq < this->_nodes.size(); iFreq++) {
          double kim = (1.0 / PI);
          kim *= (qpEnergy_spin(iState) - orbEig_spin(mState)) /
                 (std::pow(qpEnergy_spin(iState) - orbEig_spin(mState), 2) + std::pow(this->_nodes(iFreq), 2));
          kim *= this->_weights(iFreq);
          imagPart_spin(iState - this->_startOrb) -= kim * wnm_spin((iState - this->_startOrb) * states + mState, iFreq);
        }
      }
    };
  }
  return imagPart;
}

template class GW_ContourDeformation<Options::SCF_MODES::RESTRICTED>;
template class GW_ContourDeformation<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity