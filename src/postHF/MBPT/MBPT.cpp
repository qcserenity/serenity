/**
 * @file MBPT.cpp
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
#include "postHF/MBPT/MBPT.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "grid/GaussLegendre.h"
#include "io/FormattedOutputStream.h"
#include "math/LaplaceMinimaxWrapper.h"
#include "math/diis/DIIS.h"
#include "math/linearAlgebra/MatrixFunctions.h" //Matrix sqrt
#include "misc/HelperFunctions.h"               //Sparse prjections from sparse maps.
#include "parameters/Constants.h"
#include "postHF/LRSCF/Tools/RIIntegrals.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <iomanip>
namespace Serenity {

template<Options::SCF_MODES SCFMode>
MBPT<SCFMode>::MBPT(std::shared_ptr<LRSCFController<SCFMode>> lrscf, GWTaskSettings settings,
                    std::vector<std::shared_ptr<SystemController>> envSystemController,
                    std::shared_ptr<RIIntegrals<SCFMode>> riInts, int startOrb, int endOrb)
  : _lrscfController(lrscf),
    _settings(settings),
    _env(envSystemController),
    _riInts(riInts),
    _orbEig(_lrscfController->getEigenvalues()),
    _startOrb(startOrb),
    _endOrb(endOrb),
    _Jia(nullptr),
    _Jpq(nullptr),
    _Jia_transformed(nullptr),
    _Jpq_transformed(nullptr),
    _nOcc(_lrscfController->getNOccupied()),
    _nVirt(_lrscfController->getNVirtual()) {
  // Error and Warnings
  if (settings.nafThresh != 0 && _env.size() > 0) {
    throw SerenityError("NAF in combination with multiple subsystems not supported yet!");
  }

  // Only for methods relying on numerical integration
  if (settings.gwtype != Options::GWALGORITHM::ANALYTIC || settings.mbpttype == Options::MBPT::RPA) {
    if (!_riInts) {
      throw SerenityError("RI Integrals need to be initialized for dRPA/GW!");
    }
    if (_endOrb != _riInts->getPEnd() || _startOrb != _riInts->getPStart()) {
      throw SerenityError(
          "RI integrals use different orbital space for QP energies than specified for the GW calculation!");
    }
    // Initialize orbital energy differences
    _eia = this->calculateEia(_orbEig, _nOcc, _nVirt);
    // Initialize Grid
    if (_settings.integrationPoints > 0) {
      auto integrator = GaussLegendre(_settings.integrationPoints);
      Eigen::VectorXd nodes = integrator.getGridPoints();
      Eigen::VectorXd weights = integrator.getWeights();
      // Transform integral limits
      _nodes = ((1.0 + nodes.array()) / (1.0 - nodes.array())).eval();
      _weights = ((2.0 * weights.array()) / ((1.0 - nodes.array()).array().square())).eval();
    }
    _Jia = riInts->getJiaPtr();
    _nAux = riInts->getNTransformedAuxBasisFunctions();
    _Jpq = riInts->getJpqPtr();
    // needs to be at least initialized
    _Jia_transformed = std::make_shared<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>(1, 1);
    _Jpq_transformed = std::make_shared<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>(1, 1);
    // Calculate environment screening contributions and transformed jia
    if (_settings.environmentScreening) {
      _envResponse = this->environmentRespose();
      auto auxtrafo = riInts->getAuxTrafoPtr();
      auto metric = riInts->getAuxMetricPtr();
      auto metricsqrt = mSqrt_Sym(*metric);
      // V + V^{1/2}\PI^\env V^{1/2}
      Eigen::MatrixXd coulAuxbasis = (*metric) + metricsqrt * _envResponse * metricsqrt;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(coulAuxbasis);
      _eValues = eigenvalueSolver.eigenvalues();
      Eigen::MatrixXd eigenVectors = eigenvalueSolver.eigenvectors();

      auto& jia_transformed = *_Jia_transformed;
      auto& jia = *_Jia;

      for_spin(jia_transformed, jia) {
        jia_transformed_spin = jia_spin * (*auxtrafo) * eigenVectors;
      };
      // tranformation of jpq in case of GW
      if (_Jpq) {
        auto& jpq_transformed = *_Jpq_transformed;
        auto& jpq = *_Jpq;
        for_spin(jpq_transformed, jpq) {
          jpq_transformed_spin = jpq_spin * (*auxtrafo) * eigenVectors;
        };
      }
    }
  }
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd>
MBPT<SCFMode>::calculateEia(SpinPolarizedData<SCFMode, Eigen::VectorXd>& orbEigenValues,
                            SpinPolarizedData<SCFMode, unsigned int>& nOcc, SpinPolarizedData<SCFMode, unsigned int>& nVirt) {
  SpinPolarizedData<SCFMode, Eigen::VectorXd> e_ia;
  for_spin(orbEigenValues, nOcc, nVirt, e_ia) {
    e_ia_spin.resize(nOcc_spin * nVirt_spin);
    for (unsigned int ia = 0; ia < nOcc_spin * nVirt_spin; ++ia) {
      unsigned int i = floor(ia / nVirt_spin);
      unsigned int a = nOcc_spin + ia - i * nVirt_spin;
      e_ia_spin(ia) = orbEigenValues_spin(a) - orbEigenValues_spin(i);
    }
  };
  return e_ia;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd MBPT<SCFMode>::environmentRespose() {
  Timings::takeTime("MBPT -  W_env calculation");
  OutputControl::mOut << " Calculate Environment contributions ... \n" << std::endl;
  Eigen::MatrixXd envPQ = Eigen::MatrixXd::Zero(_nAux, _nAux);
  for (auto& env : _env) {
    auto riInts_J =
        RIIntegrals<SCFMode>(env, LIBINT_OPERATOR::coulomb, 0.0, true, 0, 0, _settings.nafThresh, _riInts->getGeo());
    auto orbEig_J = env->template getActiveOrbitalController<SCFMode>()->getEigenvalues();
    auto nOcc_J = env->template getNOccupiedOrbitals<SCFMode>();
    auto nVirt_J = env->template getNVirtualOrbitalsTruncated<SCFMode>();
    auto e_ia_J = this->calculateEia(orbEig_J, nOcc_J, nVirt_J);
    auto& jia_J = *(riInts_J.getJiaPtr());
    for_spin(e_ia_J, jia_J) {
      double prefactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
      Eigen::VectorXd chi_temp = prefactor * (-2.0 / e_ia_J_spin.array());
      envPQ += jia_J_spin.transpose() * chi_temp.asDiagonal() * jia_J_spin;
    };
  }
  OutputControl::mOut << " \n ... done\n" << std::endl;
  Timings::timeTaken("MBPT -  W_env calculation");
  return envPQ;
}

template<Options::SCF_MODES SCFMode>
Eigen::SparseMatrix<double> MBPT<SCFMode>::calculateTransformation(Eigen::MatrixXd& transformation,
                                                                   Eigen::MatrixXd& environmentResponse) {
  Eigen::SparseMatrix<double> projector;
  // Symmetrization
  environmentResponse = (-0.5 * (environmentResponse.transpose() + environmentResponse)).eval();
  // Custom Cholesky decomposition
  auto matCopy = environmentResponse; // TwoC is the matrix to decompose
  auto vec = environmentResponse;
  vec.setZero();
  Eigen::VectorXd diag = matCopy.diagonal();
  double thresh = 1e-10;
  unsigned int index = 0;
  unsigned int counter = 0;
  while (true) {
    if (diag.maxCoeff(&index) <= thresh)
      break;

    double factor = 1.0 / std::sqrt(diag[index]);
    vec.col(counter) = factor * matCopy.col(index);
    matCopy -= vec.col(counter) * vec.col(counter).transpose();
    matCopy.col(index).setZero();
    matCopy.row(index).setZero();
    diag = matCopy.diagonal();
    counter++;
  }
  if (diag.minCoeff() < -1e-10) {
    OutputControl::mOut << " ***** WARNING: GW-Cholesky Decomposition: Matrix not positive definite? ***** " << std::endl;
  }
  // Sparse Map
  Eigen::VectorXi neglible = Eigen::VectorXi::Constant(vec.cols(), 1);
  neglible.segment(0, counter).array() = 0;
  projector = constructProjectionMatrix(neglible);
  transformation = (vec * projector).eval();
  return projector;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::MatrixXd> MBPT<SCFMode>::calculateWnmComplex() {
  OutputControl::mOut << "\n Calculate W_nm(iw) ....... ";
  OutputControl::mOut.flush();
  Timings::takeTime("MBPT -      W_nm calculation");
  Eigen::MatrixXd unit = Eigen::MatrixXd::Identity(_nAux, _nAux);
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> W_nm;
  // needed integrals
  auto& jia = *(this->_Jia);
  auto& jpq = *(this->_Jpq);
  auto& jia_transformed = *(this->_Jia_transformed);
  auto& jpq_transformed = *(this->_Jpq_transformed);
  for_spin(W_nm, jpq, jpq_transformed) {
    W_nm_spin.resize(jpq_spin.rows(), _nodes.size());
    W_nm_spin.setZero();
    for (unsigned int iFreq = 0; iFreq < _nodes.size(); iFreq++) {
      double prefactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
      Eigen::MatrixXd dielectric = unit;
      for_spin(_eia, jia, jia_transformed) {
        Eigen::VectorXd chi_temp =
            prefactor * (-2.0 * _eia_spin.array()) / (std::pow(_nodes(iFreq), 2) + _eia_spin.array().square());
        if (_env.size() == 0) {
          dielectric.noalias() -= (jia_spin.transpose() * chi_temp.asDiagonal() * jia_spin).eval();
        }
        else {
          auto dielectricScreen = (jia_transformed_spin.transpose() * chi_temp.asDiagonal() * jia_transformed_spin).eval();
          dielectric.noalias() -= dielectricScreen * _eValues.asDiagonal();
        }
      };

      // other screening contribution
      if (_env.size() > 0) {
        Eigen::MatrixXd transformed;
        transformed = dielectric.householderQr().solve(jpq_transformed_spin.transpose());
        W_nm_spin.col(iFreq) = (jpq_transformed_spin * _eValues.asDiagonal() * transformed).diagonal() -
                               (jpq_spin * jpq_spin.transpose()).diagonal();
      }
      else {
        Eigen::LLT<Eigen::MatrixXd> llt;
        dielectric = llt.compute(dielectric).matrixL();
        if (llt.info() != Eigen::Success) {
          OutputControl::mOut << "Cholesky decomposition failed!" << std::endl;
        }
        Eigen::MatrixXd transformed;
        transformed = dielectric.triangularView<Eigen::Lower>().solve(jpq_spin.transpose());
        W_nm_spin.col(iFreq) =
            (transformed.transpose() * transformed).diagonal() - (jpq_spin * jpq_spin.transpose()).diagonal();
      }
    }
  };
  Timings::timeTaken("MBPT -      W_nm calculation");
  OutputControl::mOut << " done" << std::endl;
  return W_nm;
}
/* This Laplace-transform implementation is still experimental!
   A manual entry and more documentation will be added after further verification*/
template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::MatrixXd> MBPT<SCFMode>::calculateWnmComplexLT() {
  OutputControl::mOut << "\n Calculate W_nm(iw) ....... \n";
  OutputControl::mOut.flush();
  Timings::takeTime("MBPT -      W_nm calculation");
  Eigen::MatrixXd unit = Eigen::MatrixXd::Identity(_nAux, _nAux);
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> W_nm;
  // needed integrals
  auto& jia = *(this->_Jia);
  auto& jpq = *(this->_Jpq);
  auto& jia_transformed = *(this->_Jia_transformed);
  auto& jpq_transformed = *(this->_Jpq_transformed);
  double min = 1000;
  double max = 0.0;
  double freqMin = _nodes(1);
  for_spin(_orbEig, _nOcc) {
    double LUMO = _orbEig_spin(_nOcc_spin);
    double HOMO = _orbEig_spin(_nOcc_spin - 1);
    double minSpin = std::pow((LUMO - HOMO), 2) + std::pow(freqMin, 2);

    double HUMO = _orbEig_spin(_orbEig_spin.size() - 1);
    double LOMO = _orbEig_spin(0);
    double maxSpin = std::pow((HUMO - LOMO), 2);

    min = std::min(minSpin, min);
    max = std::max(maxSpin, max);
  };
  Eigen::VectorXd roots(0);
  Eigen::VectorXd weights(0);
  max = std::max(1e4, max);
  getMinimaxRoots(roots, weights, min, max, _settings.ltconv);

  auto intermediate = std::vector<Eigen::MatrixXd>(roots.size(), Eigen::MatrixXd::Zero(_nAux, _nAux));
  for_spin(_eia, jia, jia_transformed) {
    for (unsigned int m = 0; m < roots.size(); m++) {
      Eigen::VectorXd prefactor = _eia_spin;
      prefactor = -1.0 * prefactor.array().square().matrix() * roots(m);
      prefactor = prefactor.array().exp().matrix();
      prefactor = weights(m) * prefactor.cwiseProduct(_eia_spin);
      if (_env.size() == 0) {
        intermediate[m] += (jia_spin.transpose() * prefactor.asDiagonal() * jia_spin).eval();
      }
      else {
        intermediate[m] +=
            (jia_transformed_spin.transpose() * prefactor.asDiagonal() * jia_transformed_spin * _eValues.asDiagonal()).eval();
      }
    }
  };
  for_spin(W_nm, jpq, jpq_transformed) {
    W_nm_spin.resize(jpq_spin.rows(), _nodes.size());
    W_nm_spin.setZero();
    for (unsigned int iFreq = 0; iFreq < _nodes.size(); iFreq++) {
      double prefactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
      Eigen::MatrixXd dielectric = unit;
      for (unsigned int m = 0; m < roots.size(); m++) {
        dielectric -= -2.0 * prefactor * std::exp(-1.0 * std::pow(_nodes(iFreq), 2) * roots(m)) * intermediate[m];
      }
      if (_env.size() == 0) {
        Eigen::LLT<Eigen::MatrixXd> llt;
        dielectric = llt.compute(dielectric).matrixL();
        if (llt.info() != Eigen::Success) {
          OutputControl::mOut << "Cholesky decomposition failed!" << std::endl;
        }
        Eigen::MatrixXd transformed;
        transformed = dielectric.triangularView<Eigen::Lower>().solve(jpq_spin.transpose());
        W_nm_spin.col(iFreq) =
            (transformed.transpose() * transformed).diagonal() - (jpq_spin * jpq_spin.transpose()).diagonal();
      }
      else {
        Eigen::MatrixXd transformed;
        transformed = dielectric.householderQr().solve(jpq_transformed_spin.transpose());
        W_nm_spin.col(iFreq) = (jpq_transformed_spin * _eValues.asDiagonal() * transformed).diagonal() -
                               (jpq_spin * jpq_spin.transpose()).diagonal();
      }
    }
  };
  Timings::timeTaken("MBPT -      W_nm calculation");
  OutputControl::mOut << " done" << std::endl;
  return W_nm;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd MBPT<SCFMode>::calculatePiOmega(std::complex<double>& frequency, Eigen::MatrixXd& jia, Eigen::VectorXd& e_ia) {
  Eigen::VectorXcd prefactor(e_ia.size());
  prefactor = 1.0 / (frequency - e_ia.array());
  prefactor = prefactor.array() - 1.0 / (frequency + e_ia.array());
  Eigen::MatrixXd piPQ = jia.transpose() * prefactor.real().asDiagonal() * jia;
  if (SCFMode == Options::SCF_MODES::RESTRICTED)
    piPQ = 2.0 * piPQ;
  return piPQ;
}

template<Options::SCF_MODES SCFMode>
bool MBPT<SCFMode>::convergenCheck(SpinPolarizedData<SCFMode, Eigen::VectorXd>& new_qp,
                                   SpinPolarizedData<SCFMode, Eigen::VectorXd>& old_qp, std::shared_ptr<DIIS> diis,
                                   bool evGW, unsigned int iteration) {
  if (evGW) {
    OutputControl::mOut << "\n EV-GW Iteration:  " << iteration << std::endl;
  }
  else {
    OutputControl::mOut << "\n QP Iteration:  " << iteration << std::endl;
  }
  if (_settings.diis) {
    Timings::takeTime("MBPT -                  DIIS");
    for_spin(new_qp, old_qp) {
      new_qp_spin = ((1.0 - _settings.damping) * new_qp_spin + _settings.damping * old_qp_spin).eval();
      auto eigenspace = new_qp_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb);
      auto difference = eigenspace - old_qp_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb);
      diis->optimize(eigenspace, difference);
      new_qp_spin.segment(this->_startOrb, this->_endOrb - this->_startOrb) = eigenspace;
    };
    Timings::timeTaken("MBPT -                  DIIS");
  }
  this->printQPInfo(new_qp, old_qp);
  // check for convergence
  bool converged = true;
  unsigned spincounter = 0;
  for_spin(new_qp, old_qp, _nOcc) {
    double temp = std::abs(new_qp_spin(_nOcc_spin - 1) - old_qp_spin(_nOcc_spin - 1));
    temp += std::abs(new_qp_spin(_nOcc_spin) - old_qp_spin(_nOcc_spin));
    temp *= 0.5;
    if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
      if (spincounter == 0) {
        OutputControl::mOut << "\n (HOMO-LUMO) Alpha Gap deviation  : " << temp << " a.u." << std::endl;
      }
      else {
        OutputControl::mOut << " (HOMO-LUMO) Beta  Gap deviation  : " << temp << " a.u." << std::endl;
      }
    }
    else {
      OutputControl::mOut << "\n HOMO-LUMO Gap deviation  : " << temp << " a.u." << std::endl;
    }
    spincounter++;
    converged = (temp < _settings.ConvergenceThreshold) ? true : false;
  };
  if (converged) {
    OutputControl::mOut << "\n *** Convergence criteria reached: ";
    if (evGW) {
      OutputControl::mOut << "evGW Finished ***" << std::endl;
    }
    else {
      OutputControl::mOut << "QP Iteration Finished ***" << std::endl;
    }
  }
  return converged;
}

template<Options::SCF_MODES SCFMode>
void MBPT<SCFMode>::printQPInfo(SpinPolarizedData<SCFMode, Eigen::VectorXd>& new_qp,
                                SpinPolarizedData<SCFMode, Eigen::VectorXd>& old_qp) {
  // Write output for qp iterations
  printf("%3s %3s %15s %15s %15s \n", "", " # ", " QP old / eV ", " QP new / eV ", "   Delta E / eV ");
  printf("%3s %3s %15s %15s %15s \n", "", "---", "----------------", "--------------", "--------------");
  unsigned spincounter = 0;
  for_spin(new_qp, old_qp) {
    if (spincounter > 0) {
      std::cout << "    ----------------------------------------------------  " << std::endl;
    }
    for (int i = this->_startOrb; i < this->_endOrb; i++) {
      printf("%3s %3d %+15.5f %+15.5f %+15.5f \n", "", (i + 1), old_qp_spin(i) * HARTREE_TO_EV,
             new_qp_spin(i) * HARTREE_TO_EV, (new_qp_spin(i) - old_qp_spin(i)) * HARTREE_TO_EV);
    }
    spincounter++;
  };
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, double> MBPT<SCFMode>::calculateFermiLevel() {
  SpinPolarizedData<SCFMode, double> fermi_Level;
  unsigned spincounter = 0;
  for_spin(_orbEig, fermi_Level, _nOcc) {
    fermi_Level_spin = (_orbEig_spin(_nOcc_spin) + _orbEig_spin(_nOcc_spin - 1)) / 2.0;
    if (SCFMode == Options::SCF_MODES::UNRESTRICTED) {
      if (spincounter == 0) {
        OutputControl::mOut << " Fermi - Level Alpha           : " << fermi_Level_spin * HARTREE_TO_EV << " eV" << std::endl;
      }
      else {
        OutputControl::mOut << " Fermi - Level Beta            : " << fermi_Level_spin * HARTREE_TO_EV << " eV" << std::endl;
      }
    }
    else {
      OutputControl::mOut << " Fermi - Level            : " << fermi_Level_spin * HARTREE_TO_EV << " eV" << std::endl;
    }
    spincounter++;
  };
  return fermi_Level;
}

template<Options::SCF_MODES SCFMode>
void MBPT<SCFMode>::shiftByGap(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEner,
                               SpinPolarizedData<SCFMode, Eigen::VectorXd>& startOrbEner) {
  OutputControl::nOut
      << " ** All not included orbitals are shifted by the GAP of the highest and lowest included orbital! ** " << std::endl;
  // change orbital energy for the states, which are not included only using exchange
  for_spin(qpEner, startOrbEner, _nOcc, _nVirt) {
    double gap_occ = qpEner_spin(_startOrb) - startOrbEner_spin(_startOrb);
    double gap_lumo = qpEner_spin(_endOrb - 1) - startOrbEner_spin(_endOrb - 1);
    for (int iState = 0; iState < int(_nOcc_spin + _nVirt_spin); iState++) {
      if (iState < _startOrb || iState > _endOrb - 1) {
        if ((unsigned int)iState < _nOcc_spin) {
          qpEner_spin(iState) = startOrbEner_spin(iState) + gap_occ;
        }
        else {
          qpEner_spin(iState) = startOrbEner_spin(iState) + gap_lumo;
        }
      }
    }
  };
}

template<Options::SCF_MODES SCFMode>
double MBPT<SCFMode>::calculateRPACorrelationEnergy() {
  // Integration
  Timings::takeTime("MBPT -    RPA Corr. Ener.");
  double dRPA_correlation = 0.0;
  const Eigen::MatrixXd unit = Eigen::MatrixXd::Identity(_nAux, _nAux);
  auto& jia = *_Jia;
#pragma omp parallel for schedule(dynamic) reduction(+ : dRPA_correlation)
  for (unsigned int i = 0; i < _settings.integrationPoints; i++) {
    Eigen::MatrixXd pi_PQ = Eigen::MatrixXd::Zero(_nAux, _nAux);
    double prefactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
    for_spin(jia, _eia) {
      Eigen::VectorXd chi_temp =
          prefactor * (-2.0 * _eia_spin.array()) / (std::pow(_nodes(i), 2) + _eia_spin.array().square());
      pi_PQ += jia_spin.transpose() * chi_temp.asDiagonal() * jia_spin;
    };
    // If additional environmental screening is taken into account
    if (_env.size() > 0) {
      throw SerenityError("Environmental-Screened RPA not implemented!");
    }
    dRPA_correlation += (1.0 / (2.0 * PI)) * _weights(i) * (log((unit - pi_PQ.real()).determinant()) + pi_PQ.real().trace());
  }
  Timings::timeTaken("MBPT -    RPA Corr. Ener.");
  return dRPA_correlation;
}

template class MBPT<Options::SCF_MODES::RESTRICTED>;
template class MBPT<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
