/**
 * @file RIExchangeSigmavector.cpp
 *
 * @date Oct 07, 2021
 * @author Niklas Niemeyer
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
#include "postHF/LRSCF/Sigmavectors/RI/RIExchangeSigmavector.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/BasisFunctionMapper.h" //Construct joined fitting basis.
#include "dft/Functional.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/RI_J_IntegralControllerFactory.h"
#include "integrals/looper/TwoElecThreeCenterCalculator.h"
#include "io/FormattedOutputStream.h"
#include "math/linearAlgebra/MatrixFunctions.h"
#include "misc/Timing.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
RIExchangeSigmavector<SCFMode>::RIExchangeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                      std::vector<Eigen::MatrixXd> b, const std::vector<int> pm,
                                                      bool densFitK, bool densFitLRK)
  : ExchangeSigmavector<SCFMode>(lrscf, b, pm, densFitK, densFitLRK) {
}

template<Options::SCF_MODES SCFMode>
RIExchangeSigmavector<SCFMode>::RIExchangeSigmavector(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf,
                                                      const std::vector<int> pm, bool densFitK, bool densFitLRK)
  : ExchangeSigmavector<SCFMode>(lrscf, pm, densFitK, densFitLRK) {
  throw SerenityError("AO representation of the RIExchangeSigmavector is not implemented yet!");
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>
RIExchangeSigmavector<Options::SCF_MODES::RESTRICTED>::calcF(
    unsigned I, unsigned J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>> densityMatrices) {
  bool rpaScreen = this->_lrscf[I]->getLRSCFSettings().rpaScreening;
  this->setParameters(I, J, rpaScreen);
  if (rpaScreen && I != J) {
    throw SerenityError("Coupled BSE not implemented!");
  }

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::RESTRICTED>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Determine prescreening threshold.
  double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                          this->_lrscf[J]->getSysSettings().basis.integralThreshold);

  if (prescreeningThreshold == 0) {
    double prescreenTresholdI = this->_lrscf[I]->getBasisController()->getPrescreeningThreshold();
    double prescreenTresholdJ = this->_lrscf[J]->getBasisController()->getPrescreeningThreshold();
    prescreeningThreshold = std::min(prescreenTresholdI, prescreenTresholdJ);
  }

  double maxDens = 0.0;
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      maxDens = std::max(maxDens, (*densityMatrices)[iSet][iGuess].array().abs().maxCoeff());
    }
  }

  // System I.
  auto basisControllerI = this->_lrscf[I]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_I = basisControllerI->getNBasisFunctions();
  unsigned no_I = this->_lrscf[I]->getNOccupied();
  auto& C_I = this->_lrscf[I]->getCoefficients();
  const Eigen::MatrixXd& S_I = this->_lrscf[I]->getSys()->getOneElectronIntegralController()->getOverlapIntegrals();

  // System J.
  auto basisControllerJ = this->_lrscf[J]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_J = basisControllerJ->getNBasisFunctions();
  unsigned no_J = this->_lrscf[J]->getNOccupied();
  unsigned nv_J = this->_lrscf[J]->getNVirtual();
  auto& C_J = this->_lrscf[J]->getCoefficients();

  // Get half-transformed density matrix for each guess vector.
  std::vector<std::vector<Eigen::MatrixXd>> C_vj(_nSet);
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    C_vj[iSet] = std::vector<Eigen::MatrixXd>(_nGuess, Eigen::MatrixXd(nb_J, no_J));
    for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
      C_vj[iSet][iGuess] =
          C_J.middleCols(no_J, nv_J) * Eigen::Map<Eigen::MatrixXd>(_b[iSet][J].col(iGuess).data(), nv_J, no_J);
    }
  }

  // Prepare half-transformed Fock matrix for each guess vector.
  std::vector<std::vector<std::vector<Eigen::MatrixXd>>> F_mui(this->_nThreads);
  for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
    F_mui[iThread] = std::vector<std::vector<Eigen::MatrixXd>>(_nSet);
    for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
      F_mui[iThread][iSet] = std::vector<Eigen::MatrixXd>(_nGuess, Eigen::MatrixXd::Zero(nb_I, no_I));
    }
  }

  // Simple switch for TDA <-> TDDFT.
  bool TDDFT = false;
  for (int p : _pm) {
    if (p != 0) {
      TDDFT = true;
    }
  }

  // Conventional exchange.
  if (_densFitK && (_hfExchangeRatio > 0.0 || rpaScreen)) {
    OutputControl::dOut << " **** Performing RI-K ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector:   RI-K");

    // Determine common basis purpose.
    auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
    if (this->_lrscf[I]->getLRSCFSettings().densFitK == Options::DENS_FITS::ACD) {
      basispurpose = Options::BASIS_PURPOSES::ATOMIC_CHOLESKY;
    }
    else if (this->_lrscf[I]->getLRSCFSettings().densFitK == Options::DENS_FITS::ACCD) {
      basispurpose = Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY;
    }

    // Initialize aux basis controller.
    auto auxBasisControllerI = this->_lrscf[I]->getSys()->getBasisController(basispurpose);
    auto auxBasisControllerJ = this->_lrscf[J]->getSys()->getBasisController(basispurpose);

    // Bases for looper.
    std::shared_ptr<BasisController> superBas = nullptr;
    std::shared_ptr<BasisController> superAuxBas = nullptr;

    if (I == J) {
      superBas = basisControllerI;
      superAuxBas = auxBasisControllerI;
    }
    else {
      BasisFunctionMapper basisMapper(basisControllerI);
      BasisFunctionMapper auxBasisMapper(auxBasisControllerI);
      superBas = basisMapper.getCombinedBasis(basisControllerJ);
      superAuxBas = auxBasisMapper.getCombinedBasis(auxBasisControllerJ);
    }
    unsigned nx = superAuxBas->getNBasisFunctions();
    unsigned nxs = superAuxBas->getReducedNBasisFunctions();
    auto& auxShells = superAuxBas->getBasis();

    auto riints = RI_J_IntegralControllerFactory::getInstance().produce(superBas, superAuxBas);
    auto invM = _lrscf[I]->getInverseMetric();
    if (rpaScreen == true && invM == nullptr) {
      Eigen::MatrixXd invSqrt = riints->getInverseMSqrt();
      auto auxScreen = _lrscf[I]->getScreeningAuxMatrix();
      invSqrt = (invSqrt * (*auxScreen) * invSqrt).eval();
      invM = std::make_shared<Eigen::MatrixXd>(invSqrt);
      _lrscf[I]->setInverseMetric(invM);
    }
    else if (rpaScreen == false && invM == nullptr) {
      invM = std::make_shared<Eigen::MatrixXd>(riints->getInverseM());
      _lrscf[I]->setInverseMetric(invM);
    }

    TwoElecThreeCenterCalculator intCalculator(LIBINT_OPERATOR::coulomb, 0, basisControllerJ, basisControllerI,
                                               superAuxBas, prescreeningThreshold, maxDens);

    // Calculate (Q|ji) once.
    Eigen::MatrixXd Qji = Eigen::MatrixXd(no_J * no_I, nx);

    // Calculate (Q|\tilde ji) for each set and guess vector.
    std::vector<std::vector<Eigen::MatrixXd>> Qji_t(_nSet);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        Qji_t[iSet] = std::vector<Eigen::MatrixXd>(_nGuess, Eigen::MatrixXd(no_J * no_I, nx));
      }
    }

// First contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> kappaJ_lambdaI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::MatrixXd kappaJ_iI = kappaJ_lambdaI * C_I.middleCols(0, no_I);
        Eigen::Map<Eigen::MatrixXd> ji(Qji.col(P_all).data(), no_J, no_I);
        ji = C_J.middleCols(0, no_J).transpose() * kappaJ_iI;
        if (TDDFT) {
          for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
            for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
              Eigen::Map<Eigen::MatrixXd> ji_t(Qji_t[iSet][iGuess].col(P_all).data(), no_J, no_I);
              ji_t = C_vj[iSet][iGuess].middleCols(0, no_J).transpose() * kappaJ_iI;
            }
          }
        }
      }
    }

    // Contract with inverse metric matrix.
    Qji *= (*invM);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
          Qji_t[iSet][iGuess] *= (*invM);
        }
      }
    }

// Second contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> nuJ_muI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::Map<Eigen::MatrixXd> ji(Qji.col(P_all).data(), no_J, no_I);
        for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
          for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
            Eigen::Ref<Eigen::MatrixXd> F = F_mui[iThread][iSet][iGuess];
            if (!TDDFT) {
              // Contribution of A matrix.
              F.noalias() += nuJ_muI.transpose() * C_vj[iSet][iGuess] * ji;
            }
            else {
              // Contribution of (A pm B) matrix.
              Eigen::Map<Eigen::MatrixXd> ji_t(Qji_t[iSet][iGuess].col(P_all).data(), no_J, no_I);
              F.noalias() += nuJ_muI.transpose() * (C_vj[iSet][iGuess] * ji + _pm[iSet] * C_J.middleCols(0, no_J) * ji_t);
            }
          }
        }
      }
    }
    Timings::timeTaken("LRSCF -   Sigmavector:   RI-K");

    // Add HF exchange to Fock matrix.
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    Eigen::MatrixXd temp = S_I * C_I.middleCols(0, no_I);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Eigen::Ref<Eigen::MatrixXd> Fx = F_mui[0][iSet][iGuess];
        for (unsigned iThread = 1; iThread < this->_nThreads; ++iThread) {
          Fx.noalias() += F_mui[iThread][iSet][iGuess];
          F_mui[iThread][iSet][iGuess].setZero();
        }
        (*fock)[iSet][iGuess] += _hfExchangeRatio * temp * Fx.transpose();
        F_mui[0][iSet][iGuess].setZero();
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  } /* Exchange */

  // Long-range exchange.
  if (_densFitLRK && _lrExchangeRatio > 0.0) {
    OutputControl::dOut << " **** Performing RI-ErfK ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector: RI-LRK");

    // Determine common basis purpose.
    auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
    if (this->_lrscf[I]->getLRSCFSettings().densFitLRK == Options::DENS_FITS::ACD) {
      basispurpose = Options::BASIS_PURPOSES::ERF_ATOMIC_CHOLESKY;
    }
    else if (this->_lrscf[I]->getLRSCFSettings().densFitLRK == Options::DENS_FITS::ACCD) {
      basispurpose = Options::BASIS_PURPOSES::ERF_ATOMIC_COMPACT_CHOLESKY;
    }
    // Initialize aux basis controller
    auto auxBasisControllerI = this->_lrscf[I]->getBasisController(basispurpose);
    auto auxBasisControllerJ = this->_lrscf[J]->getBasisController(basispurpose);

    // Bases for looper.
    std::shared_ptr<BasisController> superBas = nullptr;
    std::shared_ptr<BasisController> superAuxBas = nullptr;

    if (I == J) {
      superBas = basisControllerI;
      superAuxBas = auxBasisControllerI;
    }
    else {
      BasisFunctionMapper basisMapper(basisControllerI);
      BasisFunctionMapper auxBasisMapper(auxBasisControllerI);
      superBas = basisMapper.getCombinedBasis(basisControllerJ);
      superAuxBas = auxBasisMapper.getCombinedBasis(auxBasisControllerJ);
    }
    unsigned nx = superAuxBas->getNBasisFunctions();
    unsigned nxs = superAuxBas->getReducedNBasisFunctions();
    auto& auxShells = superAuxBas->getBasis();

    RI_J_IntegralController riints(superBas, superAuxBas, nullptr, LIBINT_OPERATOR::erf_coulomb, _mu);

    auto invM = _lrscf[I]->getInverseErfMetric();
    if (invM == nullptr) {
      invM = std::make_shared<Eigen::MatrixXd>(riints.getInverseM());
      _lrscf[I]->setInverseErfMetric(invM);
    }

    TwoElecThreeCenterCalculator intCalculator(LIBINT_OPERATOR::erf_coulomb, _mu, basisControllerJ, basisControllerI,
                                               superAuxBas, prescreeningThreshold, maxDens);
    // Calculate (Q|ji) once.
    Eigen::MatrixXd Qji = Eigen::MatrixXd(no_J * no_I, nx);

    // Calculate (Q|\tilde ji) for each set and guess vector.
    std::vector<std::vector<Eigen::MatrixXd>> Qji_t(_nSet);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        Qji_t[iSet] = std::vector<Eigen::MatrixXd>(_nGuess, Eigen::MatrixXd(no_J * no_I, nx));
      }
    }

// First contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> kappaJ_lambdaI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::MatrixXd kappaJ_iI = kappaJ_lambdaI * C_I.middleCols(0, no_I);
        Eigen::Map<Eigen::MatrixXd> ji(Qji.col(P_all).data(), no_J, no_I);
        ji = C_J.middleCols(0, no_J).transpose() * kappaJ_iI;
        if (TDDFT) {
          for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
            for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
              Eigen::Map<Eigen::MatrixXd> ji_t(Qji_t[iSet][iGuess].col(P_all).data(), no_J, no_I);
              ji_t = C_vj[iSet][iGuess].middleCols(0, no_J).transpose() * kappaJ_iI;
            }
          }
        }
      }
    }

    // Contract with inverse metric matrix.
    Qji *= (*invM);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
          Qji_t[iSet][iGuess] *= (*invM);
        }
      }
    }

// Second contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> nuJ_muI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::Map<Eigen::MatrixXd> ji(Qji.col(P_all).data(), no_J, no_I);
        for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
          for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
            Eigen::Ref<Eigen::MatrixXd> F = F_mui[iThread][iSet][iGuess];
            if (!TDDFT) {
              // Contribution of A matrix.
              F += nuJ_muI.transpose() * C_vj[iSet][iGuess] * ji;
            }
            else {
              // Contribution of (A pm B) matrix.
              Eigen::Map<Eigen::MatrixXd> ji_t(Qji_t[iSet][iGuess].col(P_all).data(), no_J, no_I);
              F += nuJ_muI.transpose() * (C_vj[iSet][iGuess] * ji + _pm[iSet] * C_J.middleCols(0, no_J) * ji_t);
            }
          }
        }
      }
    }
    Timings::timeTaken("LRSCF -   Sigmavector: RI-LRK");

    // Add LR exchange to Fock matrix.
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    Eigen::MatrixXd temp = S_I * C_I.middleCols(0, no_I);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Eigen::Ref<Eigen::MatrixXd> Fx = F_mui[0][iSet][iGuess];
        for (unsigned iThread = 1; iThread < this->_nThreads; ++iThread) {
          Fx.noalias() += F_mui[iThread][iSet][iGuess];
        }
        (*fock)[iSet][iGuess] += _lrExchangeRatio * temp * Fx.transpose();
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  } /* LR-Exchange */

  return fock;
}

template<>
std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>
RIExchangeSigmavector<Options::SCF_MODES::UNRESTRICTED>::calcF(
    unsigned int I, unsigned int J,
    std::unique_ptr<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>> densityMatrices) {
  bool rpaScreen = this->_lrscf[I]->getLRSCFSettings().rpaScreening;
  this->setParameters(I, J, rpaScreen);
  if (rpaScreen && I != J) {
    throw SerenityError("Coupled BSE not implemented!");
  }

  // Set dimensions for Fock like matrices.
  auto fock = std::make_unique<std::vector<std::vector<MatrixInBasis<Options::SCF_MODES::UNRESTRICTED>>>>(this->_nSet);
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      (*fock)[iSet].emplace_back(this->_lrscf[I]->getBasisController());
    }
  }

  // Determine prescreening threshold.
  double prescreeningThreshold = std::min(this->_lrscf[I]->getSysSettings().basis.integralThreshold,
                                          this->_lrscf[J]->getSysSettings().basis.integralThreshold);

  if (prescreeningThreshold == 0) {
    double prescreenTresholdI = this->_lrscf[I]->getBasisController()->getPrescreeningThreshold();
    double prescreenTresholdJ = this->_lrscf[J]->getBasisController()->getPrescreeningThreshold();
    prescreeningThreshold = std::min(prescreenTresholdI, prescreenTresholdJ);
  }

  double maxDens = 0.0;
  for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
    for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
      maxDens = std::max(maxDens, (*densityMatrices)[iSet][iGuess].total().array().abs().maxCoeff());
    }
  }

  // System I.
  auto basisControllerI = this->_lrscf[I]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_I = basisControllerI->getNBasisFunctions();
  unsigned no_I_a = this->_lrscf[I]->getNOccupied().alpha;
  unsigned no_I_b = this->_lrscf[I]->getNOccupied().beta;
  auto& C_I = this->_lrscf[I]->getCoefficients();

  const Eigen::MatrixXd& S_I = this->_lrscf[I]->getSys()->getOneElectronIntegralController()->getOverlapIntegrals();

  // System J.
  auto basisControllerJ = this->_lrscf[J]->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  unsigned nb_J = basisControllerJ->getNBasisFunctions();
  unsigned no_J_a = this->_lrscf[J]->getNOccupied().alpha;
  unsigned nv_J_a = this->_lrscf[J]->getNVirtual().alpha;
  unsigned no_J_b = this->_lrscf[J]->getNOccupied().beta;
  unsigned nv_J_b = this->_lrscf[J]->getNVirtual().beta;
  auto& C_J = this->_lrscf[J]->getCoefficients();

  // Get half-transformed density matrix for each guess vector.
  std::vector<std::vector<std::vector<Eigen::MatrixXd>>> C_vj(_nSet);
  for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
    C_vj[iSet] = std::vector<std::vector<Eigen::MatrixXd>>(_nGuess);
    for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
      C_vj[iSet][iGuess] = std::vector<Eigen::MatrixXd>(2);
      C_vj[iSet][iGuess][0] = C_J.alpha.middleCols(no_J_a, nv_J_a) *
                              Eigen::Map<Eigen::MatrixXd>(_b[iSet][J].col(iGuess).data(), nv_J_a, no_J_a);
      C_vj[iSet][iGuess][1] =
          C_J.beta.middleCols(no_J_b, nv_J_b) *
          Eigen::Map<Eigen::MatrixXd>(_b[iSet][J].col(iGuess).data() + (long)nv_J_a * no_J_a, nv_J_b, no_J_b);
    }
  }

  // Prepare half-transformed Fock matrix for each guess vector.
  std::vector<std::vector<std::vector<std::vector<Eigen::MatrixXd>>>> F_mui(this->_nThreads);
  for (unsigned iThread = 0; iThread < this->_nThreads; ++iThread) {
    F_mui[iThread] = std::vector<std::vector<std::vector<Eigen::MatrixXd>>>(_nSet);
    for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
      F_mui[iThread][iSet] = std::vector<std::vector<Eigen::MatrixXd>>(_nGuess);
      for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
        F_mui[iThread][iSet][iGuess] = std::vector<Eigen::MatrixXd>(2);
        F_mui[iThread][iSet][iGuess][0] = Eigen::MatrixXd::Zero(nb_I, no_I_a);
        F_mui[iThread][iSet][iGuess][1] = Eigen::MatrixXd::Zero(nb_I, no_I_b);
      }
    }
  }

  // Simple switch for TDA <-> TDDFT.
  bool TDDFT = false;
  for (int p : _pm) {
    if (p != 0) {
      TDDFT = true;
    }
  }

  // Conventional exchange.
  if (_densFitK && (_hfExchangeRatio > 0.0 || rpaScreen)) {
    OutputControl::dOut << " **** Performing RI-K ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector:   RI-K");

    // Determine common basis purpose.
    auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
    if (this->_lrscf[I]->getLRSCFSettings().densFitK == Options::DENS_FITS::ACD) {
      basispurpose = Options::BASIS_PURPOSES::ATOMIC_CHOLESKY;
    }
    else if (this->_lrscf[I]->getLRSCFSettings().densFitK == Options::DENS_FITS::ACCD) {
      basispurpose = Options::BASIS_PURPOSES::ATOMIC_COMPACT_CHOLESKY;
    }

    // Initialize aux basis controller.
    auto auxBasisControllerI = this->_lrscf[I]->getBasisController(basispurpose);
    auto auxBasisControllerJ = this->_lrscf[J]->getBasisController(basispurpose);

    // Bases for looper.
    std::shared_ptr<BasisController> superBas = nullptr;
    std::shared_ptr<BasisController> superAuxBas = nullptr;

    if (I == J) {
      superBas = basisControllerI;
      superAuxBas = auxBasisControllerI;
    }
    else {
      BasisFunctionMapper basisMapper(basisControllerI);
      BasisFunctionMapper auxBasisMapper(auxBasisControllerI);
      superBas = basisMapper.getCombinedBasis(basisControllerJ);
      superAuxBas = auxBasisMapper.getCombinedBasis(auxBasisControllerJ);
    }
    unsigned nx = superAuxBas->getNBasisFunctions();
    unsigned nxs = superAuxBas->getReducedNBasisFunctions();
    auto& auxShells = superAuxBas->getBasis();

    auto riints = RI_J_IntegralControllerFactory::getInstance().produce(superBas, superAuxBas);
    auto invM = _lrscf[I]->getInverseMetric();
    if (rpaScreen == true && invM == nullptr) {
      Eigen::MatrixXd invSqrt = riints->getInverseMSqrt();
      auto auxScreen = _lrscf[I]->getScreeningAuxMatrix();
      invSqrt = (invSqrt * (*auxScreen) * invSqrt).eval();
      invM = std::make_shared<Eigen::MatrixXd>(invSqrt);
      _lrscf[I]->setInverseMetric(invM);
    }
    else if (rpaScreen == false && invM == nullptr) {
      invM = std::make_shared<Eigen::MatrixXd>(riints->getInverseM());
      _lrscf[I]->setInverseMetric(invM);
    }

    TwoElecThreeCenterCalculator intCalculator(LIBINT_OPERATOR::coulomb, 0, basisControllerJ, basisControllerI,
                                               superAuxBas, prescreeningThreshold, maxDens);

    // Calculate (Q|ji) once.
    Eigen::MatrixXd Qji_a = Eigen::MatrixXd(no_J_a * no_I_a, nx);
    Eigen::MatrixXd Qji_b = Eigen::MatrixXd(no_J_b * no_I_b, nx);

    // Calculate (Q|\tilde ji) for each set and guess vector.
    std::vector<std::vector<std::vector<Eigen::MatrixXd>>> Qji_t(_nSet);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        Qji_t[iSet] = std::vector<std::vector<Eigen::MatrixXd>>(_nGuess);
        for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
          Qji_t[iSet][iGuess] = std::vector<Eigen::MatrixXd>(2);
          Qji_t[iSet][iGuess][0] = Eigen::MatrixXd(no_J_a * no_I_a, nx);
          Qji_t[iSet][iGuess][1] = Eigen::MatrixXd(no_J_b * no_I_b, nx);
        }
      }
    }

// First contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> kappaJ_lambdaI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::MatrixXd kappaJ_iI_a = kappaJ_lambdaI * C_I.alpha.middleCols(0, no_I_a);
        Eigen::MatrixXd kappaJ_iI_b = kappaJ_lambdaI * C_I.beta.middleCols(0, no_I_b);
        Eigen::Map<Eigen::MatrixXd> ji_a(Qji_a.col(P_all).data(), no_J_a, no_I_a);
        Eigen::Map<Eigen::MatrixXd> ji_b(Qji_b.col(P_all).data(), no_J_b, no_I_b);
        ji_a = C_J.alpha.middleCols(0, no_J_a).transpose() * kappaJ_iI_a;
        ji_b = C_J.beta.middleCols(0, no_J_b).transpose() * kappaJ_iI_b;
        if (TDDFT) {
          for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
            for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
              Eigen::Map<Eigen::MatrixXd> ji_t_a(Qji_t[iSet][iGuess][0].col(P_all).data(), no_J_a, no_I_a);
              Eigen::Map<Eigen::MatrixXd> ji_t_b(Qji_t[iSet][iGuess][1].col(P_all).data(), no_J_b, no_I_b);
              ji_t_a = C_vj[iSet][iGuess][0].middleCols(0, no_J_a).transpose() * kappaJ_iI_a;
              ji_t_b = C_vj[iSet][iGuess][1].middleCols(0, no_J_b).transpose() * kappaJ_iI_b;
            }
          }
        }
      }
    }

    // Contract with inverse metric matrix.
    Qji_a *= (*invM);
    Qji_b *= (*invM);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
          Qji_t[iSet][iGuess][0] *= (*invM);
          Qji_t[iSet][iGuess][1] *= (*invM);
        }
      }
    }

// Second contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> nuJ_muI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::Map<Eigen::MatrixXd> ji_a(Qji_a.col(P_all).data(), no_J_a, no_I_a);
        Eigen::Map<Eigen::MatrixXd> ji_b(Qji_b.col(P_all).data(), no_J_b, no_I_b);
        for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
          for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
            Eigen::Ref<Eigen::MatrixXd> F_a = F_mui[iThread][iSet][iGuess][0];
            Eigen::Ref<Eigen::MatrixXd> F_b = F_mui[iThread][iSet][iGuess][1];
            if (!TDDFT) {
              // Contribution of A matrix.
              F_a.noalias() += nuJ_muI.transpose() * C_vj[iSet][iGuess][0] * ji_a;
              F_b.noalias() += nuJ_muI.transpose() * C_vj[iSet][iGuess][1] * ji_b;
            }
            else {
              // Contribution of (A pm B) matrix.
              Eigen::Map<Eigen::MatrixXd> ji_t_a(Qji_t[iSet][iGuess][0].col(P_all).data(), no_J_a, no_I_a);
              Eigen::Map<Eigen::MatrixXd> ji_t_b(Qji_t[iSet][iGuess][1].col(P_all).data(), no_J_b, no_I_b);
              F_a.noalias() += nuJ_muI.transpose() *
                               (C_vj[iSet][iGuess][0] * ji_a + _pm[iSet] * C_J.alpha.middleCols(0, no_J_a) * ji_t_a);
              F_b.noalias() += nuJ_muI.transpose() *
                               (C_vj[iSet][iGuess][1] * ji_b + _pm[iSet] * C_J.beta.middleCols(0, no_J_b) * ji_t_b);
            }
          }
        }
      }
    }
    Timings::timeTaken("LRSCF -   Sigmavector:   RI-K");

    // Add HF exchange to Fock matrix.
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    Eigen::MatrixXd temp_a = S_I * C_I.alpha.middleCols(0, no_I_a);
    Eigen::MatrixXd temp_b = S_I * C_I.beta.middleCols(0, no_I_b);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Eigen::Ref<Eigen::MatrixXd> Fx_a = F_mui[0][iSet][iGuess][0];
        Eigen::Ref<Eigen::MatrixXd> Fx_b = F_mui[0][iSet][iGuess][1];
        for (unsigned iThread = 1; iThread < this->_nThreads; ++iThread) {
          Fx_a.noalias() += F_mui[iThread][iSet][iGuess][0];
          Fx_b.noalias() += F_mui[iThread][iSet][iGuess][1];
          F_mui[iThread][iSet][iGuess][0].setZero();
          F_mui[iThread][iSet][iGuess][1].setZero();
        }
        (*fock)[iSet][iGuess].alpha += _hfExchangeRatio * temp_a * Fx_a.transpose();
        (*fock)[iSet][iGuess].beta += _hfExchangeRatio * temp_b * Fx_b.transpose();
        F_mui[0][iSet][iGuess][0].setZero();
        F_mui[0][iSet][iGuess][1].setZero();
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  } /* Exchange */

  // Long-range exchange.
  if (_densFitLRK && _lrExchangeRatio > 0.0) {
    OutputControl::dOut << " **** Performing RI-ErfK ****" << std::endl;
    Timings::takeTime("LRSCF -   Sigmavector: RI-LRK");

    // Determine common basis purpose.
    auto basispurpose = Options::BASIS_PURPOSES::AUX_CORREL;
    if (this->_lrscf[I]->getLRSCFSettings().densFitLRK == Options::DENS_FITS::ACD) {
      basispurpose = Options::BASIS_PURPOSES::ERF_ATOMIC_CHOLESKY;
    }
    else if (this->_lrscf[I]->getLRSCFSettings().densFitLRK == Options::DENS_FITS::ACCD) {
      basispurpose = Options::BASIS_PURPOSES::ERF_ATOMIC_COMPACT_CHOLESKY;
    }

    // Initialize aux basis controller
    auto auxBasisControllerI = this->_lrscf[I]->getBasisController(basispurpose);
    auto auxBasisControllerJ = this->_lrscf[J]->getBasisController(basispurpose);

    // Bases for looper.
    std::shared_ptr<BasisController> superBas = nullptr;
    std::shared_ptr<BasisController> superAuxBas = nullptr;

    if (I == J) {
      superBas = basisControllerI;
      superAuxBas = auxBasisControllerI;
    }
    else {
      BasisFunctionMapper basisMapper(basisControllerI);
      BasisFunctionMapper auxBasisMapper(auxBasisControllerI);
      superBas = basisMapper.getCombinedBasis(basisControllerJ);
      superAuxBas = auxBasisMapper.getCombinedBasis(auxBasisControllerJ);
    }
    unsigned nx = superAuxBas->getNBasisFunctions();
    unsigned nxs = superAuxBas->getReducedNBasisFunctions();
    auto& auxShells = superAuxBas->getBasis();

    RI_J_IntegralController riints(superBas, superAuxBas, nullptr, LIBINT_OPERATOR::erf_coulomb, _mu);

    auto invM = _lrscf[I]->getInverseErfMetric();
    if (invM == nullptr) {
      invM = std::make_shared<Eigen::MatrixXd>(riints.getInverseM());
      _lrscf[I]->setInverseErfMetric(invM);
    }

    TwoElecThreeCenterCalculator intCalculator(LIBINT_OPERATOR::erf_coulomb, _mu, basisControllerJ, basisControllerI,
                                               superAuxBas, prescreeningThreshold, maxDens);

    // Calculate (Q|ji) once.
    Eigen::MatrixXd Qji_a = Eigen::MatrixXd(no_J_a * no_I_a, nx);
    Eigen::MatrixXd Qji_b = Eigen::MatrixXd(no_J_b * no_I_b, nx);

    // Calculate (Q|\tilde ji) for each set and guess vector.
    std::vector<std::vector<std::vector<Eigen::MatrixXd>>> Qji_t(_nSet);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        Qji_t[iSet] = std::vector<std::vector<Eigen::MatrixXd>>(_nGuess);
        for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
          Qji_t[iSet][iGuess] = std::vector<Eigen::MatrixXd>(2);
          Qji_t[iSet][iGuess][0] = Eigen::MatrixXd(no_J_a * no_I_a, nx);
          Qji_t[iSet][iGuess][1] = Eigen::MatrixXd(no_J_b * no_I_b, nx);
        }
      }
    }

// First contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> kappaJ_lambdaI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::MatrixXd kappaJ_iI_a = kappaJ_lambdaI * C_I.alpha.middleCols(0, no_I_a);
        Eigen::MatrixXd kappaJ_iI_b = kappaJ_lambdaI * C_I.beta.middleCols(0, no_I_b);
        Eigen::Map<Eigen::MatrixXd> ji_a(Qji_a.col(P_all).data(), no_J_a, no_I_a);
        Eigen::Map<Eigen::MatrixXd> ji_b(Qji_b.col(P_all).data(), no_J_b, no_I_b);
        ji_a = C_J.alpha.middleCols(0, no_J_a).transpose() * kappaJ_iI_a;
        ji_b = C_J.beta.middleCols(0, no_J_b).transpose() * kappaJ_iI_b;
        if (TDDFT) {
          for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
            for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
              Eigen::Map<Eigen::MatrixXd> ji_t_a(Qji_t[iSet][iGuess][0].col(P_all).data(), no_J_a, no_I_a);
              Eigen::Map<Eigen::MatrixXd> ji_t_b(Qji_t[iSet][iGuess][1].col(P_all).data(), no_J_b, no_I_b);
              ji_t_a = C_vj[iSet][iGuess][0].middleCols(0, no_J_a).transpose() * kappaJ_iI_a;
              ji_t_b = C_vj[iSet][iGuess][1].middleCols(0, no_J_b).transpose() * kappaJ_iI_b;
            }
          }
        }
      }
    }

    // Contract with inverse metric matrix.
    Qji_a *= (*invM);
    Qji_b *= (*invM);
    if (TDDFT) {
      for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
        for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
          Qji_t[iSet][iGuess][0] *= (*invM);
          Qji_t[iSet][iGuess][1] *= (*invM);
        }
      }
    }

// Second contraction.
#pragma omp parallel for schedule(dynamic)
    for (size_t iShell = 0; iShell < nxs; ++iShell) {
      unsigned iThread = omp_get_thread_num();
      auto integrals = intCalculator.calculateIntegrals(iShell, iThread);

      unsigned P_in_iShell = auxShells[iShell]->getNContracted();
      for (unsigned P = 0; P < P_in_iShell; ++P) {
        unsigned P_all = superAuxBas->extendedIndex(iShell) + P;
        Eigen::Map<Eigen::MatrixXd> nuJ_muI(integrals.col(P).data(), nb_J, nb_I);
        Eigen::Map<Eigen::MatrixXd> ji_a(Qji_a.col(P_all).data(), no_J_a, no_I_a);
        Eigen::Map<Eigen::MatrixXd> ji_b(Qji_b.col(P_all).data(), no_J_b, no_I_b);
        for (unsigned iSet = 0; iSet < _nSet; ++iSet) {
          for (unsigned iGuess = 0; iGuess < _nGuess; ++iGuess) {
            Eigen::Ref<Eigen::MatrixXd> F_a = F_mui[iThread][iSet][iGuess][0];
            Eigen::Ref<Eigen::MatrixXd> F_b = F_mui[iThread][iSet][iGuess][1];
            if (!TDDFT) {
              // Contribution of A matrix.
              F_a.noalias() += nuJ_muI.transpose() * C_vj[iSet][iGuess][0] * ji_a;
              F_b.noalias() += nuJ_muI.transpose() * C_vj[iSet][iGuess][1] * ji_b;
            }
            else {
              // Contribution of (A pm B) matrix.
              Eigen::Map<Eigen::MatrixXd> ji_t_a(Qji_t[iSet][iGuess][0].col(P_all).data(), no_J_a, no_I_a);
              Eigen::Map<Eigen::MatrixXd> ji_t_b(Qji_t[iSet][iGuess][1].col(P_all).data(), no_J_b, no_I_b);
              F_a.noalias() += nuJ_muI.transpose() *
                               (C_vj[iSet][iGuess][0] * ji_a + _pm[iSet] * C_J.alpha.middleCols(0, no_J_a) * ji_t_a);
              F_b.noalias() += nuJ_muI.transpose() *
                               (C_vj[iSet][iGuess][1] * ji_b + _pm[iSet] * C_J.beta.middleCols(0, no_J_b) * ji_t_b);
            }
          }
        }
      }
    }
    Timings::timeTaken("LRSCF -   Sigmavector: RI-LRK");

    // Add LR exchange to Fock matrix.
    Timings::takeTime("LRSCF - Add/Sym Fock Matrices");
    Eigen::MatrixXd temp_a = S_I * C_I.alpha.middleCols(0, no_I_a);
    Eigen::MatrixXd temp_b = S_I * C_I.beta.middleCols(0, no_I_b);
    for (unsigned iSet = 0; iSet < this->_nSet; ++iSet) {
      for (unsigned iGuess = 0; iGuess < this->_nGuess; ++iGuess) {
        Eigen::Ref<Eigen::MatrixXd> Fx_a = F_mui[0][iSet][iGuess][0];
        Eigen::Ref<Eigen::MatrixXd> Fx_b = F_mui[0][iSet][iGuess][1];
        for (unsigned iThread = 1; iThread < this->_nThreads; ++iThread) {
          Fx_a.noalias() += F_mui[iThread][iSet][iGuess][0];
          Fx_b.noalias() += F_mui[iThread][iSet][iGuess][1];
        }
        (*fock)[iSet][iGuess].alpha += _lrExchangeRatio * temp_a * Fx_a.transpose();
        (*fock)[iSet][iGuess].beta += _lrExchangeRatio * temp_b * Fx_b.transpose();
      }
    }
    Timings::timeTaken("LRSCF - Add/Sym Fock Matrices");
  } /* LR-Exchange */

  return fock;
}

template class RIExchangeSigmavector<Options::SCF_MODES::RESTRICTED>;
template class RIExchangeSigmavector<Options::SCF_MODES::UNRESTRICTED>;
} // namespace Serenity
