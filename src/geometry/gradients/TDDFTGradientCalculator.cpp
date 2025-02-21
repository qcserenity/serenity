/**
 * @file   TDDFTGradientCalculator.cpp
 *
 * @date   Mar 11, 2024
 * @author Anton Rikus
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
#include "geometry/gradients/TDDFTGradientCalculator.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ScalarOperatorToMatrixAdder.h"
#include "data/matrices/MatrixInBasis.h"
#include "dft/Functional.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "dft/functionals/wrappers/PartialDerivatives.h" // FunctionalData
#include "geometry/Geometry.h"
#include "geometry/GeometryAdderFactory.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "integrals/OneElectronIntegralDerivativeCalculator.h"
#include "integrals/RIIntegralDerivativeCalculator.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "io/FormattedOutputStream.h"
#include "math/Matrix.h"
#include "misc/Timing.h"
#include "postHF/LRSCF/Kernel/Kernel.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/Sigmavectors/CoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/ExchangeSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/FockSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/HyperkernelSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/KernelSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/RICoulombSigmavector.h"
#include "postHF/LRSCF/Sigmavectors/RI/RIExchangeSigmavector.h"
#include "postHF/LRSCF/Tools/LRSCFSetup.h"
#include "postHF/LRSCF/Tools/NonlinearResponseSolver.h"
#include "potentials/bundles/PotentialBundle.h"
#include "settings/ElectronicStructureOptions.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h" // LRSCFTaskSettings

namespace Serenity {

template<Options::SCF_MODES SCFMode>
TDDFTGradientCalculator<SCFMode>::TDDFTGradientCalculator(std::shared_ptr<LRSCFController<SCFMode>> lrscf, double hypThresh)
  // HF functional as default so that the HF exchange ratio is correct (1 for CIS/TDHF)
  : _lrscf(lrscf), _hypThresh(hypThresh), _func(resolveFunctional(CompositeFunctionals::XCFUNCTIONALS::HF)) {
  // do not calculate kernel and hyperkernel contributions if the system's method is not DFT
  if (_lrscf->getSys()->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    // do not calculate kernel and hyperkernel contributions if the system uses the HF functional
    if (_lrscf->getSys()->getSettings().dft.functional != CompositeFunctionals::XCFUNCTIONALS::HF) {
      _usesXC = true;
      std::vector<std::shared_ptr<SystemController>> systems = {_lrscf->getSys()};
      _kernel = std::make_shared<Kernel<SCFMode>>(systems, (std::vector<std::shared_ptr<SystemController>>){},
                                                  _lrscf->getLRSCFSettings(), 1, Options::GRID_PURPOSES::SMALL);
      std::make_shared<Kernel<SCFMode>>(systems, (std::vector<std::shared_ptr<SystemController>>){},
                                        _lrscf->getLRSCFSettings(), 1, Options::GRID_PURPOSES::SMALL);

      // At this point, clearly an exchange--correlation functional is used. Which one? There are four options, the
      // composite functional from the system settings, the custom functional from the system settings, the composite
      // functional from the LRSCFTaskSettings and the custom functional from the LRSCFTaskSettings. Since the specifity
      // to the task of calculating TDDFT gradients increases in this order, the custom LRSCF functional has the highest
      // priority. To summarize: LRSCFTaskSettings before system settings and custom before composite.

      // use the custom LRSCF functional if there is one
      if (_lrscf->getLRSCFSettings().customFunc.basicFunctionals.size()) {
        _func = Functional(_lrscf->getLRSCFSettings().customFunc);
      }
      // use the composite LRSCF functional if there is one
      else if (_lrscf->getLRSCFSettings().func != CompositeFunctionals::XCFUNCTIONALS::NONE) {
        _func = resolveFunctional(this->_lrscf->getLRSCFSettings().func);
      }
      // no LRSCF functional -> use the system settings
      else if (_lrscf->getSys()->getSettings().customFunc.basicFunctionals.size()) {
        _func = Functional(_lrscf->getSys()->getSettings().customFunc);
      }
      // system settings composite functional at last
      else {
        _func = resolveFunctional(this->_lrscf->getSys()->getSettings().dft.functional);
      }
    }
  }
  if (!_usesXC) {
    printSectionTitle("TDHF/CIS Gradients");
  }
  else {
    printSectionTitle("TDDFT/TDA Gradients");
    FunctionalLibrary<SCFMode> flib(_lrscf->getLRSCFSettings().grid.blocksize);

    _basisFunctionOnGridController = BasisFunctionOnGridControllerFactory::produce(
        _lrscf->getLRSCFSettings().grid.blocksize, _lrscf->getLRSCFSettings().grid.basFuncRadialThreshold,
        // GGA gradients require the second derivative of basis functions at each gridpoint, while for LDA gradients the
        // first derivative is sufficient
        (_kernel->isGGA() ? 2 : 1), _lrscf->getSys()->getBasisController(),
        AtomCenteredGridControllerFactory::produce(
            GeometryAdderFactory::produce({_lrscf->getSys()}, _lrscf->getLRSCFSettings().subsystemgrid),
            _lrscf->getLRSCFSettings().grid, Options::GRID_PURPOSES::SMALL));
    auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
        _basisFunctionOnGridController, _lrscf->getLRSCFSettings().grid.blockAveThreshold);
    auto densityMatrixController =
        std::make_shared<DensityMatrixController<SCFMode>>(_lrscf->getSys()->template getActiveOrbitalController<SCFMode>(),
                                                           _lrscf->getSys()->template getNOccupiedOrbitals<SCFMode>());
    _gridToMatrix = std::make_shared<ScalarOperatorToMatrixAdder<SCFMode>>(
        _basisFunctionOnGridController, this->_lrscf->getLRSCFSettings().grid.blockAveThreshold);

    auto densOnGridController =
        std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densityOnGridCalculator, densityMatrixController);
    // calculate the functional data once and then store it
    _funcData = std::make_shared<FunctionalData<SCFMode>>(
        flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _func, densOnGridController, 3));
    _hyperSig = std::make_shared<HyperkernelSigmavector<SCFMode>>(_lrscf, _funcData, _hypThresh);
  }
  // set up nonlinear Sigmacalculator to calculate matrix-vector products for the CPHF left-hand-side
  _nlSolver = [&](Eigen::MatrixXd& guessVectors, Eigen::VectorXd eigenvalues) {
    (void)eigenvalues;
    auto sigma = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(guessVectors.rows(), guessVectors.cols()));

    const std::vector<Eigen::MatrixXd>& guess = {guessVectors};

    std::unique_ptr<Sigmavector<SCFMode>> D, J, K, DFK, Ke;

    std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf = {_lrscf};
    D = std::make_unique<FockSigmavector<SCFMode>>(lrscf, guess);
    if (_lrscf->getLRSCFSettings().densFitJ == Options::DENS_FITS::NONE) {
      J = std::make_unique<CoulombSigmavector<SCFMode>>(lrscf, guess);
    }
    else {
      J = std::make_unique<RICoulombSigmavector<SCFMode>>(lrscf, guess);
    }
    std::vector<int> pm(1, 1);
    bool densFitK = (_lrscf->getLRSCFSettings().densFitK != Options::DENS_FITS::NONE);
    bool densFitLRK = (_lrscf->getLRSCFSettings().densFitLRK != Options::DENS_FITS::NONE);
    if (densFitK || densFitLRK) {
      DFK = std::make_unique<RIExchangeSigmavector<SCFMode>>(lrscf, guess, pm, densFitK, densFitLRK);
    }
    if (!densFitK || !densFitLRK) {
      K = std::make_unique<ExchangeSigmavector<SCFMode>>(lrscf, guess, pm, densFitK, densFitLRK);
    }

    double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
    (*sigma) += D->getSigma()[0];
    (*sigma) += scfFactor * 2.0 * J->getSigma()[0];
    if (K)
      (*sigma) -= K->getSigma()[0];
    if (DFK)
      (*sigma) -= DFK->getSigma()[0];
    if (_usesXC) {
      Timings::takeTime("LRSCF Gradients -        XC Contribution");
      Ke = std::make_unique<KernelSigmavector<SCFMode>>(lrscf, guess, _kernel);
      (*sigma) += scfFactor * 2.0 * Ke->getSigma()[0];
      Timings::timeTaken("LRSCF Gradients -        XC Contribution");
    }
    return sigma;
  };

  _excGradList = _lrscf->getLRSCFSettings().excGradList;
  Timings::takeTime("LRSCF Gradients -   ground-state Contrib.");
  if (_lrscf->getSys()->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    _groundstateGradient =
        _lrscf->getSys()->template getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::HF>()->getGradients();
  }
  else if (_lrscf->getSys()->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) {
    _groundstateGradient =
        _lrscf->getSys()->template getPotentials<SCFMode, Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>()->getGradients();
  }
  Timings::timeTaken("LRSCF Gradients -   ground-state Contrib.");
}

template<Options::SCF_MODES SCFMode>
MatrixInBasis<SCFMode> TDDFTGradientCalculator<SCFMode>::hPlus(MatrixInBasis<SCFMode> V) {
  MatrixInBasis<SCFMode> resultingSigma(_lrscf->getSys()->getBasisController());
  const CoefficientMatrix<SCFMode>& C = _lrscf->getCoefficients();
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf = {_lrscf};

  std::unique_ptr<Sigmavector<SCFMode>> J, K, DFK, Ke;

  if (_lrscf->getLRSCFSettings().densFitJ == Options::DENS_FITS::NONE) {
    J = std::make_unique<CoulombSigmavector<SCFMode>>(lrscf);
  }
  else {
    J = std::make_unique<RICoulombSigmavector<SCFMode>>(lrscf);
  }
  bool densFitK = (_lrscf->getLRSCFSettings().densFitK != Options::DENS_FITS::NONE);
  bool densFitLRK = (_lrscf->getLRSCFSettings().densFitLRK != Options::DENS_FITS::NONE);
  if (densFitK || densFitLRK) {
    // this doesn't work - need to look closer into RI exchange terms. calcF of RIExchangeSigmavector doesn't actually
    // use the densitymatrix provided, but relies on the (MO) guessvectors to exist
    DFK = std::make_unique<RIExchangeSigmavector<SCFMode>>(lrscf, (std::vector<int>){1}, densFitK, densFitLRK);
  }
  if (!densFitK || !densFitLRK) {
    K = std::make_unique<ExchangeSigmavector<SCFMode>>(lrscf, (std::vector<int>){1}, densFitK, densFitLRK);
  }
  if (_usesXC) {
    Timings::takeTime("LRSCF Gradients -        XC Contribution");
    Ke = std::make_unique<KernelSigmavector<SCFMode>>(lrscf, _kernel);
    Timings::timeTaken("LRSCF Gradients -        XC Contribution");
  }

  MatrixInBasis<SCFMode> P(_lrscf->getSys()->getBasisController());
  for_spin(P, C, V) {
    P_spin = C_spin * V_spin * (Eigen::MatrixXd)C_spin.transpose();
  };
  std::vector<std::vector<MatrixInBasis<SCFMode>>> pSet = {{P}};
  auto pPtr = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(pSet);
  MatrixInBasis<SCFMode> JMat = (*(J->calcF(0, 0, std::move(pPtr))))[0][0];
  MatrixInBasis<SCFMode> KMat(_lrscf->getSys()->getBasisController());
  if (K) {
    auto pPtr = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(pSet);
    KMat = (*(K->calcF(0, 0, std::move(pPtr))))[0][0];
  }
  if (DFK) {
    auto pPtr = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(pSet);
    KMat += (*(DFK->calcF(0, 0, std::move(pPtr))))[0][0];
  }
  double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
  if (_usesXC) {
    Timings::takeTime("LRSCF Gradients -        XC Contribution");
    auto pPtr = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(pSet);
    MatrixInBasis<SCFMode> KeMat = (*(Ke->calcF(0, 0, std::move(pPtr))))[0][0];
    for_spin(resultingSigma, KeMat, C) {
      resultingSigma_spin = 2.0 * scfFactor * C_spin.transpose() * KeMat_spin * C_spin;
    };
    Timings::timeTaken("LRSCF Gradients -        XC Contribution");
  }
  for_spin(resultingSigma, JMat, KMat, C) {
    resultingSigma_spin += C_spin.transpose() * (2.0 * scfFactor * JMat_spin - KMat_spin) * C_spin;
  };

  return resultingSigma;
}

template<Options::SCF_MODES SCFMode>
MatrixInBasis<SCFMode> TDDFTGradientCalculator<SCFMode>::hMinus(MatrixInBasis<SCFMode> V) {
  MatrixInBasis<SCFMode> resultingSigma(_lrscf->getSys()->getBasisController());
  const CoefficientMatrix<SCFMode>& C = _lrscf->getCoefficients();
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf = {_lrscf};
  std::unique_ptr<Sigmavector<SCFMode>> K, DFK;
  bool densFitK = (_lrscf->getLRSCFSettings().densFitK != Options::DENS_FITS::NONE);
  bool densFitLRK = (_lrscf->getLRSCFSettings().densFitLRK != Options::DENS_FITS::NONE);
  // K and DFK can coexist in case long-range exchange is handled by RI while short-range exchange is not (or vice versa)
  if (densFitK || densFitLRK) {
    DFK = std::make_unique<RIExchangeSigmavector<SCFMode>>(lrscf, (std::vector<int>){-1}, densFitK, densFitLRK);
  }
  if (!densFitK || !densFitLRK) {
    K = std::make_unique<ExchangeSigmavector<SCFMode>>(lrscf, (std::vector<int>){-1}, densFitK, densFitLRK);
  }

  MatrixInBasis<SCFMode> P(_lrscf->getSys()->getBasisController());
  for_spin(P, C, V) {
    P_spin = C_spin * V_spin * (Eigen::MatrixXd)C_spin.transpose();
  };
  std::vector<std::vector<MatrixInBasis<SCFMode>>> pSet = {{P}};
  MatrixInBasis<SCFMode> KMat(_lrscf->getSys()->getBasisController());
  if (K) {
    auto pPtr = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(pSet);
    KMat = (*(K->calcF(0, 0, std::move(pPtr))))[0][0];
  }
  if (DFK) {
    auto pPtr = std::make_unique<std::vector<std::vector<MatrixInBasis<SCFMode>>>>(pSet);
    KMat += (*(DFK->calcF(0, 0, std::move(pPtr))))[0][0];
  }
  for_spin(resultingSigma, KMat, C) {
    resultingSigma_spin = C_spin.transpose() * KMat_spin * C_spin;
  };

  return resultingSigma;
}

template<>
const Eigen::MatrixXd
TDDFTGradientCalculator<RESTRICTED>::evaluateFullTwoCenterGradient(const MatrixInBasis<RESTRICTED>& PAO,
                                                                   const MatrixInBasis<RESTRICTED>& D,
                                                                   const MatrixInBasis<RESTRICTED>& xpyAO,
                                                                   const MatrixInBasis<RESTRICTED>& xmyAO) {
  Eigen::MatrixXd twoCenterGradient = Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3);
  double hfExchangeRatio = _func.getHfExchangeRatio();

  const std::function<double(const unsigned int, const unsigned int, const unsigned int&, const unsigned int&)> twoParticleDiffDens =
      [&xpyAO, &xmyAO, &D, &PAO, &hfExchangeRatio](const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                                   const unsigned int& d) {
        // PAO and D are already summed over spins, while xpyAO is not
        // factor 2 instead of 4 because restricted X+Y is normalized to sqrt(2) * unrestricted X+Y
        return (PAO(a, b) * D(c, d) + 2 * xpyAO(a, b) * xpyAO(c, d) -
                hfExchangeRatio * 0.25 * (PAO(a, c) * D(b, d) + PAO(a, d) * D(b, c)) -
                0.5 * hfExchangeRatio *
                    (xpyAO(a, d) * xpyAO(c, b) + xpyAO(a, c) * xpyAO(b, d) - xmyAO(a, d) * xmyAO(c, b) +
                     xmyAO(a, c) * xmyAO(b, d)));
      };

  std::vector<unsigned int> atomIndicesOfBasis = _lrscf->getSys()->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  std::vector<Eigen::MatrixXd> coulombGradientPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(twoCenterGradient.rows(), 3));
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _lrscf->getSys()->getAtomCenteredBasisController(),
                                    _lrscf->getSys()->getAtomCenteredBasisController()->getPrescreeningThreshold());
  auto const looperFunction = [&atomIndicesOfBasis, &twoParticleDiffDens, &coulombGradientPriv](
                                  const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                  const unsigned int& d, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    std::vector<unsigned> atomMaps = {atomIndicesOfBasis[a], atomIndicesOfBasis[b], atomIndicesOfBasis[c],
                                      atomIndicesOfBasis[d]};
    double symmetrizedFactor = twoParticleDiffDens(a, b, c, d) + twoParticleDiffDens(a, b, d, c) +
                               twoParticleDiffDens(b, a, c, d) + twoParticleDiffDens(b, a, d, c) +
                               twoParticleDiffDens(c, d, a, b) + twoParticleDiffDens(c, d, b, a) +
                               twoParticleDiffDens(d, c, a, b) + twoParticleDiffDens(d, c, b, a);
    for (unsigned iCol = 0; iCol < 4; iCol++) {
      for (unsigned iDirection = 0; iDirection < 3; iDirection++) {
        coulombGradientPriv[threadID](atomMaps[iCol], iDirection) +=
            symmetrizedFactor * intValues.data()[3 * iCol + iDirection];
      }
    }
  };
  looper.loop(looperFunction, 100.0, true);
  for (Eigen::MatrixXd grad : coulombGradientPriv) {
    twoCenterGradient += grad;
  }
  return twoCenterGradient;
}

template<>
const Eigen::MatrixXd
TDDFTGradientCalculator<UNRESTRICTED>::evaluateFullTwoCenterGradient(const MatrixInBasis<UNRESTRICTED>& PAO,
                                                                     const MatrixInBasis<UNRESTRICTED>& D,
                                                                     const MatrixInBasis<UNRESTRICTED>& xpyAO,
                                                                     const MatrixInBasis<UNRESTRICTED>& xmyAO) {
  Eigen::MatrixXd twoCenterGradient = Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3);
  double hfExchangeRatio = _func.getHfExchangeRatio();

  const std::function<double(const unsigned int&, const unsigned int&, const unsigned int&, const unsigned int&)> twoParticleDiffDens =
      [&xpyAO, &xmyAO, &D, &PAO, &hfExchangeRatio](const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                                   const unsigned int& d) {
        double gamma = PAO.alpha(a, b) * D.beta(c, d) + PAO.beta(a, b) * D.alpha(c, d) +
                       xpyAO.alpha(a, b) * xpyAO.beta(c, d) + xpyAO.beta(a, b) * xpyAO.alpha(c, d);
        for_spin(PAO, D, xpyAO, xmyAO) {
          gamma += PAO_spin(a, b) * D_spin(c, d) + xpyAO_spin(a, b) * xpyAO_spin(c, d) -
                   0.5 * hfExchangeRatio *
                       (PAO_spin(b, c) * D_spin(a, d) + PAO_spin(b, d) * D_spin(a, c) +
                        xpyAO_spin(a, d) * xpyAO_spin(c, b) + xpyAO_spin(a, c) * xpyAO_spin(b, d) -
                        xmyAO_spin(a, d) * xmyAO_spin(c, b) + xmyAO_spin(a, c) * xmyAO_spin(b, d));
        };
        return gamma;
      };

  std::vector<unsigned int> atomIndicesOfBasis = _lrscf->getSys()->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  std::vector<Eigen::MatrixXd> coulombGradientPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(twoCenterGradient.rows(), 3));
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _lrscf->getSys()->getAtomCenteredBasisController(),
                                    _lrscf->getSys()->getAtomCenteredBasisController()->getPrescreeningThreshold());
  auto const looperFunction = [&atomIndicesOfBasis, &twoParticleDiffDens, &coulombGradientPriv](
                                  const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                  const unsigned int& d, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    std::vector<unsigned> atomMaps = {atomIndicesOfBasis[a], atomIndicesOfBasis[b], atomIndicesOfBasis[c],
                                      atomIndicesOfBasis[d]};
    double symmetrizedFactor = twoParticleDiffDens(a, b, c, d) + twoParticleDiffDens(a, b, d, c) +
                               twoParticleDiffDens(b, a, c, d) + twoParticleDiffDens(b, a, d, c) +
                               twoParticleDiffDens(c, d, a, b) + twoParticleDiffDens(c, d, b, a) +
                               twoParticleDiffDens(d, c, a, b) + twoParticleDiffDens(d, c, b, a);
    for (unsigned iCol = 0; iCol < 4; iCol++) {
      for (unsigned iDirection = 0; iDirection < 3; iDirection++) {
        coulombGradientPriv[threadID](atomMaps[iCol], iDirection) +=
            symmetrizedFactor * intValues.data()[3 * iCol + iDirection];
      }
    }
  };
  looper.loop(looperFunction, 100.0, true);
  for (Eigen::MatrixXd grad : coulombGradientPriv) {
    twoCenterGradient += grad;
  }
  return twoCenterGradient;
}

template<>
const Eigen::MatrixXd
TDDFTGradientCalculator<RESTRICTED>::evaluateFullExchangeGradient(const MatrixInBasis<RESTRICTED>& PAO,
                                                                  const MatrixInBasis<RESTRICTED>& D,
                                                                  const MatrixInBasis<RESTRICTED>& xpyAO,
                                                                  const MatrixInBasis<RESTRICTED>& xmyAO) {
  Eigen::MatrixXd twoCenterGradient = Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3);
  double hfExchangeRatio = _func.getHfExchangeRatio();

  const std::function<double(const unsigned int, const unsigned int, const unsigned int&, const unsigned int&)> twoParticleDiffDens =
      [&xpyAO, &xmyAO, &D, &PAO, &hfExchangeRatio](const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                                   const unsigned int& d) {
        return (-hfExchangeRatio * 0.25 * (PAO(a, c) * D(b, d) + PAO(a, d) * D(b, c)) -
                0.5 * hfExchangeRatio *
                    (xpyAO(a, d) * xpyAO(c, b) + xpyAO(a, c) * xpyAO(b, d) - xmyAO(a, d) * xmyAO(c, b) +
                     xmyAO(a, c) * xmyAO(b, d)));
      };

  std::vector<unsigned int> atomIndicesOfBasis = _lrscf->getSys()->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  std::vector<Eigen::MatrixXd> coulombGradientPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(twoCenterGradient.rows(), 3));
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _lrscf->getSys()->getAtomCenteredBasisController(),
                                    _lrscf->getSys()->getAtomCenteredBasisController()->getPrescreeningThreshold());
  auto const looperFunction = [&atomIndicesOfBasis, &twoParticleDiffDens, &coulombGradientPriv](
                                  const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                  const unsigned int& d, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    std::vector<unsigned> atomMaps = {atomIndicesOfBasis[a], atomIndicesOfBasis[b], atomIndicesOfBasis[c],
                                      atomIndicesOfBasis[d]};
    double symmetrizedFactor = twoParticleDiffDens(a, b, c, d) + twoParticleDiffDens(a, b, d, c) +
                               twoParticleDiffDens(b, a, c, d) + twoParticleDiffDens(b, a, d, c) +
                               twoParticleDiffDens(c, d, a, b) + twoParticleDiffDens(c, d, b, a) +
                               twoParticleDiffDens(d, c, a, b) + twoParticleDiffDens(d, c, b, a);
    for (unsigned iCol = 0; iCol < 4; iCol++) {
      for (unsigned iDirection = 0; iDirection < 3; iDirection++) {
        coulombGradientPriv[threadID](atomMaps[iCol], iDirection) +=
            symmetrizedFactor * intValues.data()[3 * iCol + iDirection];
      }
    }
  };
  looper.loop(looperFunction, 100.0, true);
  for (Eigen::MatrixXd grad : coulombGradientPriv) {
    twoCenterGradient += grad;
  }
  return twoCenterGradient;
}

template<>
const Eigen::MatrixXd
TDDFTGradientCalculator<UNRESTRICTED>::evaluateFullExchangeGradient(const MatrixInBasis<UNRESTRICTED>& PAO,
                                                                    const MatrixInBasis<UNRESTRICTED>& D,
                                                                    const MatrixInBasis<UNRESTRICTED>& xpyAO,
                                                                    const MatrixInBasis<UNRESTRICTED>& xmyAO) {
  Eigen::MatrixXd twoCenterGradient = Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3);
  double hfExchangeRatio = _func.getHfExchangeRatio();

  const std::function<double(const unsigned int&, const unsigned int&, const unsigned int&, const unsigned int&)> twoParticleDiffDens =
      [&xpyAO, &xmyAO, &D, &PAO, &hfExchangeRatio](const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                                   const unsigned int& d) {
        double gamma = 0;
        for_spin(PAO, D, xpyAO, xmyAO) {
          gamma += -0.5 * hfExchangeRatio *
                   (PAO_spin(b, c) * D_spin(a, d) + PAO_spin(b, d) * D_spin(a, c) +
                    xpyAO_spin(a, d) * xpyAO_spin(c, b) + xpyAO_spin(a, c) * xpyAO_spin(b, d) -
                    xmyAO_spin(a, d) * xmyAO_spin(c, b) + xmyAO_spin(a, c) * xmyAO_spin(b, d));
        };
        return gamma;
      };

  std::vector<unsigned int> atomIndicesOfBasis = _lrscf->getSys()->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  std::vector<Eigen::MatrixXd> coulombGradientPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(twoCenterGradient.rows(), 3));
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 1, _lrscf->getSys()->getAtomCenteredBasisController(),
                                    _lrscf->getSys()->getAtomCenteredBasisController()->getPrescreeningThreshold());
  auto const looperFunction = [&atomIndicesOfBasis, &twoParticleDiffDens, &coulombGradientPriv](
                                  const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                  const unsigned int& d, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    std::vector<unsigned> atomMaps = {atomIndicesOfBasis[a], atomIndicesOfBasis[b], atomIndicesOfBasis[c],
                                      atomIndicesOfBasis[d]};
    double symmetrizedFactor = twoParticleDiffDens(a, b, c, d) + twoParticleDiffDens(a, b, d, c) +
                               twoParticleDiffDens(b, a, c, d) + twoParticleDiffDens(b, a, d, c) +
                               twoParticleDiffDens(c, d, a, b) + twoParticleDiffDens(c, d, b, a) +
                               twoParticleDiffDens(d, c, a, b) + twoParticleDiffDens(d, c, b, a);
    for (unsigned iCol = 0; iCol < 4; iCol++) {
      for (unsigned iDirection = 0; iDirection < 3; iDirection++) {
        coulombGradientPriv[threadID](atomMaps[iCol], iDirection) +=
            symmetrizedFactor * intValues.data()[3 * iCol + iDirection];
      }
    }
  };
  looper.loop(looperFunction, 100.0, true);
  for (Eigen::MatrixXd grad : coulombGradientPriv) {
    twoCenterGradient += grad;
  }
  return twoCenterGradient;
}

template<>
const Eigen::MatrixXd
TDDFTGradientCalculator<RESTRICTED>::evaluateFullLRExchangeGradient(const MatrixInBasis<RESTRICTED>& PAO,
                                                                    const MatrixInBasis<RESTRICTED>& D,
                                                                    const MatrixInBasis<RESTRICTED>& xpyAO,
                                                                    const MatrixInBasis<RESTRICTED>& xmyAO) {
  TwoElecFourCenterIntLooper longRangeLooper(
      LIBINT_OPERATOR::erf_coulomb, 1, _lrscf->getSys()->getAtomCenteredBasisController(),
      _lrscf->getSys()->getAtomCenteredBasisController()->getPrescreeningThreshold(), _func.getRangeSeparationParameter());
  double lrExchangeRatio = _func.getLRExchangeRatio();

  const std::function<double(const unsigned int, const unsigned int, const unsigned int&, const unsigned int&)> gammalaLR =
      [&xpyAO, &xmyAO, &D, &PAO, &lrExchangeRatio](const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                                   const unsigned int& d) {
        return (-lrExchangeRatio * 0.25 * (PAO(a, c) * D(b, d) + PAO(a, d) * D(b, c)) +
                -0.5 * lrExchangeRatio *
                    (xpyAO(a, d) * xpyAO(c, b) + xpyAO(a, c) * xpyAO(b, d) - xmyAO(a, d) * xmyAO(c, b) +
                     xmyAO(a, c) * xmyAO(b, d)));
      };

  std::vector<unsigned int> atomIndicesOfBasis = _lrscf->getSys()->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  std::vector<Eigen::MatrixXd> coulombGradientPrivLR(omp_get_max_threads(),
                                                     Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3));
  Eigen::MatrixXd coulombGradientLR = Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3);
  auto const looperFunction = [&atomIndicesOfBasis, &gammalaLR, &coulombGradientPrivLR](
                                  const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                  const unsigned int& d, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    std::vector<unsigned> atomMaps = {atomIndicesOfBasis[a], atomIndicesOfBasis[b], atomIndicesOfBasis[c],
                                      atomIndicesOfBasis[d]};
    double symmetrizedFactor = gammalaLR(a, b, c, d) + gammalaLR(a, b, d, c) + gammalaLR(b, a, c, d) +
                               gammalaLR(b, a, d, c) + gammalaLR(c, d, a, b) + gammalaLR(c, d, b, a) +
                               gammalaLR(d, c, a, b) + gammalaLR(d, c, b, a);
    for (unsigned iDirection = 0; iDirection < 3; iDirection++) {
      for (unsigned iCol = 0; iCol < 4; iCol++) {
        coulombGradientPrivLR[threadID](atomMaps[iCol], iDirection) +=
            symmetrizedFactor * intValues.data()[3 * iCol + iDirection];
      }
    }
  };
  longRangeLooper.loop(looperFunction, 100.0, true);
  for (Eigen::MatrixXd grad : coulombGradientPrivLR) {
    coulombGradientLR += grad;
  }
  return coulombGradientLR;
}

template<>
const Eigen::MatrixXd
TDDFTGradientCalculator<UNRESTRICTED>::evaluateFullLRExchangeGradient(const MatrixInBasis<UNRESTRICTED>& PAO,
                                                                      const MatrixInBasis<UNRESTRICTED>& D,
                                                                      const MatrixInBasis<UNRESTRICTED>& xpyAO,
                                                                      const MatrixInBasis<UNRESTRICTED>& xmyAO) {
  TwoElecFourCenterIntLooper longRangeLooper(
      LIBINT_OPERATOR::erf_coulomb, 1, _lrscf->getSys()->getAtomCenteredBasisController(),
      _lrscf->getSys()->getAtomCenteredBasisController()->getPrescreeningThreshold(), _func.getRangeSeparationParameter());
  double lrExchangeRatio = _func.getLRExchangeRatio();

  const std::function<double(const unsigned int, const unsigned int, const unsigned int&, const unsigned int&)> gammalaLR =
      [&xpyAO, &xmyAO, &D, &PAO, &lrExchangeRatio](const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                                   const unsigned int& d) {
        double gamma = 0.0;
        for_spin(PAO, D, xpyAO, xmyAO) {
          gamma -= 0.5 * lrExchangeRatio *
                   (PAO_spin(b, c) * D_spin(a, d) + PAO_spin(b, d) * D_spin(a, c) +
                    xpyAO_spin(a, d) * xpyAO_spin(c, b) + xpyAO_spin(a, c) * xpyAO_spin(b, d) -
                    xmyAO_spin(a, d) * xmyAO_spin(c, b) + xmyAO_spin(a, c) * xmyAO_spin(b, d));
        };
        return gamma;
      };

  std::vector<unsigned int> atomIndicesOfBasis = _lrscf->getSys()->getAtomCenteredBasisController()->getAtomIndicesOfBasis();

  std::vector<Eigen::MatrixXd> coulombGradientPrivLR(omp_get_max_threads(),
                                                     Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3));
  Eigen::MatrixXd coulombGradientLR = Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3);
  auto const looperFunction = [&atomIndicesOfBasis, &gammalaLR, &coulombGradientPrivLR](
                                  const unsigned int& a, const unsigned int& b, const unsigned int& c,
                                  const unsigned int& d, const Eigen::VectorXd& intValues, const unsigned int& threadID) {
    std::vector<unsigned> atomMaps = {atomIndicesOfBasis[a], atomIndicesOfBasis[b], atomIndicesOfBasis[c],
                                      atomIndicesOfBasis[d]};
    double symmetrizedFactor = gammalaLR(a, b, c, d) + gammalaLR(a, b, d, c) + gammalaLR(b, a, c, d) +
                               gammalaLR(b, a, d, c) + gammalaLR(c, d, a, b) + gammalaLR(c, d, b, a) +
                               gammalaLR(d, c, a, b) + gammalaLR(d, c, b, a);
    for (unsigned iDirection = 0; iDirection < 3; iDirection++) {
      for (unsigned iCol = 0; iCol < 4; iCol++) {
        coulombGradientPrivLR[threadID](atomMaps[iCol], iDirection) +=
            symmetrizedFactor * intValues.data()[3 * iCol + iDirection];
      }
    }
  };
  longRangeLooper.loop(looperFunction, 1.0, true);
  for (Eigen::MatrixXd grad : coulombGradientPrivLR) {
    coulombGradientLR += grad;
  }
  return coulombGradientLR;
}

template<Options::SCF_MODES SCFMode>
const Eigen::MatrixXd
TDDFTGradientCalculator<SCFMode>::evaluateGradients(const MatrixInBasis<SCFMode>& PAO, const MatrixInBasis<SCFMode>& WAO,
                                                    const MatrixInBasis<SCFMode>& D, const MatrixInBasis<SCFMode>& xpyAO,
                                                    const MatrixInBasis<SCFMode>& xmyAO, double omega) {
  std::shared_ptr<AtomCenteredBasisController> basisController = _lrscf->getSys()->getAtomCenteredBasisController();
  const unsigned int nAtoms = _lrscf->getSys()->getNAtoms();

  OneElectronIntegralDerivativeCalculator intDerivCalculator(basisController, _lrscf->getSys()->getGeometry());
  Eigen::MatrixXd overlapGradient = -1.0 * intDerivCalculator.getOverlapDerivative(WAO.total());
  OutputControl::vOut << "  excited-state overlap contribution:\n" << overlapGradient << std::endl;
  Eigen::MatrixXd hCoreGradient = intDerivCalculator.getNucKinDerivative(PAO.total());
  OutputControl::vOut << "  excited-state Hcore contribution (nuclear + kinetic)\n" << hCoreGradient << std::endl;

  Eigen::MatrixXd vxcGradient = Eigen::MatrixXd::Zero(nAtoms, 3);
  Eigen::MatrixXd fxcGradient = Eigen::MatrixXd::Zero(nAtoms, 3);
  Eigen::MatrixXd vxcKernelGradient = Eigen::MatrixXd::Zero(nAtoms, 3);
  Eigen::MatrixXd hyperkernelGradient = Eigen::MatrixXd::Zero(nAtoms, 3);
  std::vector<unsigned int> atomIndicesOfBasis = _lrscf->getSys()->getAtomCenteredBasisController()->getAtomIndicesOfBasis();
  if (_usesXC) {
    Timings::takeTime("LRSCF Gradients -        XC Contribution");
    std::shared_ptr<SystemController> system = _lrscf->getSys();

    auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<SCFMode>>(
        _basisFunctionOnGridController, _lrscf->getLRSCFSettings().grid.blockAveThreshold);
    _gridToMatrix = std::make_shared<ScalarOperatorToMatrixAdder<SCFMode>>(
        _basisFunctionOnGridController, this->_lrscf->getLRSCFSettings().grid.blockAveThreshold);

    DensityMatrix<SCFMode> xDens(basisController);
    for_spin(xDens, xpyAO) {
      xDens_spin = xpyAO_spin + xpyAO_spin.transpose();
    };
    std::shared_ptr<DensityMatrixController<SCFMode>> xDmatController(std::make_shared<DensityMatrixController<SCFMode>>(xDens));

    auto xDensOnGridController =
        std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densityOnGridCalculator, xDmatController);

    // calculates  \rho(r) = \sum_{\mu\nu} \chi_\mu(r) \chi_\nu(r) D_{\mu\nu}  for each gridpoint and in this case
    // with symmetrized X as the density matrix
    auto& xDensOnGrid = xDensOnGridController->getDensityOnGrid();
    auto& xDensGradientOnGrid = xDensOnGridController->getDensityGradientOnGrid();

    DensityMatrix<SCFMode> pDens(basisController);
    for_spin(pDens, PAO) {
      pDens_spin = PAO_spin + PAO_spin.transpose();
    };
    std::shared_ptr<DensityMatrixController<SCFMode>> pDmatController(std::make_shared<DensityMatrixController<SCFMode>>(pDens));

    auto pDensOnGridController =
        std::make_shared<DensityMatrixDensityOnGridController<SCFMode>>(densityOnGridCalculator, pDmatController);
    auto& pDensOnGrid = pDensOnGridController->getDensityOnGrid();
    auto& pDensGradientOnGrid = pDensOnGridController->getDensityGradientOnGrid();

    vxcGradient = evaluateVxcGradientContribution(pDens);
    OutputControl::vOut << "  Vxc contribution\n" << vxcGradient << std::endl;
    vxcKernelGradient = evaluateFxcGradientContribution(D, pDensOnGrid, pDensGradientOnGrid);
    OutputControl::vOut << "  Kernel contribution from vxc\n" << vxcKernelGradient << std::endl;
    fxcGradient = (SCFMode == Options::SCF_MODES::RESTRICTED ? 2.0 : 1.0) *
                  evaluateFxcGradientContribution(xDens, xDensOnGrid, xDensGradientOnGrid);
    OutputControl::vOut << "  Kernel contribution\n" << fxcGradient << std::endl;
    hyperkernelGradient = (SCFMode == Options::SCF_MODES::RESTRICTED ? 2.0 : 1.0) *
                          evaluateKxcGradientContribution(D, xDensOnGrid, xDensGradientOnGrid);
    OutputControl::vOut << "  Hyperkernel contribution \n" << hyperkernelGradient << std::endl;

    Timings::timeTaken("LRSCF Gradients -        XC Contribution");
  } /* end _usesXC */

  Eigen::MatrixXd twoCenterGradient = Eigen::MatrixXd::Zero(nAtoms, 3);

  if (_lrscf->getLRSCFSettings().densFitJ == Options::DENS_FITS::NONE) {
    twoCenterGradient += evaluateFullTwoCenterGradient(PAO, D, xpyAO, xmyAO);
  }
  // densFitJ case, but full exchange
  else {
    std::shared_ptr<AtomCenteredBasisController> auxBasisController =
        _lrscf->getSys()->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);

    RIIntegralDerivativeCalculator riintderivs(basisController, auxBasisController);
    Eigen::MatrixXd riCoulombContr = riintderivs.getCoulombIntegralDerivatives(
        xpyAO.total() * (SCFMode == Options::SCF_MODES::RESTRICTED ? std::sqrt(2) : 1.0), D.total(), PAO.total());
    OutputControl::vOut << "  excited-state Coulomb contribution (RI):\n" << riCoulombContr << std::endl;

    twoCenterGradient += evaluateFullExchangeGradient(PAO, D, xpyAO, xmyAO);
    OutputControl::dOut << "  excited-state exchange contribution:\n" << twoCenterGradient << std::endl;
    twoCenterGradient += riCoulombContr;
  }

  if (_func.getLRExchangeRatio()) {
    Eigen::MatrixXd coulombGradientLR = evaluateFullLRExchangeGradient(PAO, D, xpyAO, xmyAO);
    OutputControl::vOut << "  excited-state long-range exchange contribution:\n" << coulombGradientLR << std::endl;
    twoCenterGradient += coulombGradientLR;
  }
  OutputControl::vOut << "  entire excited-state two-center contribution:\n" << twoCenterGradient << std::endl;

  Eigen::MatrixXd grads = twoCenterGradient + overlapGradient + hCoreGradient + _groundstateGradient + fxcGradient +
                          vxcGradient + vxcKernelGradient + hyperkernelGradient;

  const auto geom = _lrscf->getSys()->getGeometry();
  const auto& atoms = geom->getAtoms();
  OutputControl::nOut << "  ground-state gradient:" << std::endl;
  for (unsigned iAtom = 0; iAtom < nAtoms; iAtom++) {
    OutputControl::n.printf("%4s %4d %2s %+15.10f %+15.10f %+15.10f\n", "", iAtom + 1,
                            atoms[iAtom]->getAtomType()->getElementSymbol().c_str(), _groundstateGradient(iAtom, 0),
                            _groundstateGradient(iAtom, 1), _groundstateGradient(iAtom, 2));
  }
  OutputControl::nOut << std::endl;
  geom->setGradients(grads);
  OutputControl::m.printf("  Excitation energy: %.6f a.u.   State energy: %.6f a.u.\n\n", omega,
                          omega + _lrscf->getSys()->template getElectronicStructure<SCFMode>()->getEnergy());
  geom->printGradients();
  OutputControl::mOut << std::endl;

  Eigen::VectorXd totalForce = grads.colwise().sum();
  OutputControl::n.printf("  Difference to translational invariance:\n %11s %+15.10f %15.10f %15.10f\n\n", "",
                          totalForce(0), totalForce(1), totalForce(2));
  const Eigen::MatrixXd& coords = geom->getCoordinates();
  Eigen::Vector3d torque = Eigen::Vector3d::Zero(3);
  for (unsigned iAtom = 0; iAtom < nAtoms; iAtom++) {
    Eigen::Vector3d gradv = grads.row(iAtom);
    Eigen::Vector3d coordv = coords.row(iAtom);
    torque += gradv.cross(coordv);
  }
  OutputControl::n.printf("  Difference to rotational invariance:\n %11s %+15.10f %15.10f %15.10f\n\n", "", torque(0),
                          torque(1), torque(2));

  return grads;
}

// first functional derivative
template<Options::SCF_MODES SCFMode>
const Eigen::MatrixXd TDDFTGradientCalculator<SCFMode>::evaluateVxcGradientContribution(const DensityMatrix<SCFMode>& densitymat) {
  // variable for the resulting gradient contribution
  Eigen::MatrixXd gradient(Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3));
  // Get potential
  const auto p = _funcData->dFdRho;
  const auto pg = _funcData->dFdGradRho;
  // GGA case
  if (_kernel->isGGA()) {
    _gridToMatrix->addScalarOperatorToGradient(
        gradient, _lrscf->getSys()->getAtomCenteredBasisController()->getAtomBasisProjection(), densitymat, *p, *pg);
  }
  // LDA
  else {
    _gridToMatrix->addScalarOperatorToGradient(
        gradient, _lrscf->getSys()->getAtomCenteredBasisController()->getAtomBasisProjection(), densitymat, *p);
  }
  return gradient;
}

// LDA
template<>
void TDDFTGradientCalculator<RESTRICTED>::prepareFxcGradient(GridData<RESTRICTED>& scalar,
                                                             const std::shared_ptr<d2F_dRho2<RESTRICTED>> pp,
                                                             const GridData<RESTRICTED>& dens) {
  scalar.array() = dens.array() * pp->array();
}
// LDA
template<>
void TDDFTGradientCalculator<UNRESTRICTED>::prepareFxcGradient(GridData<UNRESTRICTED>& scalar,
                                                               const std::shared_ptr<d2F_dRho2<UNRESTRICTED>> pp,
                                                               const GridData<UNRESTRICTED>& dens) {
  scalar.alpha.array() = dens.alpha.array() * pp->aa.array() + dens.beta.array() * pp->ab.array();
  scalar.beta.array() = dens.alpha.array() * pp->ab.array() + dens.beta.array() * pp->bb.array();
}

// GGA
template<>
void TDDFTGradientCalculator<RESTRICTED>::prepareFxcGradient(
    GridData<RESTRICTED>& scalar, Gradient<GridData<RESTRICTED>>& gradient, const std::shared_ptr<d2F_dRho2<RESTRICTED>> pp,
    const std::shared_ptr<Gradient<DoublySpinPolarizedData<RESTRICTED, GridData<RESTRICTED>>>> pg,
    const std::shared_ptr<Hessian<DoublySpinPolarizedData<RESTRICTED, GridData<RESTRICTED>>>> gg,
    const GridData<RESTRICTED>& dens, const Gradient<GridData<RESTRICTED>>& densg) {
  scalar.array() = dens.array() * pp->array() + densg.x.array() * pg->x.array() + densg.y.array() * pg->y.array() +
                   densg.z.array() * pg->z.array();
  gradient.x.array() = dens.array() * pg->x.array() + densg.x.array() * gg->xx.array() +
                       densg.y.array() * gg->xy.array() + densg.z.array() * gg->xz.array();
  gradient.y.array() = dens.array() * pg->y.array() + densg.x.array() * gg->xy.array() +
                       densg.y.array() * gg->yy.array() + densg.z.array() * gg->yz.array();
  gradient.z.array() = dens.array() * pg->z.array() + densg.x.array() * gg->xz.array() +
                       densg.y.array() * gg->yz.array() + densg.z.array() * gg->zz.array();
}

// GGA
template<>
void TDDFTGradientCalculator<UNRESTRICTED>::prepareFxcGradient(
    GridData<UNRESTRICTED>& scalar, Gradient<GridData<UNRESTRICTED>>& gradient,
    const std::shared_ptr<d2F_dRho2<UNRESTRICTED>> pp,
    const std::shared_ptr<Gradient<DoublySpinPolarizedData<UNRESTRICTED, GridData<RESTRICTED>>>> pg,
    const std::shared_ptr<Hessian<DoublySpinPolarizedData<UNRESTRICTED, GridData<RESTRICTED>>>> gg,
    const GridData<UNRESTRICTED>& dens, const Gradient<GridData<UNRESTRICTED>>& densg) {
  scalar.alpha.array() = dens.alpha.array() * pp->aa.array() + dens.beta.array() * pp->ab.array() +
                         densg.x.alpha.array() * pg->x.aa.array() + densg.x.beta.array() * pg->x.ab.array() +
                         densg.y.alpha.array() * pg->y.aa.array() + densg.y.beta.array() * pg->y.ab.array() +
                         densg.z.alpha.array() * pg->z.aa.array() + densg.z.beta.array() * pg->z.ab.array();
  scalar.beta.array() = dens.alpha.array() * pp->ab.array() + dens.beta.array() * pp->bb.array() +
                        densg.x.alpha.array() * pg->x.ba.array() + densg.x.beta.array() * pg->x.bb.array() +
                        densg.y.alpha.array() * pg->y.ba.array() + densg.y.beta.array() * pg->y.bb.array() +
                        densg.z.alpha.array() * pg->z.ba.array() + densg.z.beta.array() * pg->z.bb.array();
  gradient.x.alpha.array() = dens.alpha.array() * pg->x.aa.array() + dens.beta.array() * pg->x.ba.array() +
                             densg.x.alpha.array() * gg->xx.aa.array() + densg.x.beta.array() * gg->xx.ab.array() +
                             densg.y.alpha.array() * gg->xy.aa.array() + densg.y.beta.array() * gg->xy.ab.array() +
                             densg.z.alpha.array() * gg->xz.aa.array() + densg.z.beta.array() * gg->xz.ab.array();
  gradient.x.beta.array() = dens.alpha.array() * pg->x.ab.array() + dens.beta.array() * pg->x.bb.array() +
                            densg.x.alpha.array() * gg->xx.ab.array() + densg.x.beta.array() * gg->xx.bb.array() +
                            densg.y.alpha.array() * gg->xy.ba.array() + densg.y.beta.array() * gg->xy.bb.array() +
                            densg.z.alpha.array() * gg->xz.ba.array() + densg.z.beta.array() * gg->xz.bb.array();
  gradient.y.alpha.array() = dens.alpha.array() * pg->y.aa.array() + dens.beta.array() * pg->y.ba.array() +
                             densg.x.alpha.array() * gg->xy.aa.array() + densg.x.beta.array() * gg->xy.ba.array() +
                             densg.y.alpha.array() * gg->yy.aa.array() + densg.y.beta.array() * gg->yy.ab.array() +
                             densg.z.alpha.array() * gg->yz.aa.array() + densg.z.beta.array() * gg->yz.ab.array();
  gradient.y.beta.array() = dens.alpha.array() * pg->y.ab.array() + dens.beta.array() * pg->y.bb.array() +
                            densg.x.alpha.array() * gg->xy.ab.array() + densg.x.beta.array() * gg->xy.bb.array() +
                            densg.y.alpha.array() * gg->yy.ba.array() + densg.y.beta.array() * gg->yy.bb.array() +
                            densg.z.alpha.array() * gg->yz.ba.array() + densg.z.beta.array() * gg->yz.bb.array();
  gradient.z.alpha.array() = dens.alpha.array() * pg->z.aa.array() + dens.beta.array() * pg->z.ba.array() +
                             densg.x.alpha.array() * gg->xz.aa.array() + densg.x.beta.array() * gg->xz.ba.array() +
                             densg.y.alpha.array() * gg->yz.aa.array() + densg.y.beta.array() * gg->yz.ba.array() +
                             densg.z.alpha.array() * gg->zz.aa.array() + densg.z.beta.array() * gg->zz.ab.array();
  gradient.z.beta.array() = dens.alpha.array() * pg->z.ab.array() + dens.beta.array() * pg->z.bb.array() +
                            densg.x.alpha.array() * gg->xz.ab.array() + densg.x.beta.array() * gg->xz.bb.array() +
                            densg.y.alpha.array() * gg->yz.ab.array() + densg.y.beta.array() * gg->yz.bb.array() +
                            densg.z.alpha.array() * gg->zz.ab.array() + densg.z.beta.array() * gg->zz.bb.array();
}

// second functional derivative
template<Options::SCF_MODES SCFMode>
const Eigen::MatrixXd
TDDFTGradientCalculator<SCFMode>::evaluateFxcGradientContribution(const DensityMatrix<SCFMode>& densitymat,
                                                                  const GridData<SCFMode>& densityOnGrid,
                                                                  const Gradient<GridData<SCFMode>>& densityOnGridGradient) {
  std::shared_ptr<AtomCenteredGridController> gridController = AtomCenteredGridControllerFactory::produce(
      GeometryAdderFactory::produce({_lrscf->getSys()}, _lrscf->getLRSCFSettings().subsystemgrid),
      _lrscf->getLRSCFSettings().grid, Options::GRID_PURPOSES::SMALL);

  Eigen::MatrixXd gradient(Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3));

  const auto pp = _funcData->d2FdRho2;

  GridData<SCFMode> scalar(gridController);
  // GGA
  if (_kernel->isGGA()) {
    Gradient<GridData<SCFMode>> gradientpart = makeGradient<GridData<SCFMode>>(gridController);
    const auto pg = _funcData->d2FdRhodGradRho;
    const auto gg = _funcData->d2FdGradRho2;
    prepareFxcGradient(scalar, gradientpart, pp, pg, gg, densityOnGrid, densityOnGridGradient);
    _gridToMatrix->addScalarOperatorToGradient(gradient,
                                               _lrscf->getSys()->getAtomCenteredBasisController()->getAtomBasisProjection(),
                                               densitymat, scalar, gradientpart);
  }
  // LDA
  else {
    prepareFxcGradient(scalar, pp, densityOnGrid);
    _gridToMatrix->addScalarOperatorToGradient(
        gradient, _lrscf->getSys()->getAtomCenteredBasisController()->getAtomBasisProjection(), densitymat, scalar);
  }

  return gradient;
}

// third functional derivative
template<Options::SCF_MODES SCFMode>
const Eigen::MatrixXd
TDDFTGradientCalculator<SCFMode>::evaluateKxcGradientContribution(const DensityMatrix<SCFMode>& densitymat,
                                                                  const GridData<SCFMode>& densityOnGrid,
                                                                  const Gradient<GridData<SCFMode>>& densityOnGridGradient) {
  std::shared_ptr<AtomCenteredGridController> gridController = AtomCenteredGridControllerFactory::produce(
      GeometryAdderFactory::produce({_lrscf->getSys()}, _lrscf->getLRSCFSettings().subsystemgrid),
      _lrscf->getLRSCFSettings().grid, Options::GRID_PURPOSES::SMALL);
  Eigen::MatrixXd gradient(Eigen::MatrixXd::Zero(_lrscf->getSys()->getNAtoms(), 3));

  GridData<SCFMode> scalar(gridController);
  Gradient<GridData<SCFMode>> gradientPart(makeGradient<GridData<SCFMode>>(gridController));
  // contractKernel does what the "prepareFxcGradient" functions do in the second functional derivative case
  _hyperSig->contractKernel(scalar, gradientPart, densityOnGrid, densityOnGridGradient);

  // GGA
  if (_kernel->isGGA()) {
    _gridToMatrix->addScalarOperatorToGradient(gradient,
                                               _lrscf->getSys()->getAtomCenteredBasisController()->getAtomBasisProjection(),
                                               densitymat, scalar, gradientPart);
  }
  // LDA
  else {
    _gridToMatrix->addScalarOperatorToGradient(
        gradient, _lrscf->getSys()->getAtomCenteredBasisController()->getAtomBasisProjection(), densitymat, scalar);
  }
  return gradient;
}

template<Options::SCF_MODES SCFMode>
std::vector<Eigen::MatrixXd> TDDFTGradientCalculator<SCFMode>::calculateGradients() {
  std::shared_ptr<AtomCenteredBasisController> basisController = _lrscf->getSys()->getAtomCenteredBasisController();
  unsigned nBasis = basisController->getNBasisFunctions();

  CoefficientMatrix<SCFMode> C = _lrscf->getCoefficients();

  auto nOcc = _lrscf->getSys()->template getNOccupiedOrbitals<SCFMode>();
  auto nVirt = _lrscf->getSys()->template getNVirtualOrbitals<SCFMode>();

  const unsigned nExcitations = _excGradList.size();
  const unsigned int nAtoms = _lrscf->getSys()->getNAtoms();
  std::vector<Eigen::MatrixXd> grads(nExcitations, Eigen::MatrixXd::Zero(nAtoms, 3));
  Options::LRSCF_TYPE type = Options::LRSCF_TYPE::ISOLATED;
  if (_lrscf->getEnvSystems().size() > 0) {
    type = Options::LRSCF_TYPE::UNCOUPLED;
  }

  // nDim ist nOcc.alpha * nVirt.alpha + nOcc.beta * nVirt.beta
  unsigned nDim = (*(_lrscf->getExcitationVectors(type)))[0].rows();
  Eigen::MatrixXd rightHandSides = Eigen::MatrixXd::Zero(nDim, nExcitations);
  std::vector<MatrixInBasis<SCFMode>> xpyAOs(nExcitations, MatrixInBasis<SCFMode>(basisController));
  std::vector<MatrixInBasis<SCFMode>> xmyAOs(nExcitations, MatrixInBasis<SCFMode>(basisController));
  std::vector<MatrixInBasis<SCFMode>> Ts(nExcitations, MatrixInBasis<SCFMode>(basisController));
  std::vector<MatrixInBasis<SCFMode>> Ws(nExcitations, MatrixInBasis<SCFMode>(basisController));

  DensityMatrix<SCFMode> D = _lrscf->getSys()->template getElectronicStructure<SCFMode>()->getDensityMatrix();
  const SPMatrix<SCFMode> fock = *(_lrscf->getMOFockMatrix());

  Timings::takeTime("LRSCF Gradients -       Matrix Formation");
  Timings::takeTime("LRSCF Gradients - Z-Vector RHS Formation");
  for (unsigned iExc = 0; iExc < nExcitations; iExc++) {
    MatrixInBasis<SCFMode>& T = Ts[iExc];
    MatrixInBasis<SCFMode> xpy(basisController);
    MatrixInBasis<SCFMode> xmy(basisController);

    unsigned iaStart = 0;
    for_spin(T, xpy, xmy, nOcc, nVirt) {
      xpy_spin = Eigen::MatrixXd::Zero(nBasis, nBasis);
      xmy_spin = Eigen::MatrixXd::Zero(nBasis, nBasis);
      xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) = Eigen::Map<Eigen::MatrixXd>(
          (*(_lrscf->getExcitationVectors(type)))[0].col(_excGradList[iExc] - 1).data() + iaStart, nVirt_spin, nOcc_spin);
      xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) = Eigen::Map<Eigen::MatrixXd>(
          (*(_lrscf->getExcitationVectors(type)))[1].col(_excGradList[iExc] - 1).data() + iaStart, nVirt_spin, nOcc_spin);
      T_spin.bottomRightCorner(nVirt_spin, nVirt_spin) =
          0.5 *
          (xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() +
           xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose());
      T_spin.topLeftCorner(nOcc_spin, nOcc_spin) =
          -0.5 *
          (xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) +
           xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin));
      iaStart += nOcc_spin * nVirt_spin;
    };

    OutputControl::dOut << "  XpYtot\n" << xpy.total() << "\n  XmYtot\n" << xmy.total() << std::endl;
    MatrixInBasis<SCFMode> hPX = hPlus(xpy);
    OutputControl::dOut << "  HpX tot\n" << hPX.total() << std::endl;
    MatrixInBasis<SCFMode> hMX = hMinus(xmy);
    OutputControl::dOut << "  HmX tot\n" << hMX.total() << "\n  T tot\n" << T.total() << std::endl;
    MatrixInBasis<SCFMode> hPT = hPlus(T);

    MatrixInBasis<SCFMode> rhs(basisController);
    for_spin(rhs, xpy, xmy, hPT, hMX, hPX, nOcc, nVirt) {
      rhs_spin = hPX_spin.bottomRightCorner(nVirt_spin, nVirt_spin) * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) +
                 hMX_spin.bottomRightCorner(nVirt_spin, nVirt_spin) * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) -
                 xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * hPX_spin.topLeftCorner(nOcc_spin, nOcc_spin) -
                 xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * hMX_spin.topLeftCorner(nOcc_spin, nOcc_spin) +
                 hPT_spin.topRightCorner(nOcc_spin, nVirt_spin).transpose();
    };
    MatrixInBasis<SCFMode>& xpyAO = xpyAOs[iExc];
    MatrixInBasis<SCFMode>& xmyAO = xmyAOs[iExc];
    MatrixInBasis<SCFMode>& W = Ws[iExc];
    for_spin(C, xpy, xmy, xpyAO, xmyAO) {
      xpyAO_spin = C_spin * xpy_spin * (Eigen::MatrixXd)C_spin.transpose();
      xmyAO_spin = C_spin * xmy_spin * (Eigen::MatrixXd)C_spin.transpose();
    };

    double omega = (*(_lrscf->getExcitationEnergies(type)))(_excGradList[iExc] - 1);

    for_spin(W, nOcc, nVirt, fock, xpy, xmy, hPX, hMX) {
      W_spin.topLeftCorner(nOcc_spin, nOcc_spin) +=
          0.5 * omega *
              (xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) +
               xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin)) -
          0.5 * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() *
              fock_spin.bottomRightCorner(nVirt_spin, nVirt_spin) * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) -
          0.5 * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() *
              fock_spin.bottomRightCorner(nVirt_spin, nVirt_spin) * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin);
      W_spin.topRightCorner(nOcc_spin, nVirt_spin) = hPX_spin.topLeftCorner(nOcc_spin, nOcc_spin).transpose() *
                                                         xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() +
                                                     hMX_spin.topLeftCorner(nOcc_spin, nOcc_spin).transpose() *
                                                         xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose();
      W_spin.bottomRightCorner(nVirt_spin, nVirt_spin) =
          0.5 * omega *
              (xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() +
               xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose()) +
          0.5 * xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * fock_spin.topLeftCorner(nOcc_spin, nOcc_spin) *
              xpy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose() +
          0.5 * xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin) * fock_spin.topLeftCorner(nOcc_spin, nOcc_spin) *
              xmy_spin.bottomLeftCorner(nVirt_spin, nOcc_spin).transpose();
    };
    if (_usesXC) {
      Timings::takeTime("LRSCF Gradients -        XC Contribution");
      MatrixInBasis<SCFMode> xDens(basisController);
      for_spin(xDens, xpyAO) {
        xDens_spin = xpyAO_spin;
        xDens_spin += xpyAO_spin.transpose();
      };
      MatrixInBasis<SCFMode> hyperKernelContr = *(_hyperSig->calcF(std::make_shared<MatrixInBasis<SCFMode>>(xDens)));
      double scfFactor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;
      for_spin(C, nOcc, nVirt, hyperKernelContr, rhs, W) {
        OutputControl::dOut << "  Hyperkernelcontribution spin:\n" << hyperKernelContr_spin << std::endl;
        rhs_spin += scfFactor * scfFactor * C_spin.middleCols(nOcc_spin, nVirt_spin).transpose() *
                    hyperKernelContr_spin * C_spin.leftCols(nOcc_spin);
        W_spin.topLeftCorner(nOcc_spin, nOcc_spin) += 0.5 * scfFactor * scfFactor * C_spin.leftCols(nOcc_spin).transpose() *
                                                      hyperKernelContr_spin * C_spin.leftCols(nOcc_spin);
        OutputControl::dOut << "  full RHS spin\n" << rhs_spin << std::endl;
      };
      OutputControl::dOut << "  Hyperkernelcontribution total:\n" << hyperKernelContr.total() << std::endl;
      OutputControl::dOut << "  full RHS total\n" << rhs.total() << std::endl;
      Timings::timeTaken("LRSCF Gradients -        XC Contribution");
    }

    unsigned ia = 0;
    for_spin(nOcc, nVirt, rhs) {
      for (unsigned i = 0; i < nOcc_spin; i++) {
        for (unsigned a = 0; a < nVirt_spin; a++) {
          // ia = i * nVirt_spin + a + iaStart
          rightHandSides(ia++, iExc) = -rhs_spin(a, i);
        }
      }
    };
  } // for iExc
  Timings::timeTaken("LRSCF Gradients - Z-Vector RHS Formation");

  // here, I choose one zero frequency and as many RHS as gradients required (but still only one set)
  std::vector<double> freq(1, 0.0);
  std::vector<Eigen::MatrixXd> myRHS = {rightHandSides};
  const auto& lrscfSettings = _lrscf->getLRSCFSettings();
  Eigen::VectorXd diagonal = LRSCFSetup<SCFMode>::getDiagonal({_lrscf});
  // use the NonlinearResponseSolver because the ResponseSolver only respects the first three RHS columns
  std::unique_ptr<NonlinearResponseSolver> responseSolver;
  responseSolver = std::make_unique<NonlinearResponseSolver>(diagonal, lrscfSettings.conv, lrscfSettings.maxCycles,
                                                             lrscfSettings.maxSubspaceDimension, freq, myRHS, _nlSolver);
  Timings::takeTime("LRSCF Gradients - Solve the Z-Vector Eq.");
  responseSolver->solve();
  Eigen::MatrixXd Z = responseSolver->getEigenvectors()[0];
  Timings::timeTaken("LRSCF Gradients - Solve the Z-Vector Eq.");
  Timings::timeTaken("LRSCF Gradients -       Matrix Formation");
  for (unsigned iExc = 0; iExc < nExcitations; iExc++) {
    Timings::takeTime("LRSCF Gradients -       Matrix Formation");
    const MatrixInBasis<SCFMode>& T = Ts[iExc];
    MatrixInBasis<SCFMode>& W = Ws[iExc];
    const MatrixInBasis<SCFMode>& xpyAO = xpyAOs[iExc];
    const MatrixInBasis<SCFMode>& xmyAO = xmyAOs[iExc];
    MatrixInBasis<SCFMode> P = T;
    unsigned iaStart = 0;
    for_spin(P, nOcc, nVirt) {
      P_spin.topRightCorner(nOcc_spin, nVirt_spin) +=
          Eigen::Map<Eigen::MatrixXd>(Z.col(iExc).data() + iaStart, nVirt_spin, nOcc_spin).transpose();
      iaStart += nOcc_spin * nVirt_spin;
    };

    _lrscf->getRelaxedDiffDensities()->push_back(P);

    OutputControl::dOut << "  P tot =\n" << P.total() << std::endl;
    MatrixInBasis<SCFMode> hPlusOfP = hPlus(P);
    OutputControl::dOut << "  HpP\n" << hPlusOfP.total() << std::endl;
    for_spin(P) {
      OutputControl::dOut << "  P spin\n" << P_spin << std::endl;
    };
    for_spin(hPlusOfP) {
      OutputControl::dOut << "  HpP spin\n" << hPlusOfP_spin << std::endl;
    };

    double omega = (*(_lrscf->getExcitationEnergies(type)))(_excGradList[iExc] - 1);

    for_spin(W, nOcc, nVirt, P, hPlusOfP, fock) {
      W_spin.topLeftCorner(nOcc_spin, nOcc_spin) += 0.5 * hPlusOfP_spin.topLeftCorner(nOcc_spin, nOcc_spin);
      // topRightCorner of P is the Z vector
      W_spin.topRightCorner(nOcc_spin, nVirt_spin) +=
          fock_spin.topLeftCorner(nOcc_spin, nOcc_spin) * P_spin.topRightCorner(nOcc_spin, nVirt_spin);
    };

    MatrixInBasis<SCFMode> PAO(basisController), WAO(basisController);
    for_spin(PAO, P, C, WAO, W) {
      PAO_spin = C_spin * P_spin * (Eigen::MatrixXd)C_spin.transpose();
      WAO_spin = C_spin * W_spin * (Eigen::MatrixXd)C_spin.transpose();
      OutputControl::dOut << "  W spin\n" << W_spin << std::endl;
    };
    OutputControl::dOut << "  W tot\n"
                        << W.total() << "\n  D tot\n"
                        << D.total() << "\n  PAO tot\n"
                        << PAO.total() << "\n  WAO tot\n"
                        << WAO.total() << "\n  xpyAO\n"
                        << xpyAO << "\n  xmyAO\n"
                        << xmyAO << std::endl;

    Timings::timeTaken("LRSCF Gradients -       Matrix Formation");
    Timings::takeTime("LRSCF Gradients -    Gradient Evaluation");
    grads[iExc] = this->evaluateGradients(PAO, WAO, D, xpyAO, xmyAO, omega);
    Timings::timeTaken("LRSCF Gradients -    Gradient Evaluation");

  } /* Loop over excitations */

  return grads;
}

template class TDDFTGradientCalculator<Options::SCF_MODES::RESTRICTED>;
template class TDDFTGradientCalculator<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */