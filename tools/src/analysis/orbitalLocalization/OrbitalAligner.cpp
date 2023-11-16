/**
 * @file OrbitalAligner.cpp
 *
 * @date Mar 14, 2019
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
#include "analysis/orbitalLocalization/OrbitalAligner.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/IAOPopulationCalculator.h" //IAO populations for rotating.
#include "basis/Basis.h"                                         //get N contracted.
#include "data/ElectronicStructure.h"                            //Access to the density matrix.
#include "data/OrbitalController.h"                              //Access to occupied orbital coefficients.
#include "data/matrices/CoefficientMatrix.h"                     //Coefficient matrix definition.
#include "geometry/Geometry.h"                                   //Check atom ordering.
#include "integrals/OneElectronIntegralController.h"             //Kinetic integrals.
#include "io/FormattedOutputStream.h"                            //Filtered output streams.
#include "math/linearAlgebra/JacobiRotation.h"                   //2x2 Rotations.
#include "math/optimizer/NewtonRaphson.h"                        //Numerical minimization.
#include "system/SystemController.h"                             //System controller definition.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
OrbitalAligner<SCFMode>::OrbitalAligner(std::shared_ptr<SystemController> systemController,
                                        std::shared_ptr<SystemController> templateSystem, unsigned int exponent,
                                        bool kineticAlign, bool replaceVirtuals)
  : _system(systemController),
    _templateSystem(templateSystem),
    _exponent(exponent),
    _kineticAlign(kineticAlign),
    _replaceVirtuals(replaceVirtuals) {
  // Some sanity checks.
  assert(_system);
  assert(_templateSystem);
  if (_system->getGeometry()->getNAtoms() != _templateSystem->getGeometry()->getNAtoms())
    throw SerenityError("Template and actual system have to contain the same atoms.");
  auto nOccSys = _system->getNOccupiedOrbitals<SCFMode>();
  auto nOccTem = _templateSystem->getNOccupiedOrbitals<SCFMode>();
  for_spin(nOccSys, nOccTem) {
    if (nOccSys_spin < nOccTem_spin)
      throw SerenityError("The occupations of the system and the template do not fit!");
  };
  auto atomsSys = _system->getGeometry()->getAtoms();
  auto atomsTem = _templateSystem->getGeometry()->getAtoms();
  for (unsigned int iAtom = 0; iAtom < atomsSys.size(); ++iAtom) {
    if (atomsSys[iAtom]->getAtomType()->getElementSymbol() != atomsTem[iAtom]->getAtomType()->getElementSymbol())
      throw SerenityError("The atoms of template and actual system have to be ordered in the same way!");
  }
  if (_exponent % 2 != 0)
    throw SerenityError("The exponent has to be even! You used " + std::to_string(_exponent));
}

template<Options::SCF_MODES SCFMode>
void OrbitalAligner<SCFMode>::localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                               SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) {
  SpinPolarizedData<SCFMode, Eigen::VectorXi> optSys;
  SpinPolarizedData<SCFMode, Eigen::VectorXi> optTem;

  auto nOrbitals = _system->getNOccupiedOrbitals<SCFMode>();
  _localizeVirtuals = false;
  for_spin(orbitalRange, nOrbitals) {
    if (orbitalRange_spin.size() > 1) {
      unsigned int maxIndex = *std::max_element(orbitalRange_spin.begin(), orbitalRange_spin.end());
      if (maxIndex >= nOrbitals_spin) {
        _localizeVirtuals = true;
        nOrbitals_spin = maxIndex + 1;
      }
    }
  };

  for_spin(optSys, optTem, orbitalRange, nOrbitals) {
    optSys_spin = Eigen::VectorXi::Zero(nOrbitals_spin);
    optTem_spin = Eigen::VectorXi::Zero(nOrbitals_spin);
    for (const auto orb : orbitalRange_spin) {
      optSys_spin(orb) = 1;
      optTem_spin(orb) = 1;
    }
  };
  alignOrbitals(orbitals, optSys, optTem, orbitalRange, maxSweeps);
}
template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd> OrbitalAligner<SCFMode>::getReferenceKineticEnergies() {
  Eigen::setNbThreads(1);
  const auto& coeff = _templateSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  auto kinIntegrals = _templateSystem->getOneElectronIntegralController()->getKinIntegrals();
  auto nOcc = _system->getNOccupiedOrbitals<SCFMode>();
  SpinPolarizedData<SCFMode, Eigen::VectorXd> kineticEnergies;
  for_spin(nOcc, coeff, kineticEnergies) {
    kineticEnergies_spin.resize(nOcc_spin);
    for (unsigned int occI = 0; occI < nOcc_spin; ++occI) {
      const auto& coeffI = coeff_spin.col(occI);
      kineticEnergies_spin[occI] = coeffI.transpose() * kinIntegrals * coeffI;
    } // for occI
  };
  Eigen::setNbThreads(0);
  // save the kinetic energies.
  return kineticEnergies;
}

template<Options::SCF_MODES SCFMode>
void OrbitalAligner<SCFMode>::addToLagrangian(const double Qii, const double Qjj, const double Qij, const double qi,
                                              const double qj, const double x, const int jRef, double& value,
                                              double& gradient, double& hessian) {
  const double P = (double)_exponent;
  double cosX = std::cos(x);
  double sinX = std::sin(x);
  double cosX2 = cosX * cosX;
  double sinX2 = sinX * sinX;
  double sinCosX = sinX * cosX;
  double iDiff = (qi - Qii * cosX2 - 2.0 * Qij * sinCosX - Qjj * sinX2);
  double jDiff = (jRef < 0) ? 0.0 : (qj - Qjj * cosX2 + 2.0 * Qij * sinCosX - Qii * sinX2);
  // Update value
  double iDiffPowE2 = 1;
  double jDiffPowE2 = 1;
  for (unsigned int p = 0; p < _exponent - 2; ++p) {
    iDiffPowE2 *= iDiff;
    jDiffPowE2 *= jDiff;
  }
  const double iDiffPowE1 = iDiffPowE2 * iDiff;
  const double jDiffPowE1 = jDiffPowE2 * jDiff;
  const double iDiffPowE = iDiffPowE1 * iDiff;
  const double jDiffPowE = jDiffPowE1 * jDiff;
  value += iDiffPowE + jDiffPowE;
  // Update gradient
  const double cosX2_sinX2 = cosX2 - sinX2;
  double Qii_Qjj = Qii - Qjj;
  const double gradIDiff = 2.0 * Qii_Qjj * sinCosX - 2.0 * Qij * cosX2_sinX2;
  gradient += P * gradIDiff * (iDiffPowE1 - jDiffPowE1);
  // Update hessian
  const double gradSquare = gradIDiff * gradIDiff;
  const double nablaIDiff = 2.0 * Qii_Qjj * cosX2_sinX2 + 4.0 * sinCosX * Qij;
  const double tmp1 = (P - 1) * (iDiffPowE2 + jDiffPowE2) * gradSquare;
  const double tmp2 = (iDiffPowE1 - jDiffPowE1) * nablaIDiff;
  hessian += P * (tmp1 + tmp2); // 2.0 * (hI+hJ);
}

template<Options::SCF_MODES SCFMode>
void OrbitalAligner<SCFMode>::alignOrbitals(OrbitalController<SCFMode>& orbitals,
                                            const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitals,
                                            const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitalsRef,
                                            unsigned int maxSweeps, double angleThreshold) {
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange;
  for_spin(orbitalRange, unpairedOrbitals) {
    unsigned int nOcc = unpairedOrbitals_spin.rows();
    for (unsigned int i = 0; i < nOcc; ++i)
      orbitalRange_spin.push_back(i);
  };

  alignOrbitals(orbitals, unpairedOrbitals, unpairedOrbitalsRef, orbitalRange, maxSweeps, angleThreshold);
}

template<Options::SCF_MODES SCFMode>
void OrbitalAligner<SCFMode>::alignOrbitals(OrbitalController<SCFMode>& orbitals,
                                            const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitals,
                                            const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitalsRef,
                                            SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange,
                                            unsigned int maxSweeps, double angleThreshold) {
  OutputControl::mOut << "Aligning the orbitals of " << _system->getSystemName() << " to "
                      << _templateSystem->getSystemName() << std::endl;

  CoefficientMatrix<SCFMode> C = orbitals.getCoefficients();
  auto densMatCont = _system->getElectronicStructure<SCFMode>()->getDensityMatrixController();
  densMatCont->updateDensityMatrix();
  DensityMatrix<SCFMode> oldP = densMatCont->getDensityMatrix();
  /*
   * The following few lines are needed for the localization of virtual orbitals only. In the
   * case that the IAO basis set does not span the space of the virtual valence orbitals we must
   * first reconstruct the orbitals to do so. Otherwise, the orbital localization may lead to
   * non-orthogonal and non normalized virtual orbitals. In the case that "_replaceVirtuals" is
   * false, we explicitly check if the virtual valence orbital space is fine and replace the
   * virtual orbitals if required.
   */
  bool replaceVirtBeforeLoc =
      _localizeVirtuals &&
      (_replaceVirtuals || !IAOPopulationCalculator<SCFMode>::iaosSpanOrbitals(
                               C, _system->getOneElectronIntegralController()->getOverlapIntegrals(),
                               _system->template getNOccupiedOrbitals<SCFMode>(), _system->getBasisController(),
                               _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION)));
  if (replaceVirtBeforeLoc) {
    IAOPopulationCalculator<SCFMode>::reconstructVirtualValenceOrbitalsInplace(
        C, _system->getOneElectronIntegralController()->getOverlapIntegrals(),
        _system->template getNOccupiedOrbitals<SCFMode>(), _system->getBasisController(),
        _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION));
    orbitals.updateOrbitals(C, orbitals.getEigenvalues());
  }
  auto CIAO_othoA = IAOPopulationCalculator<SCFMode>::getCIAOCoefficients(_system, _localizeVirtuals);
  auto pops_ref = IAOPopulationCalculator<SCFMode>::calculateShellwiseOrbitalPopulations(_templateSystem, _localizeVirtuals);
  SPMatrix<SCFMode>& CIAO = CIAO_othoA.first;
  SPMatrix<SCFMode>& othoA = CIAO_othoA.second;
  auto minao = _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION);
  auto shells = minao->getBasis();
  auto orbToRef = getReferenceOrbitals(unpairedOrbitals, unpairedOrbitalsRef);

  Eigen::MatrixXd kinIntegrals;
  SpinPolarizedData<SCFMode, Eigen::VectorXd> kin_ref;
  if (_kineticAlign) {
    kinIntegrals = _system->getOneElectronIntegralController()->getKinIntegrals();
    kin_ref = getReferenceKineticEnergies();
  }

  for_spin(CIAO, pops_ref, orbToRef, C, othoA, kin_ref, orbitalRange) {
    unsigned int nOcc = orbToRef_spin.rows();

    if (_kineticAlign)
      kinIntegrals = (othoA_spin.transpose() * kinIntegrals * othoA_spin).eval();
    /*
     * 2x2 rotation loop
     */
    // Loop over pairs of unpaired orbitals
    unsigned int alignCycle = 0;
    while (true) {
      double largestAngle = 0.0;
      ++alignCycle;
      for (const auto i : orbitalRange_spin) {
        int iRef = orbToRef_spin[i];
        if (iRef < 0)
          continue;
        for (const auto j : orbitalRange_spin) {
          int jRef = orbToRef_spin[j];
          if (i == j)
            continue;
          unsigned int cycle = 0;
          double min;
          auto const updateFunction = [&](const Eigen::VectorXd& parameters, double& value, double& gradient,
                                          double& hessian, bool printInfo) {
            (void)printInfo;
            cycle++;
            gradient = 0.0;
            hessian = 0.0;
            value = 0.0;

            double x = parameters[0];
            if (std::fabs(x) > 100 * M_PI) {
              OutputControl::dOut << "Numerical minimization out of bounds! Break here!" << std::endl;
              return true;
            }
            // The penalty function is pi periodic. Adjust x accordingly.
            while (x > M_PI) {
              x -= M_PI;
            }
            while (x < -M_PI) {
              x += M_PI;
            }
            // Population-based criterion.
            // Loop over shells
            for (unsigned int shellIndex = 0; shellIndex < minao->getReducedNBasisFunctions(); ++shellIndex) {
              double Qii = 0.0;
              double Qjj = 0.0;
              double Qij = 0.0;
              double qi = pops_ref_spin(shellIndex, iRef);
              double qj = (jRef < 0) ? 0.0 : pops_ref_spin(shellIndex, jRef);
              // Loop over functions within the shell
              unsigned int nCont = shells[shellIndex]->getNContracted();
              unsigned int extStart = minao->extendedIndex(shellIndex);
              for (unsigned int mu = extStart; mu < extStart + nCont; ++mu) {
                Qii += CIAO_spin(mu, i) * CIAO_spin(mu, i);
                Qjj += CIAO_spin(mu, j) * CIAO_spin(mu, j);
                Qij += CIAO_spin(mu, i) * CIAO_spin(mu, j);
              } // for mu
              addToLagrangian(Qii, Qjj, Qij, qi, qj, x, jRef, value, gradient, hessian);
            } // for shellIndex

            // Orbital kinetic energy-based criterion.
            if (_kineticAlign) {
              const Eigen::VectorXd halfTransformedI = (kinIntegrals * CIAO_spin.col(i)).eval();
              const Eigen::VectorXd halfTransformedJ = (kinIntegrals * CIAO_spin.col(j)).eval();
              const double kii = CIAO_spin.col(i).transpose().eval() * halfTransformedI;
              const double kjj = CIAO_spin.col(j).transpose().eval() * halfTransformedJ;
              const double kij = CIAO_spin.col(i).transpose().eval() * halfTransformedJ;
              const double ki = kin_ref_spin[iRef];
              const double kj = (jRef < 0) ? 0.0 : kin_ref_spin[jRef];
              addToLagrangian(kii, kjj, kij, ki, kj, x, jRef, value, gradient, hessian);
            }
            min = hessian;
            bool converged = false;

            if (std::fabs(gradient) < 1e-8 && min > 0) {
              converged = true;
            }
            if (cycle == 100) {
              OutputControl::dOut << "Penalty result: " << value << " Numerical minimization did NOT CONVERGE!" << std::endl;
              OutputControl::dOut << "Continue with non-converged angle." << std::endl;
              converged = true;
            }
            return converged;
          };

          double start = 0.0;
          Eigen::VectorXd X = Eigen::VectorXd::Constant(1, start);
          auto optimizer = std::make_shared<NewtonRaphson>(X, 1e-6);
          optimizer->minimize2D(updateFunction);
          // Check whether the optimized angle is anywhere sensible.
          // This is not the most beautiful code but it may become necessary.
          // Check for extremely large angles. The penalty function is pi periodic,
          // so this should not happen.
          if (std::fabs(X[0]) > 100 * M_PI)
            X[0] = 0;
          while (X[0] > M_PI) {
            X[0] -= M_PI;
          }
          while (X[0] < -M_PI) {
            X[0] += M_PI;
          }
          // Check for nan. The BFGS managed to give me nan at some point,
          // so I better check for this.
          double angle = (X[0] == X[0]) ? X[0] : 0.0;
          if (std::fabs(angle) > largestAngle)
            largestAngle = angle;
          // Perform rotation.
          JacobiRotation::rotate(CIAO_spin.col(i), CIAO_spin.col(j), angle);
        } // for itJ<itI
      }   // for itI
      if (std::fabs(largestAngle) <= angleThreshold) {
        OutputControl::mOut << "    Converged after " << alignCycle << " orbital rotation cycles." << std::endl;
        C_spin.leftCols(nOcc) = othoA_spin * CIAO_spin;
        break;
      }
      else if (alignCycle == maxSweeps) {
        OutputControl::vOut << "    WARNING: IBO alignment procedure did not converged after " << alignCycle
                            << " orbital rotation cycles." << std::endl;
        OutputControl::vOut << "             The orbitals will still be stored for error analysis." << std::endl;
        OutputControl::vOut << "             Largest rotation angle in the last iteration: " << largestAngle << std::endl;
        OutputControl::vOut << "             This message can usually be ignored." << std::endl;
        C_spin.leftCols(nOcc) = othoA_spin * CIAO_spin;
        break;
      }
    } // while !converged
    const Eigen::MatrixXd S1 = _system->getOneElectronIntegralController()->getOverlapIntegrals();
    const double sanityCheck =
        (C_spin.leftCols(nOcc).transpose() * S1 * C_spin.leftCols(nOcc) - Eigen::MatrixXd::Identity(nOcc, nOcc))
            .array()
            .abs()
            .sum();
    if (sanityCheck > 1e-6)
      throw SerenityError("The orbitals are no longer orthogonal after orbital alignment!");
  }; /* spin loop */

  orbitals.updateOrbitals(C, orbitals.getEigenvalues());
  densMatCont->updateDensityMatrix();
  DensityMatrix<SCFMode> newP = densMatCont->getDensityMatrix();
  for_spin(newP, oldP) {
    auto diff = (oldP_spin - newP_spin).array().abs().maxCoeff();
    if (diff > 1e-9)
      OutputControl::mOut << "WARNING: Density matrix changed during orbital alignment!" << std::endl;
  };
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXi>
OrbitalAligner<SCFMode>::getReferenceOrbitals(const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitals,
                                              const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitalsRef) {
  SpinPolarizedData<SCFMode, Eigen::VectorXi> orbToRef;
  for_spin(orbToRef, unpairedOrbitals, unpairedOrbitalsRef) {
    orbToRef_spin = Eigen::VectorXi::Constant(unpairedOrbitals_spin.size(), -1);
    unsigned int iUnpaired = 0;
    for (unsigned int i = 0; i < unpairedOrbitals_spin.size(); ++i) {
      if (unpairedOrbitals_spin[i]) {
        unsigned int iRefUnpaired = 0;
        bool paired = false;
        for (unsigned int iRef = 0; iRef < unpairedOrbitalsRef_spin.size(); ++iRef) {
          if (unpairedOrbitalsRef_spin[iRef]) {
            if (iRefUnpaired == iUnpaired) {
              orbToRef_spin[i] = iRef;
              paired = true;
              break;
            } // if iRefUnpaired == iUnpaired
            ++iRefUnpaired;
          } // if !unpairedOrbitalsRef[iRef]
        }   // for iRef
        if (!paired)
          throw SerenityError("Wrong orbital pairing logic for alignment procedure.");
        ++iUnpaired;
      } // if !unpairedOrbitals_spin[i]
    }   // for i
  };
  return orbToRef;
}

template class OrbitalAligner<Options::SCF_MODES::RESTRICTED>;
template class OrbitalAligner<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
