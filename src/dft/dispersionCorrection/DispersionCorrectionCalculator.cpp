/**
 * @file   DispersionCorrectionCalculator.cpp
 *
 * @date   Nov 26, 2015
 * @author Jan Unsleber
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
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
#include "parameters/Constants.h"
/* Include Std and External Headers */
#include <cmath>
#include <iostream>

namespace Serenity {

/*
 * For references see DispersionCorrectionCalculator.h
 */

/*
 * ================================================
 *                     Energies
 * ================================================
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<>
double DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::NONE>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  return 0.0;
}
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<>
double DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::NONE>(
    std::shared_ptr<const Geometry> activeGeometry, std::shared_ptr<const Geometry> environmentGeometry,
    const CompositeFunctionals::XCFUNCTIONALS functional) {
  return 0.0;
}
#pragma GCC diagnostic pop

template<>
double DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();
  auto coordNums = calcCoordNumbers(geometry);
  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3>(functional, s6, rs6, s18, rs18, alp);
  double e6 = 0.0;
  double e8 = 0.0;
  // JU: Switching this loop structure to j<i makes a difference, it should not.
  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // Dispersion correction results only if the atom in question is not a ghost atom
        // TODO: MB: How to deal with ECP here? Do we still have to use the nucelar charge?
        if (nuclearChargeJ != 0)
          calculateD3Term(atomI, atomJ, coordNums[i], coordNums[j], rs6, rs18, alp, e6, e8);
      }
    }
  }
  return -e6 * s6 - e8 * s18;
}

void DispersionCorrectionCalculator::calculateD3Term(std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ,
                                                     const double& coordI, const double& coordJ, const double& rs6,
                                                     const double& rs18, const double& alp, double& e6, double& e8) {
  int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
  int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
  double r = distance(*atomI, *atomJ);
  if (r > DispersionData::RANGE_THRESHOLD_N2)
    return;
  double rr = DispersionData::r0ab(nuclearChargeI, nuclearChargeJ) / r;
  double damp6 = 1.0 / (1.0 + 6.0 * pow(rs6 * rr, alp));
  double c6 = getC6(atomI, atomJ, coordI, coordJ);
  double r6 = pow(r, 6);

  // r2r4 stored as sqrt
  double damp8 = 1.0 / (1.0 + 6.0 * pow(rs18 * rr, alp + 2));
  double c8 = 3.0 * c6 * DispersionData::r2r4(atomI->getAtomType()->getNuclearCharge()) *
              DispersionData::r2r4(atomJ->getAtomType()->getNuclearCharge());

  double r8 = r6 * r * r;

  e6 += c6 * damp6 / r6;
  e8 += c8 * damp8 / r8;
}

template<>
double DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3>(
    std::shared_ptr<const Geometry> activeGeometry, std::shared_ptr<const Geometry> environmentGeometry,
    const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto activeAtoms = activeGeometry->getAtoms();
  auto environmentAtoms = environmentGeometry->getAtoms();
  auto coordNumsAct = calcCoordNumbers(activeGeometry, environmentGeometry);
  auto coordNumsEnv = calcCoordNumbers(environmentGeometry, activeGeometry);
  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3>(functional, s6, rs6, s18, rs18, alp);
  double e6 = 0.0;
  double e8 = 0.0;
  // JU: Switching this loop structure to j<i makes a difference, it should not.
  for (unsigned int i = 0; i < activeGeometry->getNAtoms(); ++i) {
    auto atomI = activeAtoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < environmentGeometry->getNAtoms(); ++j) {
        auto atomJ = environmentAtoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // Dispersion correction results only if the atom in question is not a ghost atom
        // TODO: MB: How to deal with ECP here? Do we still have to use the nucelar charge?
        if (nuclearChargeJ != 0)
          calculateD3Term(atomI, atomJ, coordNumsAct[i], coordNumsEnv[j], rs6, rs18, alp, e6, e8);
      }
    }
  }
  return -e6 * s6 - e8 * s18;
}

template<>
double DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3ABC>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();

  auto coordNums = calcCoordNumbers(geometry);

  unsigned int nAtoms = geometry->getNAtoms();
  std::vector<bool> important(nAtoms * (nAtoms + 1) / 2, false);
  std::vector<double> distances(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> c6Sqrt(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> abcDamp(nAtoms * (nAtoms + 1) / 2, 0.0);

  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3>(functional, s6, rs6, s18, rs18, alp);

  double e6 = 0.0;
  double e8 = 0.0;
  double e6abc = 0.0;

  // JU: Switching this loop structure to j<i makes a difference, it should not.
  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // Dispersion correction results only if the atom in question is not a ghost atom
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          if (r > DispersionData::RANGE_THRESHOLD_N2)
            continue;
          double rr =
              DispersionData::r0ab(atomI->getAtomType()->getNuclearCharge(), atomJ->getAtomType()->getNuclearCharge()) / r;
          double damp6 = 1.0 / (1.0 + 6.0 * pow(rs6 * rr, alp));
          double c6 = getC6(atomI, atomJ, coordNums[i], coordNums[j]);
          double r6 = pow(r, 6);

          // r2r4 stored as sqrt
          double damp8 = 1.0 / (1.0 + 6.0 * pow(rs18 * rr, alp + 2));
          double c8 = 3.0 * c6 * DispersionData::r2r4(atomI->getAtomType()->getNuclearCharge()) *
                      DispersionData::r2r4(atomJ->getAtomType()->getNuclearCharge());

          double r8 = r6 * r * r;

          e6 += c6 * damp6 / r6;
          e8 += c8 * damp8 / r8;
          // store data for the 3 center correction
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            unsigned int ij = j + i * (i + 1) / 2;
            important[ij] = true;
            c6Sqrt[ij] = sqrt(c6);
            distances[ij] = r;
            double rr =
                DispersionData::r0ab(atomI->getAtomType()->getNuclearCharge(), atomJ->getAtomType()->getNuclearCharge()) / r;
            abcDamp[ij] = pow((1.0 / rr), (1.0 / 3.0));
          }
        }
      }
    }
  }
  // compute non-additive third-order energy using averaged C6
  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // Dispersion correction results only if the atom in question is not a ghost atom
        if (nuclearChargeJ != 0) {
          unsigned int ij = j + i * (i + 1) / 2;
          if (important[ij]) {
            for (unsigned int k = 0; k < j; ++k) {
              unsigned int ik = k + i * (i + 1) / 2;
              unsigned int jk = k + j * (j + 1) / 2;
              if (important[ik] && important[jk]) {
                // damping func product
                double tmp = 0.75 * abcDamp[ik] * abcDamp[jk] * abcDamp[ij];
                double damp9 = 1.0 / (1.0 + 6.0 * pow(tmp, (-alp - 2)));
                //  triple C6 coefficient (stored as sqrt)
                double c9 = c6Sqrt[ij] * c6Sqrt[ik] * c6Sqrt[jk];
                // angular terms  d is r^2

                double rij = distances[ij] * distances[ij];
                double rjk = distances[jk] * distances[jk];
                double rik = distances[ik] * distances[ik];
                double ang = 1.0 + 0.375 * (rij + rjk - rik) * (rij + rik - rjk) * (rik + rjk - rij) / (rij * rjk * rik);
                e6abc -= damp9 * c9 * ang / pow(rij * rjk * rik, 1.5);
              }
            }
          }
        }
      }
    }
  }

  return -e6 * s6 - e8 * s18 - e6abc;
}

template<>
double DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();

  auto coordNums = calcCoordNumbers(geometry);

  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(functional, s6, rs6, s18, rs18, alp);

  double e6 = 0.0;
  double e8 = 0.0;

  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // Dispersion correction results only if the atom in question is not a ghost atom
        if (nuclearChargeJ != 0)
          calculateD3BJTerm(atomI, atomJ, coordNums[i], coordNums[j], rs6, rs18, e6, e8);
      }
    }
  }
  return -e6 * s6 - e8 * s18;
}

template<>
double DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
    std::shared_ptr<const Geometry> activeGeometry, std::shared_ptr<const Geometry> environmentGeometry,
    const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto activeAtoms = activeGeometry->getAtoms();
  auto environmentAtoms = environmentGeometry->getAtoms();
  auto coordNumsAct = calcCoordNumbers(activeGeometry, environmentGeometry);
  auto coordNumsEnv = calcCoordNumbers(environmentGeometry, activeGeometry);
  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(functional, s6, rs6, s18, rs18, alp);

  double e6 = 0.0;
  double e8 = 0.0;

  for (unsigned int i = 0; i < activeGeometry->getNAtoms(); ++i) {
    auto atomI = activeAtoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < environmentGeometry->getNAtoms(); ++j) {
        auto atomJ = environmentAtoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // Dispersion correction results only if the atom in question is not a ghost atom
        if (nuclearChargeJ != 0) {
          calculateD3BJTerm(atomI, atomJ, coordNumsAct[i], coordNumsEnv[j], rs6, rs18, e6, e8);
        }
      }
    }
  }
  return -e6 * s6 - e8 * s18;
}

void DispersionCorrectionCalculator::DispersionCorrectionCalculator::calculateD3BJTerm(
    std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ, const double& coordI, const double& coordJ,
    const double& rs6, const double& rs18, double& e6, double& e8) {
  int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
  int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
  double r = distance(*atomI, *atomJ);
  if (r > DispersionData::RANGE_THRESHOLD_N2)
    return;

  double c6 = getC6(atomI, atomJ, coordI, coordJ);
  double r6 = pow(r, 6);

  // r2r4 stored as sqrt
  double c8 = 3.0 * c6 * DispersionData::r2r4(nuclearChargeI) * DispersionData::r2r4(nuclearChargeJ);
  double r8 = r6 * r * r;

  double tmp = sqrt(c8 / c6);
  e6 += c6 / (r6 + pow((rs6 * tmp + rs18), 6));
  e8 += c8 / (r8 + pow((rs6 * tmp + rs18), 8));
}

template<>
double DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJABC>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();

  auto coordNums = calcCoordNumbers(geometry);

  unsigned int nAtoms = geometry->getNAtoms();
  std::vector<bool> important(nAtoms * (nAtoms + 1) / 2, false);
  std::vector<double> distances(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> c6Sqrt(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> abcDamp(nAtoms * (nAtoms + 1) / 2, 0.0);

  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(functional, s6, rs6, s18, rs18, alp);

  double e6 = 0.0;
  double e8 = 0.0;
  double e6abc = 0.0;

  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          if (r > DispersionData::RANGE_THRESHOLD_N2)
            continue;

          double c6 = getC6(atomI, atomJ, coordNums[i], coordNums[j]);
          double r6 = pow(r, 6);

          // r2r4 stored as sqrt
          double c8 = 3.0 * c6 * DispersionData::r2r4(nuclearChargeI) * DispersionData::r2r4(nuclearChargeJ);
          double r8 = r6 * r * r;

          double tmp = sqrt(c8 / c6);
          e6 += c6 / (r6 + pow((rs6 * tmp + rs18), 6));
          e8 += c8 / (r8 + pow((rs6 * tmp + rs18), 8));

          // store data for the 3 center correction
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            unsigned int ij = j + i * (i + 1) / 2;
            important[ij] = true;
            c6Sqrt[ij] = sqrt(c6);
            distances[ij] = r;
            double rr = DispersionData::r0ab(nuclearChargeI, nuclearChargeJ) / r;
            abcDamp[ij] = pow((1.0 / rr), (1.0 / 3.0));
          }
        }
      }
    }
  }

  // compute non-additive third-order energy using averaged C6
  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          unsigned int ij = j + i * (i + 1) / 2;
          if (important[ij]) {
            for (unsigned int k = 0; k < j; ++k) {
              unsigned int ik = k + i * (i + 1) / 2;
              unsigned int jk = k + j * (j + 1) / 2;
              if (important[ik] && important[jk]) {
                // damping func product
                double tmp = 0.75 * abcDamp[ik] * abcDamp[jk] * abcDamp[ij];
                double damp9 = 1.0 / (1.0 + 6.0 * pow(tmp, (-alp - 2)));
                //  triple C6 coefficient (stored as sqrt)
                double c9 = c6Sqrt[ij] * c6Sqrt[ik] * c6Sqrt[jk];
                // angular terms  d is r^2

                double rij = distances[ij] * distances[ij];
                double rjk = distances[jk] * distances[jk];
                double rik = distances[ik] * distances[ik];
                double ang = 1.0 + 0.375 * (rij + rjk - rik) * (rij + rik - rjk) * (rik + rjk - rij) / (rij * rjk * rik);
                e6abc -= damp9 * c9 * ang / pow(rij * rjk * rik, 1.5);
              }
            }
          }
        }
      }
    }
  }
  return -e6 * s6 - e8 * s18 - e6abc;
}

double DispersionCorrectionCalculator::calcDispersionEnergyCorrection(Options::DFT_DISPERSION_CORRECTIONS dispType,
                                                                      std::shared_ptr<const Geometry> geometry,
                                                                      const CompositeFunctionals::XCFUNCTIONALS functional) {
  double disp = 0.0;
  switch (dispType) {
    case Options::DFT_DISPERSION_CORRECTIONS::NONE:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::NONE>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3ABC:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3ABC>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJ:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJABC:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJABC>(
          geometry, functional);
      break;
  }
  return disp;
}

double DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection(
    Options::DFT_DISPERSION_CORRECTIONS dispType, std::shared_ptr<const Geometry> activeGeometry,
    std::shared_ptr<const Geometry> environmentGeometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  double disp = 0.0;
  switch (dispType) {
    case Options::DFT_DISPERSION_CORRECTIONS::NONE:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::NONE>(
          activeGeometry, environmentGeometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3>(
          activeGeometry, environmentGeometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3ABC:
      throw SerenityError("D3ABC interaction type term is not implemented yet!");
      //        disp =
      //        DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3ABC>(
      //            activeGeometry,
      //            environmentGeometry,
      //            functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJ:
      disp = DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
          activeGeometry, environmentGeometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJABC:
      throw SerenityError("D3BJABC interaction type term is not implemented yet!");
      //        disp =
      //        DispersionCorrectionCalculator::calcDispersionEnergyInteractionCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJABC>(
      //            activeGeometry,
      //            environmentGeometry,
      //            functional);
      break;
  }
  return disp;
}
/*
 * ================================================
 *                     Gradients
 * ================================================
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<>
Eigen::MatrixXd DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::NONE>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  Eigen::MatrixXd grad(geometry->getNAtoms(), 3);
  grad.setZero();
  return grad;
}
#pragma GCC diagnostic pop

template<>
Eigen::MatrixXd DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();

  auto coordNums = calcCoordNumbers(geometry);

  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3>(functional, s6, rs6, s18, rs18, alp);

  unsigned int nAtoms = geometry->getNAtoms();
  Eigen::MatrixXd grad(nAtoms, 3);
  grad.setZero();
  std::vector<double> pairWiseGradCN(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> gradCN(nAtoms, 0.0);

  /*
   * Gradient of the energy w.r.t. the coordination number
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          if (r > DispersionData::RANGE_THRESHOLD_N2)
            continue;

          // Calc and read from raw data
          double c6 = getC6(atomI, atomJ, coordNums[i], coordNums[j]);
          auto deltaC6 = getDeltaC6(atomI, atomJ, coordNums[i], coordNums[j]);

          double r0ab = DispersionData::r0ab(nuclearChargeI, nuclearChargeJ);

          double r2r4Prod = DispersionData::r2r4(nuclearChargeI) * DispersionData::r2r4(nuclearChargeJ);

          // pre calculate factors and index
          unsigned int ij = j + i * (i + 1) / 2;
          double r6 = pow(r, 6);
          double t6 = pow((r / (rs6 * r0ab)), (-alp));
          double damp6 = 1.0 / (1.0 + 6.0 * t6);
          double t8 = pow((r / (rs18 * r0ab)), (-(alp + 2)));
          double damp8 = 1.0 / (1.0 + 6.0 * t8);

          /*
           * the following comments are taken from the original implementation
           */
          // derivative of 1/R^6 term w.r.t. drij : (dR^(-6)/drij)* fdamp * c6
          double tmp1 = s6 * 6.0 * damp6 * c6 / (r * r6);
          // derivative of 1/R^8 term w.r.t. drij : (dR^(-6)/drij)* fdamp * c8
          double tmp2 = s18 * 6.0 * c6 * r2r4Prod * damp8 / (r6 * r * r * r);

          // add to pairwise gradient
          pairWiseGradCN[ij] -= tmp1;
          pairWiseGradCN[ij] -= 4.0 * tmp2;

          // now the gradient of the damping function
          pairWiseGradCN[ij] += tmp1 * alp * t6 * damp6;             // C6 part
          pairWiseGradCN[ij] += 3.0 * tmp2 * (alp + 2) * t8 * damp8; // C8 part

          // now the last term of the 2-body gradient, i.e., the derivative dEdisp/dCNi = (dEdisp/dC6ij) * (dCij/dCNi)
          // //+ the dCNj term
          double dc6_rest = s6 * damp6 / r6 + 3.0 * s18 * r2r4Prod * damp8 / (r6 * r * r); // dEdisp/dC6ij

          //  for each atom: sum up the atomic contribution to (dEdisp/dC6ij)*(dC6ij/dCNi) and save for later
          gradCN[i] += dc6_rest * deltaC6.first;
          gradCN[j] += dc6_rest * deltaC6.second;
        }
      }
    }
  }

  /*
   * Coordination number gradient w.r.t. displacement
   *
   * and
   *
   * Gradient correction
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          /*
           * Original comments in place
           */
          // now compute derivative of coordination number w.r.t. displacement dCNi/drij
          double deltaCN = 0.0;
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            double rcovij = DispersionData::rcov(nuclearChargeI) + DispersionData::rcov(nuclearChargeJ);
            double tmp1 = exp(-DispersionData::K1 * (rcovij / r - 1.0));
            double tmp2 = (tmp1 + 1.0) * (tmp1 + 1.0) * r * r;
            deltaCN = -DispersionData::K1 * rcovij * tmp1 / tmp2;
          }

          // the total gradient w.r.t. the distance: dEdisp/drij
          double totalGrad = pairWiseGradCN[j + i * (i + 1) / 2] + deltaCN * (gradCN[i] + gradCN[j]);

          // compute the gradient in a particular Cartesian direction
          //   and add to gradient, e.g., (dEdisp/drij) * (drij/dxi)
          double gradX = totalGrad * (atomJ->getX() - atomI->getX()) / r;
          double gradY = totalGrad * (atomJ->getY() - atomI->getY()) / r;
          double gradZ = totalGrad * (atomJ->getZ() - atomI->getZ()) / r;

          grad(i, 0) += gradX;
          grad(i, 1) += gradY;
          grad(i, 2) += gradZ;
          grad(j, 0) -= gradX;
          grad(j, 1) -= gradY;
          grad(j, 2) -= gradZ;
        }
      }
    }
  }
  return grad;
}

template<>
Eigen::MatrixXd DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3ABC>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();

  auto coordNums = calcCoordNumbers(geometry);

  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3>(functional, s6, rs6, s18, rs18, alp);

  unsigned int nAtoms = geometry->getNAtoms();
  Eigen::MatrixXd grad(nAtoms, 3);
  grad.setZero();
  std::vector<double> pairWiseGradCN(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> gradCN(nAtoms, 0.0);

  std::vector<bool> important(nAtoms * (nAtoms + 1) / 2, false);
  std::vector<double> distances(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> c6saves(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> r3abc(nAtoms * (nAtoms + 1) / 2, 0.0);
  Eigen::MatrixXd dc6ij(nAtoms * (nAtoms + 1) / 2, nAtoms * (nAtoms + 1) / 2);
  dc6ij.setZero();

  /*
   * Gradient of the energy w.r.t. the coordination number
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          if (r > DispersionData::RANGE_THRESHOLD_N2)
            continue;

          // Calc and read from raw data
          double c6 = getC6(atomI, atomJ, coordNums[i], coordNums[j]);
          auto deltaC6 = getDeltaC6(atomI, atomJ, coordNums[i], coordNums[j]);

          double r0ab = DispersionData::r0ab(nuclearChargeI, nuclearChargeJ);

          double r2r4Prod = DispersionData::r2r4(nuclearChargeI) * DispersionData::r2r4(nuclearChargeJ);

          // pre calculate factors and index
          unsigned int ij = j + i * (i + 1) / 2;
          double r6 = pow(r, 6);
          double t6 = pow((r / (rs6 * r0ab)), (-alp));
          double damp6 = 1.0 / (1.0 + 6.0 * t6);
          double t8 = pow((r / (rs18 * r0ab)), (-(alp + 2)));
          double damp8 = 1.0 / (1.0 + 6.0 * t8);

          /*
           * the following comments are taken from the original implementation
           */
          // derivative of 1/R^6 term w.r.t. drij : (dR^(-6)/drij)* fdamp * c6
          double tmp1 = s6 * 6.0 * damp6 * c6 / (r * r6);
          // derivative of 1/R^8 term w.r.t. drij : (dR^(-6)/drij)* fdamp * c8
          double tmp2 = s18 * 6.0 * c6 * r2r4Prod * damp8 / (r6 * r * r * r);

          // add to pairwise gradient
          pairWiseGradCN[ij] -= tmp1;
          pairWiseGradCN[ij] -= 4.0 * tmp2;

          // now the gradient of the damping function
          pairWiseGradCN[ij] += tmp1 * alp * t6 * damp6;             // C6 part
          pairWiseGradCN[ij] += 3.0 * tmp2 * (alp + 2) * t8 * damp8; // C8 part

          // now the last term of the 2-body gradient, i.e., the derivative dEdisp/dCNi = (dEdisp/dC6ij) * (dCij/dCNi)
          // //+ the dCNj term
          double dc6_rest = s6 * damp6 / r6 + 3.0 * s18 * r2r4Prod * damp8 / (r6 * r * r); // dEdisp/dC6ij

          //  for each atom: sum up the atomic contribution to (dEdisp/dC6ij)*(dC6ij/dCNi) and save for later
          gradCN[i] += dc6_rest * deltaC6.first;
          gradCN[j] += dc6_rest * deltaC6.second;

          // save data for the abc correction
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            important[ij] = true;
            dc6ij(j, i) = deltaC6.second;
            dc6ij(i, j) = deltaC6.first;
            c6saves[ij] = c6;
            distances[ij] = r;
            r3abc[ij] = pow(r / r0ab, 1.0 / 3.0); // geometric mean
          }
        }
      }
    }
  }

  /*
   * Gradient of the three body correction term
   */
  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          unsigned int ij = j + i * (i + 1) / 2;
          if (important[ij]) {
            for (unsigned int k = 0; k < j; ++k) {
              unsigned int ik = k + i * (i + 1) / 2;
              unsigned int jk = k + j * (j + 1) / 2;
              if (important[ik] && important[jk]) {
                double r2ij = distances[ij] * distances[ij];
                double r2jk = distances[jk] * distances[jk];
                double r2ik = distances[ik] * distances[ik];
                // calculate C9 coefficient
                double c9 = sqrt(c6saves[jk] * c6saves[ik] * c6saves[ij]);
                // set-up damping function
                double tmp = pow(0.75 * r3abc[jk] * r3abc[ik] * r3abc[ij], (-alp - 2));
                double damp = 1.0 / (1.0 + 6.0 * tmp);
                double rijk3 = r2ij * r2jk * r2ik;
                double rijk3Pow = pow(rijk3, 1.5);
                // compute angular term
                double ang = 1.0 + 0.375 * (r2ij + r2jk - r2ik) * (r2ij + r2ik - r2jk) * (r2ik + r2jk - r2ij) /
                                       (r2ij * r2jk * r2ik);
                ang /= rijk3Pow;

                // now gradients
                double dfdmp = -2.0 * (alp + 2) * tmp * damp * damp; // derivative of damping function

                double dang = 0.0;
                // j-k
                // calculate the derivatives of each part w.r.t. rjk
                dang = 3.0 * pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * pow(r2ij, 2);
                dang *= r2jk;
                dang -= 5.0 * pow((r2ik - r2ij), 2) * (r2ik + r2ij);
                dang += pow(r2jk, 3) + pow(r2jk, 2) * (r2ik + r2ij);
                dang *= (-0.375);
                dang /= (distances[jk] * rijk3 * rijk3Pow);

                pairWiseGradCN[jk] += c9 * ang * dfdmp / distances[jk] - dang * c9 * damp;

                // i-k
                // calculate the derivatives of each part w.r.t. rik
                dang = 3.0 * pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * pow(r2ij, 2);
                dang *= r2ik;
                dang -= 5.0 * pow((r2jk - r2ij), 2) * (r2jk + r2ij);
                dang += pow(r2ik, 3) + pow(r2ik, 2) * (r2jk + r2ij);
                dang *= (-0.375);
                dang /= (distances[ik] * rijk3 * rijk3Pow);

                pairWiseGradCN[ik] += c9 * ang * dfdmp / distances[ik] - dang * c9 * damp;

                // i-j
                // calculate the derivatives of each part w.r.t. rij
                dang = 3.0 * pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * pow(r2ik, 2);
                dang *= r2ij;
                dang -= 5.0 * pow((r2jk - r2ik), 2) * (r2jk + r2ik);
                dang += pow(r2ij, 3) + pow(r2ij, 2) * (r2jk + r2ik);
                dang *= (-0.375);
                dang /= (distances[ij] * rijk3 * rijk3Pow);

                pairWiseGradCN[ij] += c9 * ang * dfdmp / distances[ij] - dang * c9 * damp;

                // calculate rest* dc9/dcn(iat)  and sum it up for every atom ijk (i.., dc6i[...])
                double dc6_rest = ang * damp;

                gradCN[k] += dc6_rest * ((dc6ij(k, j) / c6saves[jk]) + (dc6ij(k, i) / c6saves[ik])) * (-0.5 * c9);
                gradCN[j] += dc6_rest * ((dc6ij(j, k) / c6saves[jk]) + (dc6ij(j, i) / c6saves[ij])) * (-0.5 * c9);
                gradCN[i] += dc6_rest * ((dc6ij(i, k) / c6saves[ik]) + (dc6ij(i, j) / c6saves[ij])) * (-0.5 * c9);
              }
            }
          }
        }
      }
    }
  }

  /*
   * Coordination number gradient w.r.t. displacement
   *
   * and
   *
   * Gradient correction
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          /*
           * Original comments in place
           */
          // now compute derivative of coordination number w.r.t. displacement dCNi/drij
          double deltaCN = 0.0;
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            double rcovij = DispersionData::rcov(nuclearChargeI) + DispersionData::rcov(nuclearChargeJ);
            double tmp1 = exp(-DispersionData::K1 * (rcovij / r - 1.0));
            double tmp2 = (tmp1 + 1.0) * (tmp1 + 1.0) * r * r;
            deltaCN = -DispersionData::K1 * rcovij * tmp1 / tmp2;
          }

          // the total gradient w.r.t. the distance: dEdisp/drij
          double totalGrad = pairWiseGradCN[j + i * (i + 1) / 2] + deltaCN * (gradCN[i] + gradCN[j]);

          // compute the gradient in a particular Cartesian direction
          //   and add to gradient, e.g., (dEdisp/drij) * (drij/dxi)
          double gradX = totalGrad * (atomJ->getX() - atomI->getX()) / r;
          double gradY = totalGrad * (atomJ->getY() - atomI->getY()) / r;
          double gradZ = totalGrad * (atomJ->getZ() - atomI->getZ()) / r;

          grad(i, 0) += gradX;
          grad(i, 1) += gradY;
          grad(i, 2) += gradZ;
          grad(j, 0) -= gradX;
          grad(j, 1) -= gradY;
          grad(j, 2) -= gradZ;
        }
      }
    }
  }
  return grad;
}

template<>
Eigen::MatrixXd DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();

  auto coordNums = calcCoordNumbers(geometry);

  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(functional, s6, rs6, s18, rs18, alp);

  unsigned int nAtoms = geometry->getNAtoms();
  Eigen::MatrixXd grad(nAtoms, 3);
  grad.setZero();
  std::vector<double> pairWiseGradCN(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> gradCN(nAtoms, 0.0);

  /*
   * Gradient of the energy w.r.t. the coordination number
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          if (r > DispersionData::RANGE_THRESHOLD_N2)
            continue;

          // Calc and read from raw data
          double c6 = getC6(atomI, atomJ, coordNums[i], coordNums[j]);
          auto deltaC6 = getDeltaC6(atomI, atomJ, coordNums[i], coordNums[j]);

          double r2r4Prod = DispersionData::r2r4(nuclearChargeI) * DispersionData::r2r4(nuclearChargeJ);

          // pre calculate factors and index
          double r0 = rs6 * sqrt(3.0 * r2r4Prod) + rs18;
          unsigned int ij = j + i * (i + 1) / 2;
          double r5 = pow(r, 5);
          // t6, t8 and tmp1 are intermediates for BJ-damping
          double t6 = r5 * r + pow(r0, 6);
          double t8 = r5 * r * r * r + pow(r0, 8);
          double tmp1 = 6.0 * c6 * r5;

          /*
           * the following comments are taken from the original implementation
           */
          // derivative of 1/R^6 + damping term w.r.t. drij
          pairWiseGradCN[ij] -= s6 * tmp1 / (t6 * t6);
          // derivative of 1/R^8 + damping term w.r.t. drij
          pairWiseGradCN[ij] -= 4.0 * s18 * r2r4Prod * r * r * tmp1 / (t8 * t8);

          // now the last term of the 2-body gradient, i.e., the derivative dEdisp/dCNi = (dEdisp/dC6ij) * (dCij/dCNi)
          // //+ the dCNj term
          double dc6_rest = s6 / t6 + 3.0 * s18 * r2r4Prod / t8; // dEdisp/dC6ij

          //  for each atom: sum up the atomic contribution to (dEdisp/dC6ij)*(dC6ij/dCNi) and save for later
          gradCN[i] += dc6_rest * deltaC6.first;
          gradCN[j] += dc6_rest * deltaC6.second;
        }
      }
    }
  }

  /*
   * Coordination number gradient w.r.t. displacement
   *
   * and
   *
   * Gradient correction
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          /*
           * Original comments in place
           */
          // now compute derivative of coordination number w.r.t. displacement dCNi/drij
          double deltaCN = 0.0;
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            double rcovij = DispersionData::rcov(nuclearChargeI) + DispersionData::rcov(nuclearChargeJ);
            double tmp1 = exp(-DispersionData::K1 * (rcovij / r - 1.0));
            double tmp2 = (tmp1 + 1.0) * (tmp1 + 1.0) * r * r;
            deltaCN = -DispersionData::K1 * rcovij * tmp1 / tmp2;
          }

          // the total gradient w.r.t. the distance: dEdisp/drij
          double totalGrad = pairWiseGradCN[j + i * (i + 1) / 2] + deltaCN * (gradCN[i] + gradCN[j]);
          // compute the gradient in a particular Cartesian direction
          //   and add to gradient, e.g., (dEdisp/drij) * (drij/dxi)
          double gradX = totalGrad * (atomJ->getX() - atomI->getX()) / r;
          double gradY = totalGrad * (atomJ->getY() - atomI->getY()) / r;
          double gradZ = totalGrad * (atomJ->getZ() - atomI->getZ()) / r;

          grad(i, 0) += gradX;
          grad(i, 1) += gradY;
          grad(i, 2) += gradZ;
          grad(j, 0) -= gradX;
          grad(j, 1) -= gradY;
          grad(j, 2) -= gradZ;
        }
      }
    }
  }
  return grad;
}

template<>
Eigen::MatrixXd DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJABC>(
    std::shared_ptr<const Geometry> geometry, const CompositeFunctionals::XCFUNCTIONALS functional) {
  auto atoms = geometry->getAtoms();

  auto coordNums = calcCoordNumbers(geometry);

  double s6, rs6, s18, rs18, alp;
  DispersionData::getFunctionalParameters<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(functional, s6, rs6, s18, rs18, alp);

  unsigned int nAtoms = geometry->getNAtoms();
  Eigen::MatrixXd grad(nAtoms, 3);
  grad.setZero();
  std::vector<double> pairWiseGradCN(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> gradCN(nAtoms, 0.0);

  std::vector<bool> important(nAtoms * (nAtoms + 1) / 2, false);
  std::vector<double> distances(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> c6saves(nAtoms * (nAtoms + 1) / 2, 0.0);
  std::vector<double> r3abc(nAtoms * (nAtoms + 1) / 2, 0.0);
  Eigen::MatrixXd dc6ij(nAtoms * (nAtoms + 1) / 2, nAtoms * (nAtoms + 1) / 2);
  dc6ij.setZero();

  /*
   * Gradient of the energy w.r.t. the coordination number
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          if (r > DispersionData::RANGE_THRESHOLD_N2)
            continue;

          // Calc and read from raw data
          double c6 = getC6(atomI, atomJ, coordNums[i], coordNums[j]);
          auto deltaC6 = getDeltaC6(atomI, atomJ, coordNums[i], coordNums[j]);

          double r2r4Prod = DispersionData::r2r4(nuclearChargeI) * DispersionData::r2r4(nuclearChargeJ);

          // pre calculate factors and index
          double r0 = rs6 * sqrt(3.0 * r2r4Prod) + rs18;
          unsigned int ij = j + i * (i + 1) / 2;
          double r5 = pow(r, 5);
          // t6, t8 and tmp1 are intermediates for BJ-damping
          double t6 = r5 * r + pow(r0, 6);
          double t8 = r5 * r * r * r + pow(r0, 8);
          double tmp1 = 6.0 * c6 * r5;

          /*
           * the following comments are taken from the original implementation
           */
          // derivative of 1/R^6 + damping term w.r.t. drij
          pairWiseGradCN[ij] -= s6 * tmp1 / (t6 * t6);
          // derivative of 1/R^8 + damping term w.r.t. drij
          pairWiseGradCN[ij] -= 4.0 * s18 * r2r4Prod * r * r * tmp1 / (t8 * t8);

          // now the last term of the 2-body gradient, i.e., the derivative dEdisp/dCNi = (dEdisp/dC6ij) * (dCij/dCNi)
          // //+ the dCNj term
          double dc6_rest = s6 / t6 + 3.0 * s18 * r2r4Prod / t8; // dEdisp/dC6ij

          //  for each atom: sum up the atomic contribution to (dEdisp/dC6ij)*(dC6ij/dCNi) and save for later
          gradCN[i] += dc6_rest * deltaC6.first;
          gradCN[j] += dc6_rest * deltaC6.second;

          // save data for the abc correction
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            important[ij] = true;
            dc6ij(j, i) = deltaC6.second;
            dc6ij(i, j) = deltaC6.first;
            c6saves[ij] = c6;
            distances[ij] = r;
            r0 = DispersionData::r0ab(nuclearChargeI, nuclearChargeJ);
            r3abc[ij] = pow(r / r0, 1.0 / 3.0); // geometric mean
          }
        }
      }
    }
  }

  /*
   * Gradient of the three body correction term
   */
  for (unsigned int i = 0; i < geometry->getNAtoms(); ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          unsigned int ij = j + i * (i + 1) / 2;
          if (important[ij]) {
            for (unsigned int k = 0; k < j; ++k) {
              unsigned int ik = k + i * (i + 1) / 2;
              unsigned int jk = k + j * (j + 1) / 2;
              if (important[ik] && important[jk]) {
                double r2ij = distances[ij] * distances[ij];
                double r2jk = distances[jk] * distances[jk];
                double r2ik = distances[ik] * distances[ik];
                // calculate C9 coefficient
                double c9 = sqrt(c6saves[jk] * c6saves[ik] * c6saves[ij]);
                // set-up damping function
                double tmp = pow(0.75 * r3abc[jk] * r3abc[ik] * r3abc[ij], (-alp - 2));
                double damp = 1.0 / (1.0 + 6.0 * tmp);
                double rijk3 = r2ij * r2jk * r2ik;
                double rijk3Pow = pow(rijk3, 1.5);
                // compute angular term
                double ang = 1.0 + 0.375 * (r2ij + r2jk - r2ik) * (r2ij + r2ik - r2jk) * (r2ik + r2jk - r2ij) /
                                       (r2ij * r2jk * r2ik);
                ang /= rijk3Pow;

                // now gradients
                double dfdmp = -2.0 * (alp + 2) * tmp * damp * damp; // derivative of damping function

                double dang = 0.0;
                // j-k
                // calculate the derivatives of each part w.r.t. rjk
                dang = 3.0 * pow(r2ik, 2) + 2.0 * r2ik * r2ij + 3.0 * pow(r2ij, 2);
                dang *= r2jk;
                dang -= 5.0 * pow((r2ik - r2ij), 2) * (r2ik + r2ij);
                dang += pow(r2jk, 3) + pow(r2jk, 2) * (r2ik + r2ij);
                dang *= (-0.375);
                dang /= (distances[jk] * rijk3 * rijk3Pow);

                pairWiseGradCN[jk] += c9 * ang * dfdmp / distances[jk] - dang * c9 * damp;

                // i-k
                // calculate the derivatives of each part w.r.t. rik
                dang = 3.0 * pow(r2jk, 2) + 2.0 * r2jk * r2ij + 3.0 * pow(r2ij, 2);
                dang *= r2ik;
                dang -= 5.0 * pow((r2jk - r2ij), 2) * (r2jk + r2ij);
                dang += pow(r2ik, 3) + pow(r2ik, 2) * (r2jk + r2ij);
                dang *= (-0.375);
                dang /= (distances[ik] * rijk3 * rijk3Pow);

                pairWiseGradCN[ik] += c9 * ang * dfdmp / distances[ik] - dang * c9 * damp;

                // i-j
                // calculate the derivatives of each part w.r.t. rij
                dang = 3.0 * pow(r2jk, 2) + 2.0 * r2jk * r2ik + 3.0 * pow(r2ik, 2);
                dang *= r2ij;
                dang -= 5.0 * pow((r2jk - r2ik), 2) * (r2jk + r2ik);
                dang += pow(r2ij, 3) + pow(r2ij, 2) * (r2jk + r2ik);
                dang *= (-0.375);
                dang /= (distances[ij] * rijk3 * rijk3Pow);

                pairWiseGradCN[ij] += c9 * ang * dfdmp / distances[ij] - dang * c9 * damp;

                // calculate rest* dc9/dcn(iat)  and sum it up for every atom ijk (i.., dc6i[...])
                double dc6_rest = ang * damp;

                gradCN[k] += dc6_rest * ((dc6ij(k, j) / c6saves[jk]) + (dc6ij(k, i) / c6saves[ik])) * (-0.5 * c9);
                gradCN[j] += dc6_rest * ((dc6ij(j, k) / c6saves[jk]) + (dc6ij(j, i) / c6saves[ij])) * (-0.5 * c9);
                gradCN[i] += dc6_rest * ((dc6ij(i, k) / c6saves[ik]) + (dc6ij(i, j) / c6saves[ij])) * (-0.5 * c9);
              }
            }
          }
        }
      }
    }
  }

  /*
   * Coordination number gradient w.r.t. displacement
   *
   * and
   *
   * Gradient correction
   */
  for (unsigned int i = 0; i < nAtoms; ++i) {
    auto atomI = atoms[i];
    int nuclearChargeI = atomI->getAtomType()->getNuclearCharge();
    // no dispersion correction for ghost atoms.
    if (nuclearChargeI != 0) {
      for (unsigned int j = 0; j < i; ++j) {
        auto atomJ = atoms[j];
        int nuclearChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // no dispersion correction for ghost atoms.
        if (nuclearChargeJ != 0) {
          double r = distance(*atomI, *atomJ);
          /*
           * Original comments in place
           */
          // now compute derivative of coordination number w.r.t. displacement dCNi/drij
          double deltaCN = 0.0;
          if (r < DispersionData::RANGE_THRESHOLD_N3) {
            double rcovij = DispersionData::rcov(nuclearChargeI) + DispersionData::rcov(nuclearChargeJ);
            double tmp1 = exp(-DispersionData::K1 * (rcovij / r - 1.0));
            double tmp2 = (tmp1 + 1.0) * (tmp1 + 1.0) * r * r;
            deltaCN = -DispersionData::K1 * rcovij * tmp1 / tmp2;
          }

          // the total gradient w.r.t. the distance: dEdisp/drij
          double totalGrad = pairWiseGradCN[j + i * (i + 1) / 2] + deltaCN * (gradCN[i] + gradCN[j]);
          // compute the gradient in a particular Cartesian direction
          //   and add to gradient, e.g., (dEdisp/drij) * (drij/dxi)
          double gradX = totalGrad * (atomJ->getX() - atomI->getX()) / r;
          double gradY = totalGrad * (atomJ->getY() - atomI->getY()) / r;
          double gradZ = totalGrad * (atomJ->getZ() - atomI->getZ()) / r;

          grad(i, 0) += gradX;
          grad(i, 1) += gradY;
          grad(i, 2) += gradZ;
          grad(j, 0) -= gradX;
          grad(j, 1) -= gradY;
          grad(j, 2) -= gradZ;
        }
      }
    }
  }
  return grad;
}

Eigen::MatrixXd
DispersionCorrectionCalculator::calcDispersionGradientCorrection(Options::DFT_DISPERSION_CORRECTIONS dispType,
                                                                 std::shared_ptr<const Geometry> geometry,
                                                                 const CompositeFunctionals::XCFUNCTIONALS functional) {
  Eigen::MatrixXd grad;
  switch (dispType) {
    case Options::DFT_DISPERSION_CORRECTIONS::NONE:
      grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::NONE>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3:
      grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3ABC:
      grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3ABC>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJ:
      grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJABC:
      grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJABC>(
          geometry, functional);
      break;
  }
  return grad;
}

/*
 * ================================================
 *                     Hessian
 * ================================================
 */

template<Options::DFT_DISPERSION_CORRECTIONS SCFMode>
std::vector<Eigen::MatrixXd>
DispersionCorrectionCalculator::calcDispersionHessianCorrection(std::shared_ptr<const Geometry> geometry,
                                                                const CompositeFunctionals::XCFUNCTIONALS functional) {
  const double delta = 0.05;
  const unsigned int nAtoms = geometry->getNAtoms();
  std::vector<Eigen::MatrixXd> hessian;

  // make a deep copy of the geometry, because the atoms will be moved
  //   and we want the integrals and gird to stay alive
  std::vector<std::shared_ptr<Atom>> atomCopies(0, nullptr);
  for (auto atom : geometry->getAtoms()) {
    atomCopies.push_back(
        std::make_shared<Atom>(atom->getAtomType()->getElementSymbol(), atom->getX(), atom->getY(), atom->getZ()));
  }
  auto geometryCopy = std::make_shared<const Geometry>(atomCopies);

  /*
   * XX,XY,XZ
   */
  Eigen::MatrixXd xx(nAtoms, nAtoms);
  Eigen::MatrixXd xy(nAtoms, nAtoms);
  Eigen::MatrixXd xz(nAtoms, nAtoms);
  xx.setZero();
  xy.setZero();
  xz.setZero();
  for (unsigned int i = 0; i != nAtoms; ++i) {
    geometryCopy->getAtoms()[i]->addToX(delta);
    auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<SCFMode>(geometryCopy, functional);
    for (unsigned int j = 0; j != nAtoms; ++j) {
      xx(i, j) += grad(j, 0);
      xy(i, j) += grad(j, 1);
      xz(i, j) += grad(j, 2);
    }
    geometryCopy->getAtoms()[i]->addToX(-2.0 * delta);
    grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<SCFMode>(geometryCopy, functional);
    for (unsigned int j = 0; j != nAtoms; ++j) {
      xx(i, j) -= grad(j, 0);
      xy(i, j) -= grad(j, 1);
      xz(i, j) -= grad(j, 2);
      xx(i, j) /= 2.0 * delta;
      xy(i, j) /= 2.0 * delta;
      xz(i, j) /= 2.0 * delta;
    }
    geometryCopy->getAtoms()[i]->addToX(delta);
  }
  hessian.push_back(xx);
  hessian.push_back(xy);
  hessian.push_back(xz);
  /*
   * YY,YZ
   */
  Eigen::MatrixXd yy(nAtoms, nAtoms);
  Eigen::MatrixXd yz(nAtoms, nAtoms);
  yy.setZero();
  yz.setZero();
  for (unsigned int i = 0; i != nAtoms; ++i) {
    geometryCopy->getAtoms()[i]->addToY(delta);
    auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<SCFMode>(geometryCopy, functional);
    for (unsigned int j = 0; j != nAtoms; ++j) {
      yy(i, j) += grad(j, 1);
      yz(i, j) += grad(j, 2);
    }
    geometryCopy->getAtoms()[i]->addToY(-2.0 * delta);
    grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<SCFMode>(geometryCopy, functional);
    for (unsigned int j = 0; j != nAtoms; ++j) {
      yy(i, j) -= grad(j, 1);
      yz(i, j) -= grad(j, 2);
      yy(i, j) /= 2.0 * delta;
      yz(i, j) /= 2.0 * delta;
    }
    geometryCopy->getAtoms()[i]->addToY(delta);
  }
  hessian.push_back(yy);
  hessian.push_back(yz);
  /*
   * ZZ
   */
  Eigen::MatrixXd zz(nAtoms, nAtoms);
  zz.setZero();
  for (unsigned int i = 0; i != nAtoms; ++i) {
    geometryCopy->getAtoms()[i]->addToZ(delta);
    auto grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<SCFMode>(geometryCopy, functional);
    for (unsigned int j = 0; j != nAtoms; ++j) {
      zz(i, j) += grad(j, 0);
    }
    geometryCopy->getAtoms()[i]->addToZ(-2.0 * delta);
    grad = DispersionCorrectionCalculator::calcDispersionGradientCorrection<SCFMode>(geometryCopy, functional);
    for (unsigned int j = 0; j != nAtoms; ++j) {
      zz(i, j) -= grad(j, 0);
      zz(i, j) /= 2.0 * delta;
    }
    geometryCopy->getAtoms()[i]->addToZ(delta);
  }
  hessian.push_back(zz);
  return hessian;
}

std::vector<Eigen::MatrixXd>
DispersionCorrectionCalculator::calcDispersionHessianCorrection(Options::DFT_DISPERSION_CORRECTIONS dispType,
                                                                std::shared_ptr<const Geometry> geometry,
                                                                const CompositeFunctionals::XCFUNCTIONALS functional) {
  std::vector<Eigen::MatrixXd> hessian;
  switch (dispType) {
    case Options::DFT_DISPERSION_CORRECTIONS::NONE:
      hessian = DispersionCorrectionCalculator::calcDispersionHessianCorrection<Options::DFT_DISPERSION_CORRECTIONS::NONE>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3:
      hessian = DispersionCorrectionCalculator::calcDispersionHessianCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3ABC:
      hessian = DispersionCorrectionCalculator::calcDispersionHessianCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3ABC>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJ:
      hessian = DispersionCorrectionCalculator::calcDispersionHessianCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJ>(
          geometry, functional);
      break;
    case Options::DFT_DISPERSION_CORRECTIONS::D3BJABC:
      hessian = DispersionCorrectionCalculator::calcDispersionHessianCorrection<Options::DFT_DISPERSION_CORRECTIONS::D3BJABC>(
          geometry, functional);
      break;
  }
  return hessian;
}

/*
 * ================================================
 *                 Private Helper
 * ================================================
 */

double DispersionCorrectionCalculator::getC6(std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ,
                                             const double& nCoordI, const double& nCoordJ) {
  double c6abMem = 0.0;
  double rsum = 0.0;
  double csum = 0.0;
  int zI = atomI->getAtomType()->getNuclearCharge();
  int zJ = atomJ->getAtomType()->getNuclearCharge();
  for (unsigned int i = 0; i <= DispersionData::maxNC6(zI); ++i) {
    for (unsigned int j = 0; j <= DispersionData::maxNC6(zJ); ++j) {
      double c6ab = DispersionData::c6ab(zI, zJ, i, j, 0);
      if (c6ab > 0.0) {
        c6abMem = c6ab;
        double tmp1 = pow((DispersionData::c6ab(zI, zJ, i, j, 1) - nCoordI), 2);
        tmp1 += pow((DispersionData::c6ab(zI, zJ, i, j, 2) - nCoordJ), 2);
        tmp1 = exp(DispersionData::K3 * tmp1);

        rsum += tmp1;
        csum += tmp1 * c6ab;
      }
    }
  }
  return (rsum > 0.0) ? csum / rsum : c6abMem;
}

std::pair<double, double> DispersionCorrectionCalculator::getDeltaC6(std::shared_ptr<Atom> atomI, std::shared_ptr<Atom> atomJ,
                                                                     double& nCoordI, double& nCoordJ) {
  double rsum = 0.0;
  double csum = 0.0;
  double deltaRsumI = 0.0;
  double deltaRsumJ = 0.0;
  double deltaCsumI = 0.0;
  double deltaCsumJ = 0.0;

  int zI = atomI->getAtomType()->getNuclearCharge();
  int zJ = atomJ->getAtomType()->getNuclearCharge();
  for (unsigned int i = 0; i <= DispersionData::maxNC6(zI); ++i) {
    for (unsigned int j = 0; j <= DispersionData::maxNC6(zJ); ++j) {
      double c6ab = DispersionData::c6ab(zI, zJ, i, j, 0);
      if (c6ab > 0.0) {
        double tmpI = DispersionData::c6ab(zI, zJ, i, j, 1) - nCoordI;
        double tmpJ = DispersionData::c6ab(zI, zJ, i, j, 2) - nCoordJ;
        double sumPowIJ = tmpI * tmpI + tmpJ * tmpJ;
        double expIJ = exp(DispersionData::K3 * sumPowIJ);
        rsum += expIJ;
        csum += expIJ * c6ab;

        deltaCsumI -= c6ab * expIJ * 2.0 * DispersionData::K3 * tmpI;
        deltaCsumJ -= c6ab * expIJ * 2.0 * DispersionData::K3 * tmpJ;

        deltaRsumI -= expIJ * 2.0 * DispersionData::K3 * tmpI;
        deltaRsumJ -= expIJ * 2.0 * DispersionData::K3 * tmpJ;
      }
    }
  }
  return (rsum > 0.0) ? std::make_pair(((deltaCsumI * rsum) - (deltaRsumI * csum)) / (rsum * rsum),
                                       ((deltaCsumJ * rsum) - (deltaRsumJ * csum)) / (rsum * rsum))
                      : std::make_pair(0.0, 0.0);
}

std::vector<double> DispersionCorrectionCalculator::calcCoordNumbers(std::shared_ptr<const Geometry> geometry,
                                                                     std::shared_ptr<const Geometry> environmentGeometry) {
  auto atoms = geometry->getAtoms();
  std::vector<double> coordNumbers;
  auto totalAtoms = atoms;
  if (environmentGeometry)
    totalAtoms.insert(totalAtoms.end(), atoms.begin(), atoms.end());

  for (auto atomI : atoms) {
    double coordNum = 0.0;
    int atomChargeI = atomI->getAtomType()->getNuclearCharge();
    // only not ghost atoms can have a coordination number > 0
    if (atomChargeI != 0) {
      for (auto atomJ : totalAtoms) {
        int atomChargeJ = atomJ->getAtomType()->getNuclearCharge();
        // Only if the atom is not a ghost atom, it can add to the coordination number.
        if (atomI != atomJ and atomChargeJ != 0) {
          double tmp = DispersionData::rcov(atomI->getAtomType()->getNuclearCharge());
          tmp += DispersionData::rcov(atomJ->getAtomType()->getNuclearCharge());
          tmp /= distance(*atomI, *atomJ);
          // original comment: counting function exponential has a better long-range behavior than MHGs inverse damping
          coordNum += 1.0 / (1.0 + exp(-1.0 * (DispersionData::K1) * (tmp - 1.0)));
        } // if I!=J and atomChargeJ!=0
      }   // for atomJ
    }     // if atomChargeI!=0
    coordNumbers.push_back(coordNum);
  } // for atomI
  return coordNumbers;
}

} /* namespace Serenity */
