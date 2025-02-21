/**
 * @file MultipoleMomentTask.cpp
 *
 * @date Nov 23, 2015
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
#include "tasks/MultipoleMomentTask.h"
/* Include Serenity Internal Headers */
#include "analysis/multipoles/MultipoleMomentCalculator.h"
#include "analysis/multipoles/NumericalDipoleMomentCalculator.h"
#include "geometry/Geometry.h"
#include "io/FormattedOutput.h"
#include "parameters/Constants.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <vector>

namespace Serenity {

MultipoleMomentTask::MultipoleMomentTask(std::vector<std::shared_ptr<SystemController>> activeSystems)
  : _activeSystems(activeSystems) {
}

void MultipoleMomentTask::run() {
  _origin = Point(0, 0, 0);
  if (settings.origin == Options::GAUGE_ORIGIN::COM) {
    printSmallCaption("Choosing the (supersystem) center of mass as the origin");
    Geometry supersystemGeometry;
    for (auto sys : _activeSystems) {
      supersystemGeometry += *sys->getGeometry();
    }
    supersystemGeometry.deleteIdenticalAtoms();
    _origin = supersystemGeometry.getCenterOfMass();
  }

  if (settings.highestOrder > 2)
    throw SerenityError("Only dipoles and quadrupoles are available.");

  printf("\n");
  printf("%4s Origin chosen as:\n", "");
  printf("%4s (%+3.5f, %+3.5f, %+3.5f)\n\n", " ", _origin.getX(), _origin.getY(), _origin.getZ());
  printf("%4s (Will be used for the origin correction of multipoles for charged systems!)\n\n", "");

  Eigen::Vector3d totalDipole = Eigen::Vector3d::Zero(3);
  Eigen::MatrixXd totalQuadrupole = Eigen::MatrixXd::Zero(2, 3);
  for (auto system : _activeSystems) {
    Eigen::Vector3d dipole = Eigen::Vector3d::Zero(3);
    Eigen::MatrixXd quadrupole = Eigen::MatrixXd::Zero(2, 3);
    calculateMultipoleMoments(dipole, quadrupole, system);
    totalDipole += dipole;
    totalQuadrupole += quadrupole;
  } // for system
  if (settings.printTotal) {
    const double absTotalDipole = totalDipole.norm();
    /*
     * Note: Serenity is compatible with the program SNF, which
     * searches for the string 'Total Dipole Moment' in
     * Serenity's output. Thus, please make sure that at least this
     * substring remains intact!
     */
    printf("\n");
    printSmallCaption("Total Dipole Moment (origin corrected)");
    printf("%11s x %11s y %11s z\n", "", "", "");
    printf("%4s %+9.10f %+9.10f %+9.10f a.u.\n", "", totalDipole[0], totalDipole[1], totalDipole[2]);
    printf("%4s %+9.10f %+9.10f %+9.10f Debye\n\n", "", totalDipole[0] * AU_TO_DEBYE, totalDipole[1] * AU_TO_DEBYE,
           totalDipole[2] * AU_TO_DEBYE);
    printf("%4s |Dipole Moment|: %+9.10f a.u.   %+9.10f Debye\n", "", absTotalDipole, absTotalDipole * AU_TO_DEBYE);
    if (settings.highestOrder > 1) {
      const double trace = (totalQuadrupole(0, 0) + totalQuadrupole(1, 1) + totalQuadrupole(1, 2)) / 3.0;
      printf("\n\n");
      printSmallCaption("Total Quadrupole Moment (origin corrected)");
      printf("%11s x %11s y %11s z\n", "", "", "");
      printf("%2sx%1s %+9.10f %+9.10f %+9.10f a.u.\n", "", "", totalQuadrupole(0, 0), totalQuadrupole(1, 0),
             totalQuadrupole(0, 1));
      printf("%2sy%1s %+9.10f %+9.10f %+9.10f a.u.\n", "", "", totalQuadrupole(1, 0), totalQuadrupole(1, 1),
             totalQuadrupole(0, 2));
      printf("%2sz%1s %+9.10f %+9.10f %+9.10f a.u.\n\n", "", "", totalQuadrupole(0, 1), totalQuadrupole(0, 2),
             totalQuadrupole(1, 2));
      printf("%4s 1/3 Trace: %+9.10f a.u.\n", "", trace);
    }
  }
}

void MultipoleMomentTask::calculateMultipoleMoments(Eigen::Vector3d& dipoleMoment, Eigen::MatrixXd& quadrupoleMoment,
                                                    std::shared_ptr<SystemController> system) {
  auto molCharge = system->getCharge();

  // Always calculate the analytical multipole moments, since we do not have a numerical implementation
  // for orders higher than 1 and it hardly takes any time.
  const auto lastSCFMode = system->getSCFMode();
  std::vector<std::vector<double>> multipoleMoment;
  if (lastSCFMode == RESTRICTED) {
    multipoleMoment =
        MultipoleMomentCalculator::calculateMultipoleMoment<Options::SCF_MODES::RESTRICTED>(system, settings.highestOrder);
  }
  else {
    multipoleMoment =
        MultipoleMomentCalculator::calculateMultipoleMoment<Options::SCF_MODES::UNRESTRICTED>(system, settings.highestOrder);
  }
  for (unsigned int i = 0; i < 3; i++) {
    dipoleMoment(i) = multipoleMoment[0][i];
  }
  if (settings.numerical) {
    std::cout << "\n Numerical Dipole Moment Requested!\n";
    if (lastSCFMode == RESTRICTED) {
      dipoleMoment = NumericalDipoleMomentCalculator::calculateDipoleMoment<Options::SCF_MODES::RESTRICTED>(system);
    }
    else {
      dipoleMoment = NumericalDipoleMomentCalculator::calculateDipoleMoment<Options::SCF_MODES::UNRESTRICTED>(system);
    }
  }

  /*
   * Calculate absolute Value of Dipole Moment in Debye
   */
  double absDipole = dipoleMoment.norm();
  /*
   * Output
   */
  if (settings.printFragments || _activeSystems.size() == 1) {
    printf("\n");
    printSmallCaption("Dipole Moment");
    printf("%11s x %11s y %11s z\n", "", "", "");
    printf("%4s %+9.10f %+9.10f %+9.10f a.u.\n", "", dipoleMoment[0], dipoleMoment[1], dipoleMoment[2]);
    printf("%4s %+9.10f %+9.10f %+9.10f Debye\n\n", "", dipoleMoment[0] * AU_TO_DEBYE, dipoleMoment[1] * AU_TO_DEBYE,
           dipoleMoment[2] * AU_TO_DEBYE);
    printf("%4s |Dipole Moment|: %+9.10f a.u.   %+9.10f Debye\n", "", absDipole, absDipole * AU_TO_DEBYE);
  }

  if (molCharge != 0) {
    printf("\n\n");
    printSmallCaption("Calculate origin-corrected Dipole Moment mu_corr=mu-charge*R_0");

    Eigen::Vector3d corrDipoles = Eigen::Vector3d::Zero(3);

    corrDipoles[0] = dipoleMoment[0] - molCharge * _origin.getX();
    corrDipoles[1] = dipoleMoment[1] - molCharge * _origin.getY();
    corrDipoles[2] = dipoleMoment[2] - molCharge * _origin.getZ();
    dipoleMoment = corrDipoles;

    /*
     * Calculate absolute Value of Dipole Moment in Debye
     */
    absDipole = corrDipoles.norm();
    /*
     * Output
     */
    if (settings.printFragments || _activeSystems.size() == 1) {
      printf("\n");
      printSmallCaption("Origin-Corrected Dipole Moment");
      printf("%11s x %11s y %11s z\n", "", "", "");
      printf("%4s %+9.10f %+9.10f %+9.10f a.u.\n", "", corrDipoles[0], corrDipoles[1], corrDipoles[2]);
      printf("%4s %+9.10f %+9.10f %+9.10f Debye\n\n", "", corrDipoles[0] * AU_TO_DEBYE, corrDipoles[1] * AU_TO_DEBYE,
             corrDipoles[2] * AU_TO_DEBYE);
      printf("%4s |Dipole Moment|: %+9.10f a.u.   %+9.10f Debye\n", "", absDipole, absDipole * AU_TO_DEBYE);
    }
  }

  if (settings.highestOrder > 1) {
    /*
     * Calculate 1/3 of trace
     */
    quadrupoleMoment(0, 0) = multipoleMoment[1][0]; // xx
    quadrupoleMoment(1, 0) = multipoleMoment[1][1]; // xy
    quadrupoleMoment(0, 1) = multipoleMoment[1][2]; // xz
    quadrupoleMoment(1, 1) = multipoleMoment[1][3]; // yy
    quadrupoleMoment(0, 2) = multipoleMoment[1][4]; // yz
    quadrupoleMoment(1, 2) = multipoleMoment[1][5]; // zz
    double trace = (multipoleMoment[1][0] + multipoleMoment[1][3] + multipoleMoment[1][5]) / 3;

    /*
     * Output
     */
    if (settings.printFragments || _activeSystems.size() == 1) {
      printf("\n\n");
      printSmallCaption("Quadrupole Moment");
      printf("%11s x %11s y %11s z\n", "", "", "");
      printf("%2sx%1s %+9.10f %+9.10f %+9.10f a.u.\n", "", "", multipoleMoment[1][0], multipoleMoment[1][1],
             multipoleMoment[1][2]);
      printf("%2sy%1s %+9.10f %+9.10f %+9.10f a.u.\n", "", "", multipoleMoment[1][1], multipoleMoment[1][3],
             multipoleMoment[1][4]);
      printf("%2sz%1s %+9.10f %+9.10f %+9.10f a.u.\n\n", "", "", multipoleMoment[1][2], multipoleMoment[1][4],
             multipoleMoment[1][5]);
      printf("%4s 1/3 Trace: %+9.10f a.u.\n", "", trace);
    }

    if (molCharge != 0 or absDipole != 0) {
      printf("\n\n");
      printSmallCaption("Calculate origin-corrected Quadrupole Moment q_corr=q-mu*R_0-R_0*mu+charge*R_0^2");

      std::vector<double> corrQuadrupoles(6, 0.0);

      corrQuadrupoles[0] = multipoleMoment[1][0] - _origin.getX() * multipoleMoment[0][0] -
                           _origin.getX() * multipoleMoment[0][0] + molCharge * _origin.getX() * _origin.getX();
      corrQuadrupoles[1] = multipoleMoment[1][1] - _origin.getX() * multipoleMoment[0][1] -
                           _origin.getY() * multipoleMoment[0][0] + molCharge * _origin.getX() * _origin.getY();
      corrQuadrupoles[2] = multipoleMoment[1][2] - _origin.getX() * multipoleMoment[0][2] -
                           _origin.getZ() * multipoleMoment[0][0] + molCharge * _origin.getX() * _origin.getZ();
      corrQuadrupoles[3] = multipoleMoment[1][3] - _origin.getY() * multipoleMoment[0][1] -
                           _origin.getY() * multipoleMoment[0][1] + molCharge * _origin.getY() * _origin.getY();
      corrQuadrupoles[4] = multipoleMoment[1][4] - _origin.getY() * multipoleMoment[0][2] -
                           _origin.getZ() * multipoleMoment[0][1] + molCharge * _origin.getY() * _origin.getZ();
      corrQuadrupoles[5] = multipoleMoment[1][5] - _origin.getZ() * multipoleMoment[0][2] -
                           _origin.getZ() * multipoleMoment[0][2] + molCharge * _origin.getZ() * _origin.getZ();
      quadrupoleMoment(0, 0) = corrQuadrupoles[0]; // xx
      quadrupoleMoment(1, 0) = corrQuadrupoles[1]; // xy
      quadrupoleMoment(0, 1) = corrQuadrupoles[2]; // xz
      quadrupoleMoment(1, 1) = corrQuadrupoles[3]; // yy
      quadrupoleMoment(0, 2) = corrQuadrupoles[4]; // yz
      quadrupoleMoment(1, 2) = corrQuadrupoles[5]; // zz
      /*
       * Calculate 1/3 of trace
       */
      double trace = (corrQuadrupoles[0] + corrQuadrupoles[3] + corrQuadrupoles[5]) / 3;

      /*
       * Output
       */
      if (settings.printFragments || _activeSystems.size() == 1) {
        printf("\n");
        printSmallCaption("Origin-Corrected Quadrupole Moment");
        printf("%11s x %11s y %11s z\n", "", "", "");
        printf("%2sx%1s %+9.10f %+9.10f %+9.10f a.u.\n", "", "", corrQuadrupoles[0], corrQuadrupoles[1], corrQuadrupoles[2]);
        printf("%2sy%1s %+9.10f %+9.10f %+9.10f a.u.\n", "", "", corrQuadrupoles[1], corrQuadrupoles[3], corrQuadrupoles[4]);
        printf("%2sz%1s %+9.10f %+9.10f %+9.10f a.u.\n\n", "", "", corrQuadrupoles[2], corrQuadrupoles[4], corrQuadrupoles[5]);
        printf("%4s 1/3 Trace: %+9.10f a.u.\n", "", trace);
      }
    }
  }
}

} /* namespace Serenity */
