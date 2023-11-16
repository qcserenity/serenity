/**
 * @file ELFCalculator.h
 *
 * @date Apr 4, 2016, rework Jun 28, 2017
 * @author Melanie Börner, rework Jan Unsleber
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

#ifndef POSTSCF_LOCALIZATIONFUNCTIONS_ELFCALCULATOR_H_
#define POSTSCF_LOCALIZATIONFUNCTIONS_ELFCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
#include "tasks/Task.h"

namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class GridData;
class SystemController;
class GridController;
class AtomCenteredBasisController;

/**
 * @class ELFCalculator ELFCalculator.h
 * @brief Calculates the kinetic energy density and the Electron Localization Function (ELF).
 *        To calculate the ELF according to Savin et al., use "printELFOnGrid".
 *        To calculate the ELF according to Tsirelson and Stash, use "print ELF_TSOnGrid".
 * @param systemController which provides you with all needed information about your system and
 *        configuration
 */
template<Options::SCF_MODES SCFMode>
class ELFCalculator {
 public:
  ELFCalculator(std::shared_ptr<SystemController> systemController);
  virtual ~ELFCalculator() = default;

  /**
   * @brief Calculates the kinetic energy density on a grid.
   *        This function must be given the gridController of the respective grid on which the kinetic energy density is
   * to be calculated.
   * @returns A vector of doubles holding the kinetic energy density on each grid point.
   */
  GridData<SCFMode> calcKineticEnergyDensity(std::shared_ptr<GridController> gridController);

  /**
   * @brief Calculates the ELF on a grid.
   *
   * The Electron Localization Function (ELF) is defined as:
   *
   * \f$ ELF(r) = \left[1+\left(\frac{\tau-\frac{1}{8}\frac{(\nabla \rho(r))^2}
   * {\rho(r)}}{\frac{3}{10}\left(3\pi^2\right)^{2/3}\rho(r)^{5/3}}\right)^2\right]^{-1}\f$
   *
   * where \f$ \rho \f$ is the electron density and the kinetic energy density \f$ \tau \f$ is:
   *
   * \f$ \tau(r) = \frac{1}{2} \sum \limits_{i=1}^n |\nabla \psi_i(r)|^2\f$
   *
   * where n is the number of electrons and \f$ \psi \f$ are the molecular orbitals.
   *
   * Ref.: A. Savin, O. Jepsen, J. Flad, O. K. Anderson, H. Preuß,
   * H. G. von Schnering. Angew. Chem. 1992, 104, 186.
   *
   * If the SCF mode is unrestricted, a SCFMode-polarized ELF for alpha and beta electrons each will be
   * returned. Note that ELF_alpha + ELF_beta does not equal the total ELF.
   * In order to calculate a total ELF for unrestricted systems, use calculateTotalELFOnGrid().
   *
   * @param gridController The grid to work on.
   * @returns The ELF on the grid for the total restricted density or both alpha and beta density separate.
   */
  GridData<SCFMode> calculateELFOnGrid(std::shared_ptr<GridController> gridController);

  /**
   * @brief Calculates the total ELF (unrestricted).
   *
   * \f$ ELF_{total} = \left[1+\left(\frac{\tau_\alpha + \tau_\beta -\frac{1}{8}\frac{(\nabla \rho_\alpha(r))^2}
   * {\rho_\alpha(r)}-\frac{1}{8}\frac{(\nabla \rho_\beta(r))^2}
   * {\rho_\beta(r)}}{\frac{3}{10}\left(6\pi^2\right)^{2/3}\left(\rho_\alpha(r)^{5/3}+\rho_\beta(r)^{5/3}\right)}\right)^2\right]^{-1}\f$
   *
   * Ref.: M. Kohout and A. Savin. Int. J. Quant. Chem. 1996, 60, 875.
   *
   * @param gridController The grid to work on.
   * @returns The total ELF on the grid given grid.
   */
  GridData<RESTRICTED> calculateTotalELFOnGrid(std::shared_ptr<GridController> gridController);

  /**
   * @brief Calculates the approximated ELF.
   *
   * The approximate ELF (ELF_TS) is defined as:
   *
   * \f$ ELF_{TS}(r) = \left[1+\left(\frac{\tau_{TS}-\frac{1}{8}\frac{(\nabla \rho(r))^2}
   * {\rho(r)}}{\frac{3}{10}\left(3\pi^2\right)^{2/3}\rho(r)^{5/3}}\right)^2\right]^{-1}\f$
   *
   * and the kinetic energy density \f$ \tau_{TS} \f$ is:
   *
   * \f$ \tau_{TS}(r) = \frac{3}{10}\left(3\pi^2\right)^{2/3}\rho(r)^{5/3} + \frac{1}{72}
   * \frac{|\nabla \rho(r)|^2}{\rho(r)} +  \frac{1}{6} \nabla^2\rho(r)\f$
   *
   * Ref.:  V. Tsirelson and A. Stash. Determination of the electron localization function from
   * electron density. Chem. Phys. Lett., 351, 142-148, 2002.
   *
   * If the SCF mode is unrestricted, a SCFMode-polarized ELF_TS for alpha and beta electrons each will be
   * returned. Note that ELF_TS_alpha + ELF_TS_beta does not equal the total ELF_TS.
   * In order to calculate a total ELF for unrestricted systems, use calculateTotalELFTSOnGrid().
   *
   * @param gridController The grid to work on.
   * @returns The ELF_TS on the grid for the total restricted density or both alpha and beta density separate.
   */
  GridData<SCFMode> calculateELFTSOnGrid(std::shared_ptr<GridController> gridController);

  /**
   * @brief Calculates the total approximated ELF (unrestricted).
   *
   * See "printTotalELFOnGrid()" and "printELF_TSOnGrid()" for further information.
   *
   * Ref.:  V. Tsirelson and A. Stash. Determination of the electron localization function from
   *        electron density. Chem. Phys. Lett., 351, 142-148, 2002.
   *
   * Ref.: M. Kohout and A. Savin. Int. J. Quant. Chem. 1996, 60, 875.
   *
   * @param gridController The grid to work on.
   * @returns The total ELF_TS on the grid given grid.
   */
  GridData<RESTRICTED> calculateTotalELFTSOnGrid(std::shared_ptr<GridController> gridController);

 private:
  std::shared_ptr<SystemController> _systemController;
  std::shared_ptr<AtomCenteredBasisController> _basisController;
};
} /* namespace Serenity */
#endif /* POSTSCF_LOCALIZATIONFUNCTIONS_ELFCALCULATOR_H_ */
