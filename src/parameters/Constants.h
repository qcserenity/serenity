/**
 * @file   Constants.h
 *
 * @date   long before Feb 2014
 * @author Thomas Dresselhaus
 *
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
#ifndef CONSTANTS_H
#define CONSTANTS_H

/* Include Std and External Headers */
#include <array>

namespace Serenity {
/**
 * The maximum supported angular momentum of basis functions
 */
constexpr unsigned int AM_MAX = 6;
/**
 * Number of basis functions in a shell with angular momentum [index]
 * (Cartesian)
 */
constexpr std::array<unsigned int, 11> N_SHELL_CART = {{1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66}};
/**
 * Number of basis functions in a shell with angular momentum [index]
 * (Spherical)
 */
constexpr std::array<unsigned int, 11> N_SHELL_SPH = {{1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21}};
/**
 * basis function labels for given angular momentum
 */
constexpr std::array<char, 20> ANGMOM_TO_LABEL = {
    {'s', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z'}};
/**
 * @enum ANGULAR_QUANTUM_NUMBER
 * @brief also called azimuthal quantum number
 */
enum class ANGULAR_QUANTUM_NUMBER { s = 0, p = 1, d = 2, f = 3, g = 4, h = 5, i = 6, k = 7, m = 8, n = 9, o = 10 };
/***************************************************************************
 * The data below are taken from:                                          *
 * Peter J. Mohr, David B. Newell, Barry N. Taylor, Eite Tiesinga          *
 * "CODATA Recommended Values of the Fundamental Physical Constants: 2022" *
 * https://physics.nist.gov/constants (08.08.2024)                         *
 ***************************************************************************/
/*
 * Natural Constants
 */
// The number of particles in one mol. Unit: mol^-1
constexpr double AVOGADRO_CONSTANT = 6.02214076E+023;
constexpr double BOHR = 5.29177210544E-011;
constexpr double SPEEDOFLIGHT = 299792458.0;
constexpr double PI = 3.14159265358979;
constexpr double BOLTZMANN_CONSTANT = 1.380649E-23;
constexpr double ELEMENTARY_CHARGE = 1.602176634E-019;
constexpr double PLANCK_CONSTANT = 6.62607015E-34;
constexpr double UNIVERSALGAS_CONSTANT = 8.314462618;
constexpr double VAC_ELEC_PERMITTIVITY = 8.8541878188E-12;
constexpr double SPEEDOFLIGHT_AU =
    2 * VAC_ELEC_PERMITTIVITY * SPEEDOFLIGHT * PLANCK_CONSTANT / (ELEMENTARY_CHARGE * ELEMENTARY_CHARGE);
constexpr double ELEC_MASS = 9.1093837139E-31; // Electron mass
/*
 * Conversion factors
 */
/* Length */
///  inBohr * BOHR_TO_ANGSTROM = inAngstrom
constexpr double BOHR_TO_ANGSTROM = BOHR * 1.0E+010;
///  inAngstrom * ANGSTROM_TO_BOHR = inBohr
constexpr double ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM;

/* Energy */
/// inHARTREE * HARTREE_TO_J = inJoule
constexpr double HARTREE_TO_J = PLANCK_CONSTANT * PLANCK_CONSTANT / (4 * PI * PI * ELEC_MASS * BOHR * BOHR);
/// inHartree * HARTREE_TO_EV = inElectronVolt
constexpr double HARTREE_TO_EV = HARTREE_TO_J / ELEMENTARY_CHARGE;
/// inElectronVolt * EV_TO_HARTREE = inHartree
constexpr double EV_TO_HARTREE = 1.0 / HARTREE_TO_EV;
/// HARTREE_TO_NM / inHartree = inNanoMetre
constexpr double HARTREE_TO_NM = PLANCK_CONSTANT * SPEEDOFLIGHT / HARTREE_TO_J * 1.0E+09;
/// inHartree * HARTREE_TO_OOCM = inWaveNumber
constexpr double HARTREE_TO_OOCM = HARTREE_TO_J / (PLANCK_CONSTANT * SPEEDOFLIGHT * 100);
/// inHartree * HARTREE_TO_KJ_PER_MOL = inkJPerMol
constexpr double HARTREE_TO_KJ_PER_MOL = HARTREE_TO_J * AVOGADRO_CONSTANT / 1000;
/* The elementary charge in Coulomb */
constexpr double AU_TO_COULOMB = ELEMENTARY_CHARGE;
/* Dipole moment */
constexpr double AU_TO_COULOMB_METER = AU_TO_COULOMB * BOHR_TO_ANGSTROM * 1e-10;
/*
 * Warning: this conversion factor from the CGS unit system may be inaccurate (inherent to the
 * unit system). This factor was found on https://physics.nist.gov/cgi-bin/cuu/Value?auedm (Nov 08, 2017)
 */
/// inDebye * DEBYE_TO_AU = inAU
constexpr double DEBYE_TO_AU = 0.393430307;
/// inAU * AU_TO_DEBYE = inDebye
constexpr double AU_TO_DEBYE = 1.0 / DEBYE_TO_AU;
/// inHartree * HARTREE_TO_KCAL_PER_MOL = inKcalPerMol
constexpr double HARTREE_TO_KCAL_PER_MOL = HARTREE_TO_J * AVOGADRO_CONSTANT / 4184;

/// AU to 10^{-40 cgs} (needed for rotatory strengths)
constexpr double AU_TO_CGS = 64604.8164;

constexpr double ELEC_MASS_TO_U = ELEC_MASS * AVOGADRO_CONSTANT * 1000;
constexpr double U_TO_ELEC_MASS = 1.0 / ELEC_MASS_TO_U;

} /* namespace Serenity */
#endif /* CONSTANTS_H_ */
