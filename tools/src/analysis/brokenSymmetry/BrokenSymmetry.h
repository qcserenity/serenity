/**
 * @file   BrokenSymmetry.h
 * @author Anja Massolle
 *
 * @date   17. November 2020
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

#ifndef BROKENSYMMETRY_H_
#define BROKENSYMMETRY_H_

/* Include Serenity Internal Headers */
#include "settings/EmbeddingSettings.h"
#include "settings/LocalizationOptions.h"
#include "settings/OrthogonalizationOptions.h"
/* Include Std and External Headers */
#include <limits>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward Declarations */
class SystemController;
struct FreezeAndThawTaskSettings;

/**
 * @class BrokenSymmetry BrokenSymmetry.h
 * @brief A class which handles Broken-Symmetry calculations.
 */
class BrokenSymmetry {
 public:
  /**
   * @brief Constructor of the Broken-Symmetry Class.
   *
   * @param hsSystemController The system controller of the high spin state
   * @param embedding The embedding settings
   * @param tsOrtho if true the non-additive kinetic energy is evaluated via orthogonalized supersystem orbitals
   * @param allOrtho if true all energy contributions are evaluated via orthogonalized supersystem orbitals
   * @param orthogonalizationScheme the orthogonalization scheme used for the orthogonalization
   * @param bsSystemController if specified these systems are loaded from disk as broken-symmetry systems. E.g.
   *                           for a different energy evaluation.
   */
  BrokenSymmetry(std::vector<std::shared_ptr<SystemController>> hsSystemController, EmbeddingSettings embedding = {},
                 bool tsOrtho = false, bool allOrtho = false,
                 Options::ORTHOGONALIZATION_ALGORITHMS orthogonalizationScheme = Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN,
                 std::vector<std::shared_ptr<SystemController>> bsSystemController = {});
  /**
   * @brief Default destructor.
   */
  virtual ~BrokenSymmetry() = default;
  /**
   * @brief performs a standard BS-DFT calculation
   * @param nA number of unpaired electrons of spin site 1
   * @param nB number of unpaired electrons of spin site 2
   * @param locType the localization procedure applied for localizing the SONOs
   * @param noThreshold Threshold for the assignment of the NO orbitals as SONO, UONO and DONO
   */
  void bsDFT(unsigned int nA, unsigned int nB,
             Options::ORBITAL_LOCALIZATION_ALGORITHMS locType = Options::ORBITAL_LOCALIZATION_ALGORITHMS::PIPEK_MEZEY,
             double noThreshold = 0.2);
  /**
   * @brief performs a parallel FDE calculation for the HS and BS state
   */
  void bsFDE();
  /**
   * @brief performs a FAT calculation for the HS and BS state
   * @param fatSettings Settings for the freeze and thaw task
   */
  void bsFAT(FreezeAndThawTaskSettings fatSettings);
  /**
   * @brief use the isolated (spin) densities for the energy evaluation of the HS and BS state
   */
  void bsIsolated();

  /**
   * @brief Getter for the density based S2 value of the HS state
   */
  double getS2HS();

  /**
   * @brief Getter for the density based S2 value of the BS state
   */
  double getS2BS();

  /**
   * @brief Getter for the overlap of the magnetic orbitals
   */
  double getSab();
  /**
   * @brief Getter for the orbital based S2 value of the HS state
   */
  double getS2uhfHS();
  /**
   * @brief Getter for the orbital based S2 value of the BS state
   */
  double getS2uhfBS();
  /**
   * @brief Getter for the magnetic exchange coupling constant J1, see manual
   */
  double getJ1();
  /**
   * @brief Getter for the magnetic exchange coupling constant J2, see manual
   */
  double getJ2();
  /**
   * @brief Getter for the magnetic exchange coupling constant J3, see manual
   */
  double getJ3();
  /**
   * @brief Getter for the magnetic exchange coupling constant J3 where the S2 values are evaluated based on the
   * orbitals, see manual
   */
  double getJ3UHF();
  /**
   * @brief Getter for the magnetic exchange coupling constant J4, see manual
   */
  double getJ4();

 private:
  // The system Controller holding the HS state of the spin sites
  std::vector<std::shared_ptr<SystemController>> _hsSystemController;
  // The embedding Settings for sDFT calculations
  EmbeddingSettings _embedding;
  // if true the non-additive kinetic energy is evaluated from orthogonal supersystem orbitals
  bool _tsOrtho;
  // if true all energy contributions are evaluated from orthogonal supersystem orbitals
  bool _allOrtho;
  // The orthogonalization scheme used for tsOrtho or allOrtho
  Options::ORTHOGONALIZATION_ALGORITHMS _orthogonalizationScheme;
  // The system Controller holding the BS states of the spin sites (if loaded from disk)
  std::vector<std::shared_ptr<SystemController>> _bsSystemController;
  // S*S of the HS state evaluated with the spin density
  double _s2HS = std::numeric_limits<double>::infinity();
  // S*S of the BS state evaluated with the spin density
  double _s2BS = std::numeric_limits<double>::infinity();
  // Overlap of the magnetic orbitals after a corresponding orbital transformation
  double _sAB = std::numeric_limits<double>::infinity();
  // S*S of the HS state evaluated with the orbitals
  double _s2UHFhs = std::numeric_limits<double>::infinity();
  // S*S of the BS state evaluated with the orbitals
  double _s2UHFbs = std::numeric_limits<double>::infinity();
  // J1
  double _j1 = std::numeric_limits<double>::infinity();
  // J2
  double _j2 = std::numeric_limits<double>::infinity();
  // J3 where S*S is evaluated based on the spin density
  double _j3 = std::numeric_limits<double>::infinity();
  // J3 where S*S is evaluated based on orbitals
  double _j3UHF = std::numeric_limits<double>::infinity();
  // J4
  double _j4 = std::numeric_limits<double>::infinity();
  // evaluates the non-additive kinetic energy from orthogonal supersystem orbitals
  void tsOrtho();
  // evaluates all energy contributions from orthogonal supersystem orbitals
  void allOrtho();
  // calculates the magnetic exchange coupling constants
  void evalJ();
};
} /*namespace Serenity*/

#endif /* TASKS_BROKENSYMMETRY_H_ */