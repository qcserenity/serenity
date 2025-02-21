/**
 * @file Solvents.h
 *
 * @author Moritz Bensberg
 * @date May 19, 2020
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

#ifndef SOLVATION_SOLVENTS_H_
#define SOLVATION_SOLVENTS_H_

/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {

class AtomType;
struct PCMSettings;
namespace Options {
enum class PCM_SOLVENTS;
}

/**
 * @class Solvents Solvents.h
 * @brief Contains tabulated data for different implicit solvents.
 */
class Solvents {
 private:
  Solvents() = default;
  virtual ~Solvents() = default;

 public:
  /**
   * @brief Getter for the probe radius.
   * @param solvent The solvent type. Has to be != EXPLICIT.
   *
   * Parameters taken from:
   *   https://pcmsolver.readthedocs.io/en/stable/users/input.html#input-parameters
   * @return The probe radius.
   */
  static double getProbeRadius(Options::PCM_SOLVENTS solvent);
  /**
   * @brief Getter for the probe radius of the solvent to be used in the calculation
   *        of the cavity formation energy.
   * @param solvent The solvent.
   * @return The probe radius.
   */
  static double getCavityFormationProbeRadius(Options::PCM_SOLVENTS solvent);
  /**
   * @brief Getter for the static permittivity.
   * @param solvent The solvent type. Has to be != EXPLICIT.
   *
   * Parameters taken from:
   *   https://pcmsolver.readthedocs.io/en/stable/users/input.html#input-parameters
   * @return The static permittivity.
   */
  static double getStaticPermittivity(Options::PCM_SOLVENTS solvent);
  /**
   * @brief Print information about the implicit solvation model:
   *          Solvation model\n
   *          Solvent\n
   *          Static permittivity\n
   *          Probe radius\n
   * @param settings
   */
  static void printSolventInfo(const PCMSettings& settings);
  /**
   * @brief Returns the number density (number of particles per volume) of the solvent.
   *        at 25 Celsius.
   * @param solvent The solvent.
   * @return The number density.
   */
  static double getNumberDensity(Options::PCM_SOLVENTS solvent);
  /**
   * @brief Getter for the molar volume of the solvent (volume per mole)
   *        at 25 Celsius.
   * @param solvent The solvent.
   * @return The molar volume.
   */
  static double getMolarVolume(Options::PCM_SOLVENTS solvent);
  /**
   * @brief Getter for the molecular weight of the solvent.
   * @param solvent The solvent.
   * @return The molecular weight.
   */
  static double getMolecularWeight(Options::PCM_SOLVENTS solvent);
  /**
   * @brief Getter for the density of the solvent at 25 Celsius.
   * @param solvent The solvent.
   * @return The density.
   */
  static double getDensity(Options::PCM_SOLVENTS solvent);

 private:
  static double resolveMass(std::vector<std::pair<unsigned int, std::shared_ptr<const AtomType>>> types);
  static std::vector<std::pair<unsigned int, std::shared_ptr<const AtomType>>> getAtomTypes(Options::PCM_SOLVENTS solvent);
};

} /* namespace Serenity */

#endif /* SOLVATION_SOLVENTS_H_ */
