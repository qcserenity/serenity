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

namespace Serenity {

class PCMSettings;
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
};

} /* namespace Serenity */

#endif /* SOLVATION_SOLVENTS_H_ */
