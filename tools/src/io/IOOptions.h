/**
 * @file   IOOptions.h
 *
 * @date   Sep 9, 2014
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
#ifndef IOOPTIONS_H_
#define IOOPTIONS_H_

namespace Serenity {
/**
 * @struct IOOptions IOOptions.h
 *
 * @brief Holds the options for the output and input of data.
 */
struct IOOptions {
 public:
  /*
   * boolean switches
   */
  /**
   * @brief Boolean switch whether to print the final orbital energies.
   */
  bool printFinalOrbitalEnergies = true;

  /**
   * @brief Boolean switch whether to print the geometry.
   */
  bool printGeometry = false;

  /**
   * @brief Boolean switch whether to print the information of each SCF cycle.
   */
  bool printSCFCycleInfo = true;

  /**
   * @brief Boolean switch whether to print the SCF results.
   */
  bool printSCFResults = true;

  /**
   * @brief Boolean switch whether to print all the debug info.
   */
  bool printDebugInfos = false;
  /**
   * @brief switch whether information about integration grids shall be printed upon construction
   */
  bool printGridInfo = true;
  /**
   * @brief switch whether checks for the grid accuracy shall be printed
   */
  bool gridAccuracyCheck = false;

  /*
   * integer, level options
   */
  /**
   * @brief Print level of timing statements.
   *        0: minimal output
   *        1: normal output
   *        2: detailed output
   *        3: debug (all)
   */
  unsigned int timingsPrintLevel = 1;
};

extern IOOptions iOOptions;

} // namespace Serenity
#endif /* IOOPTIONS_H_ */
