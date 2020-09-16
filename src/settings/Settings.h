/**
 * @file   Settings.h
 *
 * @date   Mar 7, 2014
 * @author Thomas Dresselhaus
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
#ifndef SETTINGS_H_
#define SETTINGS_H_
/* Include Serenity Internal Headers */
#include "io/Filesystem.h"
#include "misc/SerenityError.h"
#include "settings/BasisOptions.h"
#include "settings/DFTOptions.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/GridOptions.h"
#include "settings/Options.h"
#include "settings/PCMSettings.h"
#include "settings/SCFOptions.h"
/* Include Std and External Headers */
#include <sys/stat.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace Serenity {
using namespace Serenity::Reflection;
/**
 * @class  Settings Settings.h
 * @brief  Holds one set of all switchable options for a single system.
 */
struct DFT {
  DFT()
    : functional(CompositeFunctionals::XCFUNCTIONALS::BP86),
      densityFitting(Options::DENS_FITS::RI),
      dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE) {
  }

 public:
  REFLECTABLE((CompositeFunctionals::XCFUNCTIONALS)functional, (Options::DENS_FITS)densityFitting,
              (Options::DFT_DISPERSION_CORRECTIONS)dispersion)
};

struct SCF {
  SCF()
    : initialguess(Options::INITIAL_GUESSES::ATOM_SCF),
      maxCycles(100),
      writeRestart(5),
      energyThreshold(1e-8),
      rmsdThreshold(1e-8),
      damping(Options::DAMPING_ALGORITHMS::SERIES),
      seriesDampingStart(0.7),
      seriesDampingEnd(0.2),
      seriesDampingStep(0.05),
      seriesDampingInitialSteps(2),
      staticDampingFactor(0.7),
      endDampErr(5e-2),
      useLevelshift(true),
      useOffDiagLevelshift(false),
      minimumLevelshift(0.0),
      diisFlush(1000),
      diisStartError(5e-2),
      diisMaxStore(10),
      diisThreshold(1e-7),
      canOrthThreshold(1.0e-7),
      useADIIS(false) {
  }

 public:
  REFLECTABLE((Options::INITIAL_GUESSES)initialguess, (unsigned int)maxCycles, (unsigned int)writeRestart,
              (double)energyThreshold, (double)rmsdThreshold, (Options::DAMPING_ALGORITHMS)damping,
              (double)seriesDampingStart, (double)seriesDampingEnd, (double)seriesDampingStep,
              (double)seriesDampingInitialSteps, (double)staticDampingFactor, (double)endDampErr, (bool)useLevelshift,
              (bool)useOffDiagLevelshift, (double)minimumLevelshift, (unsigned int)diisFlush, (double)diisStartError,
              (unsigned int)diisMaxStore, (double)diisThreshold, (double)canOrthThreshold, (bool)useADIIS)
};
struct BASIS {
  BASIS()
    : label("6-31GS"),
      auxJLabel("RI_J_WEIGEND"),
      auxCLabel(""),
      makeSphericalBasis(true),
      integralThreshold(1e-10),
      integralIncrementThresholdStart(1e-8),
      integralIncrementThresholdEnd(1e-12),
      incrementalSteps(5),
      basisLibPath(""),
      firstECP(37) {
  }

 public:
  REFLECTABLE((std::string)label, (std::string)auxJLabel, (std::string)auxCLabel, (bool)makeSphericalBasis,
              (double)integralThreshold, (double)integralIncrementThresholdStart, (double)integralIncrementThresholdEnd,
              (unsigned int)incrementalSteps, (std::string)basisLibPath, (unsigned int)firstECP)
};

struct GRID {
  GRID()
    : gridType(Options::GRID_TYPES::SSF),
      radialGridType(Options::RADIAL_GRID_TYPES::AHLRICHS),
      sphericalGridType(Options::SPHERICAL_GRID_TYPES::LEBEDEV),
      blocksize(128),
      accuracy(4),
      smallGridAccuracy(2),
      blockAveThreshold(1e-11),
      basFuncRadialThreshold(1e-9),
      weightThreshold(1e-14),
      smoothing(3),
      gridPointSorting(true) {
  }

 public:
  REFLECTABLE((Options::GRID_TYPES)gridType, (Options::RADIAL_GRID_TYPES)radialGridType,
              (Options::SPHERICAL_GRID_TYPES)sphericalGridType, (unsigned int)blocksize, (unsigned int)accuracy,
              (unsigned int)smallGridAccuracy, (double)blockAveThreshold, (double)basFuncRadialThreshold,
              (double)weightThreshold, (unsigned int)smoothing, (bool)gridPointSorting)
};
struct Settings {
 public:
  /**
   * @brief Default constructor
   *
   * This constructor holds all the default values for settings
   *   that are not part of any other block.
   * Each block has its defaults set in its own constructor.
   *
   * Add new (non-blocked) variables in this constructor.
   */
  Settings()
    : name("default"),
      identifier(""),
      path("./"),
      charge(0),
      spin(0),
      geometry(""),
      load(""),
      scfMode(Options::SCF_MODES::RESTRICTED),
      method(Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    // Set identifier
    identifier = std::to_string(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }
  /*
   * Member Variables
   */
  REFLECTABLE((std::string)name, (std::string)identifier, (std::string)path, (int)charge, (int)spin, (std::string)geometry,
              (std::string)load, (Options::SCF_MODES)scfMode, (Options::ELECTRONIC_STRUCTURE_THEORIES)method, (DFT)dft,
              (SCF)scf, (BASIS)basis, (GRID)grid, (PCMSettings)pcm)

  /**
   * @brief Constructor using text input.
   */
  Settings(std::ifstream& input);

  /**
   * @brief Setter for options.
   *
   * @param blockname
   * @param name
   * @param value
   */
  void set(std::string blockname, std::string name, std::string value);

  void printSettings();
};
} /* namespace Serenity */

#endif /* SETTINGS_H_ */
