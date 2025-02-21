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
#include "dft/functionals/BasicFunctionals.h"
#include "dft/functionals/CompositeFunctionals.h"
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
  DFT() : functional(CompositeFunctionals::XCFUNCTIONALS::BP86), dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE) {
  }

 public:
  REFLECTABLE((CompositeFunctionals::XCFUNCTIONALS)functional, (Options::DFT_DISPERSION_CORRECTIONS)dispersion)
};

struct CUSTOMFUNCTIONAL {
  CUSTOMFUNCTIONAL()
    : impl(CompositeFunctionals::IMPLEMENTATIONS::LIBXC),
      // a custom functional is interpreted to be active as soon as the basicFunctionals-vector is not empty
      basicFunctionals({}),
      mixingFactors({1.0}),
      hfExchangeRatio(0.0),
      hfCorrelRatio(0.0),
      lrExchangeRatio(0.0),
      mu(0.0),
      ssScaling(1.0),
      osScaling(1.0) {
  }

 public:
  REFLECTABLE((CompositeFunctionals::IMPLEMENTATIONS)impl, (std::vector<BasicFunctionals::BASIC_FUNCTIONALS>)basicFunctionals,
              (std::vector<double>)mixingFactors, (double)hfExchangeRatio, (double)hfCorrelRatio,
              (double)lrExchangeRatio, (double)mu, (double)ssScaling, (double)osScaling)
  // This is for enabling CUSTOMFUNC blocks in tasks.
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("CUSTOMFUNC")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};

struct SCF {
  SCF()
    : initialguess(Options::INITIAL_GUESSES::ATOM_SCF),
      maxCycles(100),
      writeRestart(5),
      energyThreshold(5e-8),
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
      diisThreshold(5e-7),
      canOrthThreshold(1.0e-7),
      useADIIS(false),
      allowNotConverged(false),
      rohf(Options::ROHF_TYPES::NONE),
      suhfLambda(0.01),
      degeneracyThreshold(0.0) {
  }

 public:
  REFLECTABLE((Options::INITIAL_GUESSES)initialguess, (unsigned int)maxCycles, (unsigned int)writeRestart,
              (double)energyThreshold, (double)rmsdThreshold, (Options::DAMPING_ALGORITHMS)damping,
              (double)seriesDampingStart, (double)seriesDampingEnd, (double)seriesDampingStep,
              (double)seriesDampingInitialSteps, (double)staticDampingFactor, (double)endDampErr, (bool)useLevelshift,
              (bool)useOffDiagLevelshift, (double)minimumLevelshift, (unsigned int)diisFlush, (double)diisStartError,
              (unsigned int)diisMaxStore, (double)diisThreshold, (double)canOrthThreshold, (bool)useADIIS,
              (bool)allowNotConverged, (Options::ROHF_TYPES)rohf, (double)suhfLambda, (double)degeneracyThreshold)
};
struct BASIS {
  BASIS()
    : label("6-31GS"),
      auxJLabel("RI_J_WEIGEND"),
      auxJKLabel(""),
      auxCLabel(""),
      makeSphericalBasis(true),
      integralThreshold(0),
      integralIncrementThresholdStart(1e-8),
      integralIncrementThresholdEnd(0),
      incrementalSteps(5),
      basisLibPath(""),
      firstECP(37),
      densFitJ(Options::DENS_FITS::RI),
      densFitK(Options::DENS_FITS::NONE),
      densFitLRK(Options::DENS_FITS::NONE),
      densFitCorr(Options::DENS_FITS::RI),
      cdThreshold(1e-6),
      extendSphericalACDShells(Options::EXTEND_ACD::SIMPLE),
      intCondition(-1),
      secondCD(1e-8),
      cdOffset(1e-2) {
  }

 public:
  REFLECTABLE((std::string)label, (std::string)auxJLabel, (std::string)auxJKLabel, (std::string)auxCLabel,
              (bool)makeSphericalBasis, (double)integralThreshold, (double)integralIncrementThresholdStart,
              (double)integralIncrementThresholdEnd, (unsigned int)incrementalSteps, (std::string)basisLibPath,
              (unsigned int)firstECP, (Options::DENS_FITS)densFitJ, (Options::DENS_FITS)densFitK,
              (Options::DENS_FITS)densFitLRK, (Options::DENS_FITS)densFitCorr, (double)cdThreshold,
              (Options::EXTEND_ACD)extendSphericalACDShells, (int)intCondition, (double)secondCD, (double)cdOffset)
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

  // This is for enabling GRID blocks in tasks.
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("GRID")) {
      visit_each(*this, v);
      return true;
    }
    return false;
  }
};

struct EFIELD {
  EFIELD()
    : use(false),
      analytical(false),
      pos1({0.0, 0.0, 0.0}),
      pos2({0.0, 0.0, 1.0}),
      distance(5e1),
      nRings(50),
      radius(1e0),
      fieldStrength(1e-3),
      nameOutput("") {
  }

 public:
  REFLECTABLE((bool)use, (bool)analytical, (std::vector<double>)pos1, (std::vector<double>)pos2, (double)distance,
              (unsigned)nRings, (double)radius, (double)fieldStrength, (std::string)nameOutput)
};

struct EXTERNALCHARGES {
  EXTERNALCHARGES() : externalChargesFile("") {
  }

 public:
  REFLECTABLE((std::string)externalChargesFile)
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
      ignoreCharge(false),
      spin(0),
      geometry(""),
      externalGridPotential(""),
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
  REFLECTABLE((std::string)name, (std::string)identifier, (std::string)path, (int)charge, (bool)ignoreCharge, (int)spin,
              (std::string)geometry, (std::string)externalGridPotential, (std::string)load, (Options::SCF_MODES)scfMode,
              (Options::ELECTRONIC_STRUCTURE_THEORIES)method)
  DFT dft;
  SCF scf;
  BASIS basis;
  GRID grid;
  EFIELD efield;
  PCMSettings pcm;
  EXTERNALCHARGES extCharges;
  CUSTOMFUNCTIONAL customFunc;

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
