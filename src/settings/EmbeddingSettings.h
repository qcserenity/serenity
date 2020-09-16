/**
 * @file EmbeddingSettings.h
 *
 * @date Aug 9, 2019
 * @author Moritz Bensberg
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

#ifndef SETTINGS_EMBEDDINGSETTINGS_H_
#define SETTINGS_EMBEDDINGSETTINGS_H_
/* Include Serenity Internal Headers */
#include "settings/DFTOptions.h"
#include "settings/EmbeddingOptions.h"
#include "settings/PCMSettings.h"
#include "settings/Reflection.h"
/* Include Std and External Headers */
#include <string>

namespace Serenity {
using namespace Serenity::Reflection;

/**
 * @brief Settings commonly used in embedding calculations.
 *        -levelShiftParameter : Shift projected orbitals by levelShiftParameters (default: 1.0e6 Eh)
 *        -naddXCFunc : The non-additive exchange correlation functional
 *        -naddKinFunc : The non-additive kinetic functional (used in for EMBEDDING_MODE::NADD_FUNC)
 *        -longRangeNaddKinFunc : The kinetic energy functional used to correct contributions from non orthogonal
 * orbitals. -embeddingMode : The type of embedding to run (e.g., level shift, potential reconstruction... see Options
 * class) -dispersion : The dispersion interaction between subsystems. -smoothFactor : Smoothing to be used in potential
 * reconstruction -potentialBasis : Basis to express the potential in during Wu-Yang reconstruction -singValThreshold :
 * Threshold for singular value decomposition in Wu-Yang Newton-Raphson step -lbDamping : Damping to be used during the
 * van Leeuwen-Baerends reconstruction -lbCycles : Maximum cycles for van Leeuwen-Baerends scheme -carterCycles :
 * Maximum cycles for Zhang-Carter scheme -borderAtomThreshold : The Mulliken population threshold used to determine if
 * an orbital is considered "distant" or not. The Mulliken population of the orbital on all not "distant" atoms has to
 * exceed this threshold in order to be included in the projector. -basisFunctionRatio : The minimum ratio of retained
 * basis functions needed in order to consider an atom to be not "distant". -truncateProjector: A flag whether the
 * projector should be truncated (default: false). See HuzinagaFDEProjectionPotential.h. -projectionTruncThreshold: The
 * projection truncation threshold (default: 1.0e+1). See HuzinagaFDEProjectionPotential.h.
 */
struct EmbeddingSettings {
  EmbeddingSettings()
    : levelShiftParameter(1e6),
      naddXCFunc(CompositeFunctionals::XCFUNCTIONALS::BP86),
      naddKinFunc(CompositeFunctionals::KINFUNCTIONALS::PW91K),
      longRangeNaddKinFunc(CompositeFunctionals::KINFUNCTIONALS::NONE),
      embeddingMode(Options::KIN_EMBEDDING_MODES::LEVELSHIFT),
      dispersion(Options::DFT_DISPERSION_CORRECTIONS::NONE),
      smoothFactor(0.0),
      potentialBasis(""),
      singValThreshold(0.0),
      lbDamping(0.995),
      lbCycles(0),
      carterCycles(0),
      borderAtomThreshold(0.02),
      basisFunctionRatio(0.0),
      truncateProjector(false),
      projecTruncThresh(1.0e+1),
      fermiShift(1.0),
      calculateMP2Correction(true) {
  }
  REFLECTABLE((double)levelShiftParameter, (CompositeFunctionals::XCFUNCTIONALS)naddXCFunc,
              (CompositeFunctionals::KINFUNCTIONALS)naddKinFunc, (CompositeFunctionals::KINFUNCTIONALS)longRangeNaddKinFunc,
              (Options::KIN_EMBEDDING_MODES)embeddingMode, (Options::DFT_DISPERSION_CORRECTIONS)dispersion,
              (double)smoothFactor, (std::string)potentialBasis, (double)singValThreshold, (double)lbDamping,
              (unsigned int)lbCycles, (unsigned int)carterCycles, (double)borderAtomThreshold, (double)basisFunctionRatio,
              (bool)truncateProjector, (double)projecTruncThresh, (double)fermiShift, (bool)calculateMP2Correction)
 public:
  PCMSettings pcm;
  /**
   * @brief Parse the settings from the visitor to this object.
   * @param v The visitor.
   * @param blockname The block name.
   * @return True if the block name corresponds to this block. Otherwise false.
   */
  bool visitSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("EMB")) {
      visit_each(*this, v);
      return true;
    }
    else if (!blockname.compare("PCM")) {
      visit_each(this->pcm, v);
      return true;
    }
    return false;
  }
};

} /* namespace Serenity */

#endif /* SETTINGS_EMBEDDINGSETTINGS_H_ */
