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
#include "settings/LocalizationOptions.h"
#include "settings/PCMSettings.h"
#include "settings/Reflection.h"
#include "settings/Settings.h"
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
 * projector should be truncated (default: false). See HuzinagaProjectionPotential.h. -projectionTruncThreshold: The
 * projection truncation threshold (default: 1.0e+1). See HuzinagaProjectionPotential.h.
 * -projecTruncThresh : Total overlap threshold for the truncation of the projection operator. Only used if
truncateProjector true and hence only useful in any kind of projection technique. By
default 1.0e+1.
 * -fermiShift : An optional shift for the Huzinaga operator. Only used if
embeddingMode FERMI_SHIFTED_HUZINAGA is used. By default 1.0.
 * -calculateMP2Correction : A flag to control the evaluation of the MP2 correction for double hybrid functionals
(default: true).
 * -fullMP2Coupling : If true, the MP2 contribution of the non-additive exchange-correlation energy for double hy-
brids captures the effect of environment orbital pairs on the active-pair amplitudes (default: false).
 * -naddXCFuncList : A list of
 * non-additive exchange correlation functionals -naddKinFuncList : A list of non-additive kinetic functionals
 *        -embeddingModeList : A list of embedding types to use (e.g., level shift, potential reconstruction... see
 * Options class)
 * -partialChargesForCoulombInt : Use partial charges to approximate the Coulomb integral as a sum.
 * -chargeModel : Partial Charge model to be used for the calculation of the partial charges in the approximation of the
Coulomb integral.
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
      calculateMP2Correction(true),
      fullMP2Coupling(false),
      naddXCFuncList({}),
      naddKinFuncList({}),
      embeddingModeList({}),
      partialChargesForCoulombInt(false),
      chargeModel(Options::POPULATION_ANALYSIS_ALGORITHMS::MULLIKEN),
      loewdinOrder(1),
      loewdinWeights({1.0, 1.0, 1.0}) {
  }
  REFLECTABLE((double)levelShiftParameter, (CompositeFunctionals::XCFUNCTIONALS)naddXCFunc,
              (CompositeFunctionals::KINFUNCTIONALS)naddKinFunc, (CompositeFunctionals::KINFUNCTIONALS)longRangeNaddKinFunc,
              (Options::KIN_EMBEDDING_MODES)embeddingMode, (Options::DFT_DISPERSION_CORRECTIONS)dispersion,
              (double)smoothFactor, (std::string)potentialBasis, (double)singValThreshold, (double)lbDamping,
              (unsigned int)lbCycles, (unsigned int)carterCycles, (double)borderAtomThreshold, (double)basisFunctionRatio,
              (bool)truncateProjector, (double)projecTruncThresh, (double)fermiShift, (bool)calculateMP2Correction,
              (bool)fullMP2Coupling, (std::vector<CompositeFunctionals::XCFUNCTIONALS>)naddXCFuncList,
              (std::vector<CompositeFunctionals::KINFUNCTIONALS>)naddKinFuncList,
              (std::vector<Options::KIN_EMBEDDING_MODES>)embeddingModeList, (bool)partialChargesForCoulombInt,
              (Options::POPULATION_ANALYSIS_ALGORITHMS)chargeModel, (unsigned int)loewdinOrder,
              (std::vector<double>)loewdinWeights)
 public:
  PCMSettings pcm;
  CUSTOMFUNCTIONAL customNaddXCFunc;
  CUSTOMFUNCTIONAL customNaddKinFunc;
  CUSTOMFUNCTIONAL customLongRangeNaddKinFunc;
  /**
   * @brief Parse the settings from the visitor to this object.
   * @param v The visitor.
   * @param blockname The block name.
   * @return True if the block name corresponds to this block. Otherwise false.
   */
  bool visitAsBlockSettings(set_visitor v, std::string blockname) {
    if (!blockname.compare("EMB")) {
      visit_each(*this, v);
      return true;
    }
    else if (this->pcm.visitAsBlockSettings(v, blockname)) {
      return true;
    }
    // visit_each instead of visitAsBlockSettings because these are custom functional blocks with different block names
    else if (!blockname.compare("CUSTOMNADDXC")) {
      visit_each(this->customNaddXCFunc, v);
      return true;
    }
    else if (!blockname.compare("CUSTOMNADDKIN")) {
      visit_each(this->customNaddKinFunc, v);
      return true;
    }
    else if (!blockname.compare("CUSTOMLONGRANGENADDKIN")) {
      visit_each(this->customLongRangeNaddKinFunc, v);
      return true;
    }
    return false;
  }
};

} /* namespace Serenity */

#endif /* SETTINGS_EMBEDDINGSETTINGS_H_ */
