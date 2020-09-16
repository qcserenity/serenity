/**
 * @file   LRSCFTask.h
 *
 * @date   Aug 17, 2016
 * @author M. Boeckers
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
#ifndef LRSCFTASK_H_
#define LRSCFTASK_H_
/* Include Serenity Internal Headers */
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/EmbeddingSettings.h"
#include "settings/LRSCFOptions.h"
#include "tasks/Task.h"

namespace Serenity {
/* Forward declarations */
class SystemController;

template<Options::SCF_MODES SCFMode>
class LRSCFController;
using namespace Serenity::Reflection;
struct LRSCFTaskSettings {
  LRSCFTaskSettings()
    : nEigen(3),
      pseudoDensityThreshold(1.0e-8),
      convThresh(1.0e-5),
      maxCycles(100),
      maxSubspaceDimension(1e+9),
      responseType(Options::RESPONSE_PROBLEM::TDA),
      tda(false),
      dominantThresh(10),
      superSystemGrid(true),
      func(CompositeFunctionals::XCFUNCTIONALS::NONE),
      analysis(true),
      setAlpha({}),
      excludeAlpha({}),
      setBeta({}),
      excludeBeta({}),
      besleyAtoms(0),
      besleyCutoff({}),
      excludeProjection(false),
      energyInclusion({}),
      energyExclusion({}),
      uncoupledSubspace({}),
      fullFDEc(false),
      loadType(Options::LRSCF_TYPE::UNCOUPLED),
      localMO(false),
      gauge(Options::GAUGE::LENGTH),
      gaugeOrigin({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()}),
      frequencies({}),
      frequencyRange({}),
      damping(0.0),
      couplingPattern({}),
      localVirtualOrbitals(0.0),
      envVirtualOrbitals(0.0),
      saveResponseMatrix(false) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::NONE;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  };
  REFLECTABLE((unsigned int)nEigen, (double)pseudoDensityThreshold, (double)convThresh, (unsigned int)maxCycles,
              (unsigned int)maxSubspaceDimension, (Options::RESPONSE_PROBLEM)responseType, (bool)tda,
              (double)dominantThresh, (bool)superSystemGrid, (CompositeFunctionals::XCFUNCTIONALS)func, (bool)analysis,
              (std::vector<unsigned int>)setAlpha, (std::vector<unsigned int>)excludeAlpha, (std::vector<unsigned int>)setBeta,
              (std::vector<unsigned int>)excludeBeta, (unsigned int)besleyAtoms, (std::vector<double>)besleyCutoff,
              (bool)excludeProjection, (std::vector<double>)energyInclusion, (std::vector<double>)energyExclusion,
              (std::vector<unsigned int>)uncoupledSubspace, (bool)fullFDEc, (Options::LRSCF_TYPE)loadType,
              (bool)localMO, (Options::GAUGE)gauge, (std::vector<double>)gaugeOrigin, (std::vector<double>)frequencies,
              (std::vector<double>)frequencyRange, (double)damping, (std::vector<unsigned int>)couplingPattern,
              (double)localVirtualOrbitals, (double)envVirtualOrbitals, (bool)saveResponseMatrix)
 public:
  EmbeddingSettings embedding;
};

/**
 * @class  LRSCFTask LRSCF.h
 * @brief  Perform a LRSCF calculation
 */
template<Options::SCF_MODES SCFMode>
class LRSCFTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem The main (active) system.
   * @param environmentSystems Environment systems.
   */
  LRSCFTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
            const std::vector<std::shared_ptr<SystemController>>& passiveSystems = {});
  /**
   * @brief Default destructor.
   */
  virtual ~LRSCFTask() = default;
  /**
   * @see Task.h
   */
  void run();
  /**
   * @brief Parse the settings to the task settings.
   * @param c The task settings.
   * @param v The visitor which contains the settings strings.
   * @param blockname A potential block name.
   */
  void visit(LRSCFTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
    }
    else if (!c.embedding.visitSettings(v, blockname)) {
      throw SerenityError((string) "Unknown block in LRSCFTaskSettings: " + blockname);
    }
  }
  /**
   * @brief The settings/keywords for the LRSCFTask: \n
   *        -nEigen: Number of eigenstates to be determined in a supermolecular or uncoupled (FDEu) subsystem response
   *                 calculation. Note that the number of excitations in a coupled subsystem response calculation
   *                 (FDEc) is determined by the number and choice of uncoupled transitions which were determined
   * before. In order to supress the calculation of eigenstates in FDEu/FDEc calculations, manually set to 0.
   *        -pseudoDensityThreshold: A prescreening threshold. Often, the matrix of guess vectors is rather sparse
   *                                 and has contributions from a few subsystems only, i.e. the density matrices of
   *                                 pure environment systems will be close to zero. If the maximum density matrix
   *                                 element of all density matrices obtained from the sets of guess vectors is lower
   *                                 than this threshold, the calculation is skipped.
   *        -convThresh: Convergence criterion for iterative eigenvalue solver. The iterative solution of the response
   *                     problem is stopped when the residual norm of all desired roots falls below this threshold.
   *        -maxCycles: Maximum number of iterations for the iterative eigenvalue solver.
   *        -maxSubspaceDimension: Maximum dimension of the subspace used in iterative eigenvalue solver.
   *        -responseType: Type of the LRSCF problem. Used to determine if a hermititian or non-
   *                       Hermitian problem is solved. Can enforce RPA type on pure TDDFT
   *        -tda: If true, use the Tamm--Dancoff-Approximation
   *        -dominantThresh: Orbital transitions with squared coefficients (multiplied by a factor of 100) larger than
   *                         dominantThresh are considered dominant and their contribution is written into the output.
   *        -superSystemGrid: Uses supersystem Grid for the evaluation of the Kernel Contribution in the subsystem TDDFT
   * case -func: IF another function should be used for the Kernel than for the system during SCF -analysis: If false,
   * LRSCF analysis and excitation / CD spectra will be suppressed -setAlpha: Set alpha reference orbitals given a set
   * of reference orbitals (e.g. from a SCF calculation) -excludeAlpha: Exclude alpha orbitals from provided reference
   *        -setBeta: Set beta reference orbitals given a set of reference orbitals. Note, that this keyword will
   *                  be ignored when performing a spin-restricted LRSCF calculation. For spin-restricted calculations,
   *                  use setAlpha.
   *        -excludeBeta: Exclude alpha orbitals from provided reference. Note, that this keyword will
   *                      be ignored when performing a spin-restricted LRSCF calculation. For spin-restricted
   *                      calculations, use exludeAlpha.
   *        -besleyAtoms: Number of Besley Atoms (the first n atoms will be taken from the xyz file)
   *        -besleyCutoff: Besley Cutoff for Occ and Virt
   *        -excludeProjection: Exclude all shifted orbitals from the set of reference orbitals.
   *        -energyInclusion: Includes all orbitals within the given energy interval.
   *        -energyExclusion: Excludes all orbitals outside the given energy interval.
   *        -uncoupledSubspace: Uncoupled subspace for the FDEc-LRSCF problem. Given a set of active subsystems, a
   * subspace of excitations vectors can be defined by uncoupledSubspace { 2 1 2 3 4 8 10} where the first number n
   * gives the number of states of that subsystem and the following n numbers define the respective transitions. For
   * this example, vectors 1 and 2 are taken from subsystem one,  and vectors 4, 8 and 10 are chosen from active
   * subsystem 2. If uncoupledSubspace is not set, all uncoupled vectors will be used to span the subspace. -fullFDEc:
   * Solve full FDEc problem using approximate solutions as initial guess -loadType: Reference states used to build
   * unitary transformation matrix for FDEc calculations -localMO: Perform supersystem TDDFT with local orbitals
   *        -gaugeOrigin: The gauge origin for dipole integrals (important for properties)
   *        -gauge: The gauge for response properties, i.e. length or velocity
   *        -frequencies: Frequencies for which dynamic polarizabilities and optical rotation are to be calculated
   *        -frequencyRange: Frequency range expects 3 parameters: <start> <stop> <stepwidth>
   *        -damping: Damping parameter for response properties (finite lifetime effects), e.g. broadens the absorption
   * spectrum -couplingPattern: Sets a certain coupling pattern: For example {1 1 0 2 2 0 0 0 1} means: Read in coupled
   * solution of 1 and 2 and uncoupled solution from 1 and perform FDEc step IMPORTANT: This needs to be done in order
   * coupled -> uncoupled With ordered numbers 1 ... 2 ... 3 etc and the numbers refer to the ordering of the systems in
   * the input (1: first system etc ...) -localVirtualOrbitals: Select canonical virtual orbitals located on the active
   * subsystem based on a modified overlap criterion -envVirtualOrbitals: Select canonical virtual orbitals located on
   * the environment subsystems based on a modified overlap criterion -saveResponseMatrix: Saves TDA response matrix if
   * requested
   *
   */
  LRSCFTaskSettings settings;

  std::vector<std::shared_ptr<SystemController>> _act;
  std::vector<std::shared_ptr<SystemController>> _env;
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;
  // If a special coupling matrix is used:
  Eigen::MatrixXi _couplingPatternMatrix;
  // Special pattern for loading references if coupled and uncoupled vectors from the same system are needed
  // follows the order of the couplingpattern
  std::vector<Options::LRSCF_TYPE> _referenceLoadingType;
};

} /* namespace Serenity */

#endif /* LRSCFTASK_H_*/
