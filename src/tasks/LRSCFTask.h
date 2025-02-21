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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the GNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef LRSCFTASK_H_
#define LRSCFTASK_H_
/* Include Serenity Internal Headers */
#include "settings/BasisOptions.h"
#include "settings/EmbeddingSettings.h"
#include "settings/LRSCFOptions.h"
#include "settings/Settings.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>

namespace Serenity {

/* Forward declarations */
class SystemController;

template<Options::SCF_MODES SCFMode>
class LRSCFController;

struct LRSCFTaskSettings {
  LRSCFTaskSettings()
    : nEigen(4),
      conv(1.0e-5),
      maxCycles(100),
      maxSubspaceDimension(1e6),
      dominantThresh(0.85),
      func(CompositeFunctionals::XCFUNCTIONALS::NONE),
      analysis(true),
      besleyAtoms(0),
      besleyCutoff({}),
      excludeProjection(false),
      uncoupledSubspace({}),
      fullFDEc(false),
      loadType(Options::LRSCF_TYPE::UNCOUPLED),
      gauge(Options::GAUGE::LENGTH),
      gaugeOrigin({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()}),
      frequencies({}),
      frequencyRange({}),
      damping(0.0),
      couplingPattern({}),
      method(Options::LR_METHOD::TDDFT),
      diis(false),
      diisStore(50),
      preopt(1e-3),
      cctrdens(false),
      ccexdens(false),
      sss(1.0),
      oss(1.0),
      nafThresh(0),
      samedensity({}),
      subsystemgrid({}),
      rpaScreening(false),
      restart(false),
      densFitJ(Options::DENS_FITS::RI),
      densFitK(Options::DENS_FITS::NONE),
      densFitLRK(Options::DENS_FITS::NONE),
      densFitCache(Options::DENS_FITS::RI),
      transitionCharges(false),
      partialResponseConstruction(false),
      grimme(false),
      adaptivePrescreening(true),
      frozenCore(false),
      frozenVirtual(0),
      coreOnly(false),
      ltconv(0),
      aocache(true),
      triplet(false),
      noCoupling(false),
      approxCoulomb({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}),
      scfstab(Options::STABILITY_ANALYSIS::NONE),
      stabroot(0),
      stabscal(0.5),
      noKernel(false),
      excGradList({}),
      hypthresh(1e-12) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::NONE;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
    // Don't use two grids by default.
    grid.smallGridAccuracy = 4;
    grid.accuracy = 4;
  };
  REFLECTABLE((unsigned int)nEigen, (double)conv, (unsigned int)maxCycles, (unsigned int)maxSubspaceDimension,
              (double)dominantThresh, (CompositeFunctionals::XCFUNCTIONALS)func, (bool)analysis, (unsigned int)besleyAtoms,
              (std::vector<double>)besleyCutoff, (bool)excludeProjection, (std::vector<unsigned int>)uncoupledSubspace,
              (bool)fullFDEc, (Options::LRSCF_TYPE)loadType, (Options::GAUGE)gauge, (std::vector<double>)gaugeOrigin,
              (std::vector<double>)frequencies, (std::vector<double>)frequencyRange, (double)damping,
              (std::vector<unsigned int>)couplingPattern, (Options::LR_METHOD)method, (bool)diis, (unsigned)diisStore,
              (double)preopt, (bool)cctrdens, (bool)ccexdens, (double)sss, (double)oss, (double)nafThresh,
              (std::vector<unsigned int>)samedensity, (std::vector<unsigned int>)subsystemgrid, (bool)rpaScreening,
              (bool)restart, (Options::DENS_FITS)densFitJ, (Options::DENS_FITS)densFitK, (Options::DENS_FITS)densFitLRK,
              (Options::DENS_FITS)densFitCache, (bool)transitionCharges, (bool)partialResponseConstruction,
              (bool)grimme, (bool)adaptivePrescreening, (bool)frozenCore, (double)frozenVirtual, (bool)coreOnly,
              (double)ltconv, (bool)aocache, (bool)triplet, (bool)noCoupling, (std::vector<double>)approxCoulomb,
              (Options::STABILITY_ANALYSIS)scfstab, (unsigned)stabroot, (double)stabscal, (bool)noKernel,
              (std::vector<unsigned int>)excGradList, (double)hypthresh)

  EmbeddingSettings embedding;
  GRID grid;
  CUSTOMFUNCTIONAL customFunc;

  /**
   * @brief Print the LRSCFTaskSettings to file.
   * @param filename The filename to print to.
   */
  void printSettings(std::basic_string<char> filename);
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
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    if (c.grid.visitAsBlockSettings(v, blockname)) {
      return;
    }
    if (c.customFunc.visitAsBlockSettings(v, blockname)) {
      return;
    }
    throw SerenityError((std::string) "Unknown block in LRSCFTaskSettings: " + blockname);
  }

  /**
   * @brief The settings/keywords for the LRSCFTask: \n
   *        -nEigen: Number of eigenstates to be determined in a supermolecular or uncoupled (FDEu) subsystem response
   *                 calculation. Note that the number of excitations in a coupled subsystem response calculation
   *                 (FDEc) is determined by the number and choice of uncoupled transitions which were determined
   * before. In order to supress the calculation of eigenstates in FDEu/FDEc calculations, manually set to 0.
   *        -conv: Convergence criterion for iterative eigenvalue solver. The iterative solution of the response
   *                     problem is stopped when the residual norm of all desired roots falls below this threshold.
   *        -maxCycles: Maximum number of iterations for the iterative eigenvalue solver.
   *        -maxSubspaceDimension: Maximum dimension of the subspace used in iterative eigenvalue solver.
   *        -dominantThresh: Orbital transitions with squared coefficients larger than
   *                         dominantThresh are considered dominant and their contribution is written into the output.
   * case -func: IF another function should be used for the Kernel than for the system during SCF -analysis: If false,
   * LRSCF analysis and excitation / CD spectra will be suppressed
   *        -besleyAtoms: Number of Besley Atoms (the first n atoms will be taken from the xyz file)
   *        -besleyCutoff: Besley Cutoff for Occ and Virt
   *        -excludeProjection: Exclude all shifted orbitals from the set of reference orbitals.
   *        -uncoupledSubspace: Uncoupled subspace for the FDEc-LRSCF problem. Given a set of active subsystems, a
   * subspace of excitations vectors can be defined by uncoupledSubspace { 2 1 2 3 4 8 10} where the first number n
   * gives the number of states of that subsystem and the following n numbers define the respective transitions. For
   * this example, vectors 1 and 2 are taken from subsystem one,  and vectors 4, 8 and 10 are chosen from active
   * subsystem 2. If uncoupledSubspace is not set, all uncoupled vectors will be used to span the subspace. -fullFDEc:
   * Solve full FDEc problem using approximate solutions as initial guess -loadType: Reference states used to build
   * unitary transformation matrix for FDEc calculations
   *        -gaugeOrigin: The gauge origin for dipole integrals (important for properties)
   *        -gauge: The gauge for response properties, i.e. length or velocity
   *        -frequencies: Frequencies for which dynamic polarizabilities and optical rotation are to be calculated
   *        -frequencyRange: Frequency range expects 3 parameters: <start> <stop> <stepwidth>
   *        -damping: Damping parameter for response properties (finite lifetime effects), e.g. broadens the absorption
   * spectrum -couplingPattern: Sets a certain coupling pattern: For example {1 1 0 2 2 0 0 0 1} means: Read in coupled
   * solution of 1 and 2 and uncoupled solution from 1 and perform FDEc step IMPORTANT: This needs to be done in order
   * coupled -> uncoupled With ordered numbers 1 ... 2 ... 3 etc and the numbers refer to the ordering of the systems in
   * the input (1: first system etc ...)
   *         -method: Determines the method to be used. The default is tddft. Also available are: tda, cc2, cisdinf,
cisd, adc2.
   *         -diis: Specifies whether the nonlinear eigenvalue solver uses a DIIS after preoptimization. If false, the
quasi-linear Davidson algorithm will be used until conv is reached.
   *         -diisStore: Specifies how many diis vectors can be stored (for CC2 ground state and nonlinear eigenvalue
solver). Default is 50.
   *         -preopt: Convergence threshold for the preoptimization of eigenvectors in nonlinear eigenvalue solvers for
CC2/ADC(2). Up to this threshold, a quasi-linear Davidson   algorithm will be used, after this a DIIS eigenvalue solver
is turned on and converged to the parameter given by conv. The default 1e-3.
   *         -cctrdens: Calculate transition moments for CC2.
   *         -ccexdens: Calculate excited-state densities and properties for CC2.
   *         -sss: Scaling parameter for same-spin contributions (CC2/ADC(2)). The default is 1.0.
   *         -oss: Scaling parameter for opposite-spin contributions (CC2/ADC(2)). The default is 1.0.
   *         - nafThresh: Truncates the three-center MO integral basis using the natural auxiliary function technique.
The default is false . Treshold for truncation. The smaller, the fewer NAFs are truncated. The default is 1e-2 .
   *         -samedensity: If two subsystems are used in the calculation with the same occupied but different virtual
orbital spaces, the keyword need to be set due to the fact that the kernel needs only to be evaluated with one of the
densities. Expects a list of arguments {Number of subsystems used for kernel evaluation}
   *         -subsystemgrid: Only includes the grid points associated with atoms of the specified subsystems in the
kernel evaluation (subsystems are numbered in the input order starting from 1).
   *         -rpaScreening: Performs the exchange integral evaluation with static RPA screening.
   *                          The default is false . Note: If environmental
subsystems are specified their screening contribution is included approximately.
   *         -restart: Tries restarting from (preferably converged) eigenpairs that the tasks looks for in the system
folder.
   *         -densFitJ: Density fitting for Coulomb sigma vectors.
   *         -densFitK: Density fitting for exchange sigma vectors.
   *         -densFitLRK: Density fitting for long-range exchange sigma vectors.
   *         -densFitCache: Density fitting for RIIntegrals, e.g. for CC2/ADC(2).
   *         -transitionCharges: Calculates transition charges and stores them on disk.
   *         -partialResponseConstruction: Exploit symmetry of FDEc.
   *         -grimme: Invokes simplified TDA/TDDFT for this task.
   *         -approxCoulomb: This keyword accepts up to two doubles which are distance thresholds. Below the first
   *          parameter, Coulomb interactions between two subsystems will be calculated either NORI or RI (using
   *          the union of the aux. bases of the two subsystems). Between the two parameters, not the union but
   *          the aux. bases alone will be used, and above the last threshold, simplified TDA will be used to
   *          calculate the Coulomb integrals. If only one threshold is given, the second is set to infinity.
   *         -triplet: Requests triplet instead of singlet excitations for a restricted reference.
   *         -scfstab: Switch for different types of stability analyses (real, nonreal, spinflip).
   *         -stabroot: Instruct to rotate the orbitals along the direction devised by the instability indexed by
stabroot (usually 1, i.e. the lowest). Another SCF needs to be run.
   *         -stabscal: Mixing parameter for the new orbitals: C_new = C_old * U with U = exp(stabscal * F).
   *         -noKernel: Can be used to explicitly turn off the numerical integration of the XC and kinetic kernels.
   *                    This keyword combined with "auxclabel MINPARKER_S/SP" in the basis block of the system
   *                    settings can be used to do TDDFT-ris calculations (10.1021/acs.jpclett.2c03698)
   *         -excGradList: A 1-based list of excited states for which gradients are to be calculated. Excited-State
gradients are calculated when this list is not empty.
   *         -hypthresh: Density threshold for TDDFT gradient calculations. If the density at a point falls below this
threshold, the third derivatives of the xc functional (sometimes called the hyperkernel) are set to zero at this point.
   */
  LRSCFTaskSettings settings;

  /**
   * @brief Returns this task's LRSCF controllers.
   * @return See above.
   */
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> getLRSCFControllers();

  /**
   * @brief Returns results of this task (excitation energies and transition moments).
   * @return A (nEigen x 6) matrix. The first column contains the excitation energies in atomic units, the second and
   * third contain oscillator strengths and the fourth, fifth and sixth column contain rotatory strengths (see
   * ExcitationSpectrum).
   */
  const Eigen::MatrixXd& getTransitions();

  /**
   * @brief Returns results of this task (response properties).
   * @return See above.
   */
  const std::vector<std::tuple<double, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d>>& getProperties();

 private:
  // Active systems.
  std::vector<std::shared_ptr<SystemController>> _act;

  // Environment systems.
  std::vector<std::shared_ptr<SystemController>> _env;

  // LRSCFController for each system.
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;

  // If a special coupling matrix is used:
  Eigen::MatrixXi _couplingPattern;

  // Special pattern for loading references if coupled and uncoupled vectors from the same system are needed
  // follows the order of the couplingpattern
  std::vector<Options::LRSCF_TYPE> _referenceLoadingType;

  // A matrix to lazily store the results of this task (excitation energies and osc./rot. strengths).
  Eigen::MatrixXd _excitations;

  // A matrix to lazily store the results of this task (response properties).
  std::vector<std::tuple<double, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d>> _properties;
};

} /* namespace Serenity */

#endif /* LRSCFTASK_H_*/
