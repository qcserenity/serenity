/**
 * @file   GWTask.h
 *
 * @date   24.03.2020
 * @author J. Toelle
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n
 */
#ifndef GWTASK_H_
#define GWTASK_H_
/* Include Serenity Internal Headers */
#include "geometry/Geometry.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/BasisOptions.h"
#include "settings/EmbeddingSettings.h"
#include "settings/MBPTOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"

namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace Serenity::Reflection;
struct GWTaskSettings {
  GWTaskSettings()
    : mbpttype(Options::MBPT::GW),
      gwtype(Options::GWALGORITHM::CD),
      linearized(false),
      qpiterations(0),
      eta(0.001),
      nVirt(10),
      nOcc(10),
      integrationPoints(128),
      padePoints(16),
      fermiShift(0.0),
      derivativeShift(0.002),
      imagShift(0.001),
      gridCutOff(-1.0),
      evGW(false),
      evGWcycles(5),
      diis(true),
      diisMaxStore(10),
      ConvergenceThreshold(1e-6),
      nafThresh(0),
      subsystemAuxillaryBasisOnly(false),
      freq({}),
      damping(0.2),
      gap(false),
      environmentScreening(true),
      ltconv(0),
      frozenCore(false),
      coreOnly(false),
      densFitCache(Options::DENS_FITS::RI) {
    embedding.naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::NONE;
    embedding.naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
    embedding.embeddingMode = Options::KIN_EMBEDDING_MODES::NONE;
  };
  REFLECTABLE((Options::MBPT)mbpttype, (Options::GWALGORITHM)gwtype, (bool)linearized, (unsigned int)qpiterations,
              (double)eta, (unsigned int)nVirt, (unsigned int)nOcc, (unsigned int)integrationPoints,
              (unsigned int)padePoints, (double)fermiShift, (double)derivativeShift, (double)imagShift, (double)gridCutOff,
              (bool)evGW, (unsigned int)evGWcycles, (bool)diis, (unsigned int)diisMaxStore, (double)ConvergenceThreshold,
              (double)nafThresh, (bool)subsystemAuxillaryBasisOnly, (std::vector<unsigned int>)hybrid,
              (std::vector<double>)freq, (double)damping, (bool)gap, (bool)environmentScreening, (double)ltconv,
              (bool)frozenCore, (bool)coreOnly, (Options::DENS_FITS)densFitCache)
 public:
  EmbeddingSettings embedding;
};

/**
 * @class  GWTask GWTask.h
 * @brief  Perform a GW calculation
 */
template<Options::SCF_MODES SCFMode>
class GWTask : public Task {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem The main (active) system.
   * @param environmentSystems Environment systems.
   */
  GWTask(const std::vector<std::shared_ptr<SystemController>>& activeSystems,
         const std::vector<std::shared_ptr<SystemController>>& passiveSystems = {});
  /**
   * @brief Default destructor.
   */
  virtual ~GWTask() = default;
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
  void visit(GWTaskSettings& c, set_visitor v, std::string blockname) {
    if (!blockname.compare("")) {
      visit_each(c, v);
      return;
    }
    if (c.embedding.visitAsBlockSettings(v, blockname))
      return;
    // If reached, the blockname is unknown.
    throw SerenityError((std::string) "Unknown block in GWTaskSettings: " + blockname);
  }
  /**
   * @brief The settings/keywords for GWTask:\n
   *  -mbpttype:           The type of MBPT calculation (GW/dRPA)
   *  -gwtype:             The GW algorithm (Analytic/Contour deformation)
   *  -linearized:         Whether the linearized solution is calculated
   *  -qpiterations:       The number of QP iterations for the G0W0 calculation
   *  -eta:                Imaginary shift parameter
   *  -nVirt:              The number virtual orbitals included in GW calculation
   *  -nOcc:               The number occupied orbitals included in GW calculation
   *  -integrationPoints:  The number of integration points for CD-GW and dRPA
   *  -padePoints:         The number of points used in the pade approximation for analytic continuation
   *  -fermiShift:         The initial fermi shift of the HOMO/LUMO for analytic continuation
   *  -derivativeShift:    The shift in for the evaluation of the numerical derviation of the self-energy (for
   * linearization)
   *  -imagShift:          Additional imaginary shift for numerical derivative if the derivative is larger than one
   *  -gridCutOff:         Gridcutoff to extend the grid of subsystem with grid points of environment subsystems
   *  -evGW:               Whether an evGW calculation is performed
   *  -evGWcycles:         evGW cycles
   *  -diis:               Whether a DIIS is used for convergence acceleration in qpiterations or evgw cycles
   *  -diisMaxStore:       The numbers of DIIS vectors to be stored
   *  -ConvergenceThreshold: The HOMO-LUMO gap convergence threshold for qpiterations/evGW cycles
   *  -naf:                Whether natural-auxiliary functions should be used
   *  -nafThresh:          The threshold for the naf functions
   *  -subsystemAuxillaryBasisOnly: Whether subsystem screening contributions should be calculated with subsystem
   * auxiliary basis only
   *  -freq:               Start, end, stepsize for real axes frquency of the self-energy (only working for GW-Analytic)
   *  -damping:            Damping factor for convergence acceleration -gap: Whether to shift
   * occupied and virtual orbitals not included in the GW caclulation by the gap of change of the highest/lowest
   * included occupied/virtual orbital
   *  -environmentScreening: Whether environmental screening is included in an embedded GW/dRPA calculation
   *  -ltconv:             Convergence criterion for num. int with Laplace transformation.
   *  -frozenCore:         Whether frozenCore approximation is used or not
   *  -coreOnly:           Whether only core orbitals are taken into account
   */
  GWTaskSettings settings;
  /// @brief The active systems
  std::vector<std::shared_ptr<SystemController>> _act;
  /// @brief The environment/embedding system
  std::vector<std::shared_ptr<SystemController>> _env;
  /// @brief GW orbital correlation energies
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _correlationEnergies;
  /// @brief RPA orbital correlation energy
  double _correlationEnergy;

 private:
  /// @brief prints out all needed information of a GW calculation
  void printGWResult(SpinPolarizedData<SCFMode, Eigen::VectorXd>& orbEigenValues,
                     SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEnergy,
                     SpinPolarizedData<SCFMode, Eigen::VectorXd>& vxc_energies,
                     SpinPolarizedData<SCFMode, Eigen::VectorXd>& x_energies,
                     SpinPolarizedData<SCFMode, Eigen::VectorXd>& correlation, SpinPolarizedData<SCFMode, Eigen::VectorXd>& z,
                     SpinPolarizedData<SCFMode, Eigen::VectorXd>& dsigma_de, int& start, int& end);
  /// @brief determination of the structure for the RI integrals in the GW calculation
  std::shared_ptr<Geometry> superMolGeo();
};

} /* namespace Serenity */

#endif /* GWTASK_H_*/