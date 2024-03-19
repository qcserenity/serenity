/**
 * @file   ConvergenceController.h
 *
 * @date   28. Dezember 2013, 18:13
 * @author Thomas Dresselhaus, M. Boeckers
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
#ifndef CONVERGENCECONTROLLER_H
#define CONVERGENCECONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory> //smart ptr.

namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class Damper;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
class DIIS;
class ADIIS;
class EnergyComponentController;
class OneElectronIntegralController;
struct Settings;
/**
 * @class ConvergenceController ConvergenceController.h
 *
 * @brief The basic class for managing the convergence of an SCF.
 *
 * On the one hand convergence is accelerated/ensured by this class and on the other hand it checks
 * whether convergence has been reached. F\*P\*S - S\*P\*F is used as the DIIS error vector.
 */
template<Options::SCF_MODES SCFMode>
class ConvergenceController {
 public:
  /**
   * @param settings                    The system settings.
   * @param densityMatrix               The density matrix, used in convergence acceleration and error measurement.
   * @param orbitalController           Access to the underlying orbitals making up the density matrix.
   * @param oneIntController            OneElectronIntegralController, allows to access the overlap matrix.
   * @param energyComponentController   Used to get the total energy.
   */
  ConvergenceController(const Settings& settings, std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrix,
                        std::shared_ptr<OrbitalController<SCFMode>> orbitalController,
                        std::shared_ptr<OneElectronIntegralController> oneIntController,
                        std::shared_ptr<EnergyComponentController> energyComponentController);
  virtual ~ConvergenceController() = default;
  /**
   * @brief Getter for levelshift information.
   * @return Returns all the information needed for the levelshift in the orbital updater (energy and number of occupied
   * orbitals).
   */
  std::pair<Eigen::VectorXd, SpinPolarizedData<SCFMode, Eigen::VectorXd>> getLevelshift();
  /**
   * Optimizes the Fock matrix to produce way better orbitals in the next SCF cycle.
   * @param F The Fock matrix to be updated.
   */
  void accelerateConvergence(FockMatrix<SCFMode>& F, DensityMatrix<SCFMode> D);
  /**
   * Returns true if convergence has been reached.
   */
  bool checkConvergence();

 private:
  const Settings& _settings;
  std::shared_ptr<DensityMatrixController<SCFMode>> _dmatContr;
  std::shared_ptr<OrbitalController<SCFMode>> _orbitalController;
  std::shared_ptr<DensityMatrix<SCFMode>> _oldP;
  const std::shared_ptr<OneElectronIntegralController> _oneIntController;
  const std::shared_ptr<EnergyComponentController> _energyComponentController;
  double _oldEnergy;
  double _oldOneElEnergy;
  double _diisConvMeasure;
  double _rmsdOfDensity;
  std::unique_ptr<Eigen::MatrixXd> _orthoS;
  std::shared_ptr<Damper<SCFMode>> _damping;
  std::shared_ptr<DIIS> _diis;
  std::shared_ptr<ADIIS> _adiis;
  unsigned int _diisZoneStart;
  timespec _time;
  std::string _mode;
  unsigned int _cycle = 0;
  bool _first = true;
  // Get the number of convergence criteria met.
  unsigned int getNConverged();
  // Print the SCF cycle info.
  void printCycleInfo();
  /**
   *
   * @param A possibly modified overlap matrix. This is needed for EDA calculations.
   * @return The error vector [F,P].
   */
  MatrixInBasis<SCFMode> calcFPSminusSPF(FockMatrix<SCFMode>& F);

  double calcRMSDofDensity();

  bool _levelShiftInLastIteration = false;

  const unsigned int _nNeccessaryToConverge = 2;
};

} /* namespace Serenity */
#endif /* CONVERGENCECONTROLLER_H */
