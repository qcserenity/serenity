/**
 * @file   ConvergenceController.h
 *
 * @date   28. Dezember 2013, 18:13
 * @author Thomas Dresselhaus, M. Boeckers
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
#ifndef CONVERGENCECONTROLLER_H
#define	CONVERGENCECONTROLLER_H
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "settings/Options.h"
#include "settings/Settings.h"
#include "data/SpinPolarizedData.h"
#include "data/OrbitalController.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declarations */
template<Options::SCF_MODES T>class Damper;
template<Options::SCF_MODES T>class DensityMatrixController;
class DIIS;
class ADIIS;
class EnergyComponentController;
class OneElectronIntegralController;
/**
 * @class ConvergenceController ConvergenceController.h
 *
 * @brief The basic class for managing the convergence of an SCF.
 *
 * On the one hand convergence is accelerated/ensured by this class and on the other hand it checks
 * whether convergence has been reached. F\*P\*S - S\*P\*F is used as the DIIS error vector.
 * 
 */
template<Options::SCF_MODES T>
class ConvergenceController {
public:
  /**
   * @param fockMatrix          The Fock matrix (potential matrix) used in convergence acceleration
   *                            and error measurement
   * @param densityMatrix       Also used for error measurement.
   * @param overlapIntegrals    Also used for error measurement.
   * @param forceExactFockBuild a magical switch that gives the possibility that SOMEWHERE, at SOME
   *                            POINT in the code an exact fock build can be enforced (instead of
   *                            incremental)
   * @param electronicEnergy    Needed for convergence check.
   * @param deltaEConvThreshold Convergence threshold for the (absolute) energy difference.
   */
  ConvergenceController(
      const Settings& settings,
      std::shared_ptr< DensityMatrixController<T> > densityMatrix,
      std::shared_ptr<OrbitalController<T> > orbitalController,
      std::shared_ptr<OneElectronIntegralController> oneIntController,
      std::shared_ptr<EnergyComponentController> energyComponentController);
  virtual ~ConvergenceController() = default;
  /**
   * @brief Getter for levelshift information.
   * @return Returns all the information needed for the levelshift in the orbital updater (energy and number of occupied orbitals).
   */
  std::pair<Eigen::VectorXd,SpinPolarizedData<T, Eigen::VectorXd > > getLevelshift();
  /**
   * Optimizes the fock matrix to produce way better orbitals in the next SCF cycle.
   * @param F The Fock matrix to be updated.
   */
  void accelerateConvergence(FockMatrix<T>& F);

  void setOrthoS(Eigen::MatrixXd& orthoS) {
    _orthoS.reset(new Eigen::MatrixXd (orthoS));
  }
  /**
   * Returns true if convergence has been reached.
   */
  bool checkConvergence();
  void reinitDIIS();
private:
  const Settings& _settings;
  std::shared_ptr< DensityMatrixController<T> > _dmatContr;
  std::shared_ptr<OrbitalController<T> > _orbitalController;
  std::shared_ptr<DensityMatrix<T> > _oldP;
  const std::shared_ptr<OneElectronIntegralController> _oneIntController;
  const std::shared_ptr<EnergyComponentController> _energyComponentController;
  double _oldEnergy;
  double _oldOneElEnergy;
  double _diisConvMeasure;
  double _rmsdOfDensity;
  std::unique_ptr<Eigen::MatrixXd> _orthoS;
  std::shared_ptr<Damper<T> > _damping;
  std::shared_ptr<DIIS> _diis;
  std::shared_ptr<ADIIS> _adiis;
  unsigned int _diisZoneStart;
  timespec _time;
  std::string _mode;
  unsigned int _cycle = 0;
  bool _first = true;
  /**
   *
   * @param A possibly modified overlap matrix. This is needed for EDA calculations.
   * @return The error vector [F,P].
   */
  SpinPolarizedData<T, Eigen::MatrixXd > calcFPSminusSPF(FockMatrix<T>& F);

  double calcRMSDofDensity();
};

} /* namespace Serenity */
#endif	/* CONVERGENCECONTROLLER_H */
