/**
 * @file CCSD.h
 *
 * @date Mar 31, 2016
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

#ifndef ELECTRONICSTRUCTURECALCULATIONS_POSTHF_CC_CCSD_H_
#define ELECTRONICSTRUCTURECALCULATIONS_POSTHF_CC_CCSD_H_

/* Include Serenity Internal Headers */
#include "math/Matrix.h"
#include "math/RegularRankFourTensor.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;

/**
 * @class CCSD CCSD.h
 * @brief A class that calculates the CCSD energy.
 */
class CCSD {
 public:
  /**
   * @brief Constructor.
   * @param systemController    The system of interest.
   */
  CCSD(std::shared_ptr<SystemController> systemController, double normThreshold, unsigned int maxCycles);
  /**
   * @brief Default Destructor.
   */
  virtual ~CCSD() = default;

  /**
   * @brief Converges the CCSD amplitudes and returns the energies.
   * @return A pair containing the MP2 and the CCSD energy corrections (in that order).
   */
  std::pair<double, double> calculateElectronicEnergyCorrections();

 private:
  /**
   * @brief Handles AO->MO transformation and storage.
   * TODO: this could be more efficient also not all of the MO integrals
   *       need to be stored.
   */
  void prepERIS();

  /**
   * @brief Initializes singles and double amplitudes.
   *        Single amplitudes are set to zero and double amplitudes are initialized
   *        with the MP2 amplitudes.
   */
  void initializeAmplitudes();

  /**
   * @brief Updates F_oo F_ov and F_vv.
   */
  void updateF();

  /**
   * @brief  Runs one update of the singles and the doubles amplitudes.
   * @return The norm of the change in the amplitudes.
   */
  double updateAmplitudes();

  /**
   * @brief Calculates the energy correction associated with the current amplitudes.
   * @return The CCSD energy correction.
   */
  double calculateCCSDEnergyCorrection();

 protected:
  std::unique_ptr<RegularRankFourTensor<double>> _eris;
  /// @brief Singles Amplitudes
  std::unique_ptr<Matrix<double>> _t1;
  /// @brief Doubles Amplitudes
  std::unique_ptr<Matrix<Matrix<double>>> _t2;

  std::shared_ptr<SystemController> _systemController;
  bool _converged;

 private:
  double _MP2EnergyCorrection;

  std::unique_ptr<Matrix<double>> _foo;
  std::unique_ptr<Matrix<double>> _fov;
  std::unique_ptr<Matrix<double>> _fvv;
  double _normThreshold;
  unsigned int _maxCycles;
  // DIIS _diis(6,);
};

} /* namespace Serenity */

#endif /* ELECTRONICSTRUCTURECALCULATIONS_POSTHF_CC_CCSD_H_ */
