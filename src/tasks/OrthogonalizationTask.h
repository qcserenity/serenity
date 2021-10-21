/**
 * @file OrthogonalizationTask.h
 *
 * @date Aug 06, 2020
 * @author Anja Massolle
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

#ifndef TASKS_ORTHOGONALIZATIONTASK_H_
#define TASKS_ORTHOGONALIZATIONTASK_H_

/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "settings/Options.h"
#include "settings/OrthogonalizationOptions.h"
#include "settings/Reflection.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
using namespace Serenity::Reflection;

struct OrthogonalizationTaskSettings {
  OrthogonalizationTaskSettings()
    : orthogonalizationScheme(Options::ORTHOGONALIZATION_ALGORITHMS::LOEWDIN), maxIterations(50) {
  }

  REFLECTABLE((Options::ORTHOGONALIZATION_ALGORITHMS)orthogonalizationScheme, (unsigned int)maxIterations)
};
/**
 * @class  OrthogonalizationTask OrthogonalizationTask.h
 * @brief  Orthogonalize Molecular Orbitals. Possible orthogonalization schemes are:\n
 *  - Löwdins symmetric orthogonalization, according to: Löwdin, Per‐Olov, J. Chem. Phys. 18, 1950, 365-375 and \n
 *                                                      Löwdin, Per-Olov, Advances in quantum chemistry. Vol. 5.
 * Academic Press, 1970. 185-199.
 *  - Pipeks iterative orthogonalization scheme, according to: Pipek, J. Int. J. Quantum Chem. 1985,
 * 27, 527–546.
 *  - Broers orthogonalization scheme based on corresponding orbitals, according to: Broer, R. Int. J.
 * Quantum Chem. 1993, 45, 587-590.
 */
template<Options::SCF_MODES SCFMode>
class OrthogonalizationTask : public Task {
 public:
  /**
   * @brief Task which orthogonalize orbitals for a supersystem.
   * @param superSystem Supersystem which also holds the orthogonalized orbitals after the transformation.
   * @param systemController subsystems used to create supersystem orbitals.
   */
  OrthogonalizationTask(std::vector<std::shared_ptr<SystemController>> systemController,
                        std::shared_ptr<SystemController> superSystem = nullptr);
  /**
   * @brief Default destructor.
   */
  virtual ~OrthogonalizationTask() = default;
  /**
   * @see Task
   */
  void run();
  /**
   * @brief Checks if the orbitals are orthogonal.
   * @return the result as a boolian
   */
  bool isOrtho() {
    return _ortho;
  }

  /**
   * @brief The settings/keywords for the Orthogonalization Task.\n
   * -orthogonalizationScheme the applied orthogonalization procedure
   * -maxIterations Maximum number of iterations performed in the Pipek scheme
   */
  OrthogonalizationTaskSettings settings;

 private:
  // the subsystems
  std::vector<std::shared_ptr<SystemController>> _systemController;
  // the supersystem
  std::shared_ptr<SystemController> _superSystem;
  // function which normalizes the orbitals
  CoefficientMatrix<SCFMode> normalize(CoefficientMatrix<SCFMode> coeff);
  // calculates the matrix T used for the pipek algorithm
  Eigen::MatrixXd calcT(Eigen::MatrixXd coeff, unsigned int nOcc);
  // calculates epsilon used for the pipek algorithm
  double calcEpsilon(unsigned int nOcc, Eigen::MatrixXd coeff);
  // calcualtes gamma used for the pipek algorithm
  double calcGamma(unsigned int l, unsigned int n, double epsilon);
  // checks if the orbitals are orthogonal
  bool checkOrtho(CoefficientMatrix<SCFMode> coeff, SpinPolarizedData<SCFMode, unsigned int> nOcc, bool print = false);
  // true if the orbitals are orthogonal
  bool _ortho = false;
  // AO overlap matrix
  Eigen::MatrixXd _sAO;
};

} /* namespace Serenity */

#endif /* TASKS_ORTHOGONALIZATIONTASKTASK_H_ */