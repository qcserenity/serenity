/**
 * @file   Kernel.h
 *
 * @date   Mar 19, 2017
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

#ifndef LRSCF_KERNEL
#define LRSCF_KERNEL

/* Include Serenity Internal Headers */
#include "data/DoublySpinPolarizedData.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "dft/functionals/wrappers/PartialDerivatives.h"
#include "math/Derivatives.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "system/SystemController.h"

namespace Serenity {
/**
 * @class Kernel Kernel.h
 * @brief Calculates and stores the derivatives needed for construction of the exchange-correlation and
 *        non-additive kernels. Holds objects to calculate
 *        \f[
 *          f_{xc}^{IJ} = \frac{\delta^2 E_{xc}}{\delta \rho_{tot}^2} + \frac{\delta^2 T}{\delta \rho_{tot}^2} +
 * \delta_{IJ} \frac{\delta^2 E_{xc}}{\delta \rho_{I}^2}
 *                        - \delta_{IJ} \frac{\delta^2 E_{xc}}{\delta \rho_{I}^2} - \frac{\delta^2 T}{\delta \rho_{I}^2}
 *        \f]
 *        for subsystems I and J. Here, the terms \f$ \frac{\delta^2 E_{xc}}{\delta \rho_{tot}^2} \f$, \f$
 * \frac{\delta^2 T}{\delta \rho_{tot}^2} \f$, \f$ - \delta_{IJ} \frac{\delta^2 E_{xc}}{\delta \rho_{I}^2} \f$ and \f$
 * \frac{\delta^2 T}{\delta \rho_{I}^2} \f$ are evaluated using the specified non-additive density functionals. The \f$+
 * \delta_{IJ} \frac{\delta^2 E_{xc}}{\delta \rho_{I}^2} \f$ is evaluated using the exchange-correlation functional
 * specified in the DFT settings of the respective systems. The second functional derivative of a GGA type density
 * functional of the general type \f$ E = \int \epsilon(\rho,\nabla \rho) d\mathbf{r}\f$ evaluates to \f[ \frac{\delta^2
 * E}{\delta \rho^2} = \frac{\partial^2 \epsilon}{\rho^2} - \nabla \frac{\partial^2 \epsilon}{\partial \rho \partial
 * \nabla \rho} + \sum_{ij} \nabla_i \nabla_j \frac{\partial^2 \epsilon}{\partial \nabla_i \rho \partial \nabla_j \rho}
 * \; . \f] For the matrix elements of the kernel in the adiabatic approximantion, which are needed for the response
 * matrix, we use integration by parts to move the derivative operators which than solely act on the basis function
 * products. It is thus only necessary to store \f$ \frac{\partial^2 \epsilon}{\rho^2} \f$, \f$ \frac{\partial^2
 * \epsilon}{\partial \rho \partial \nabla \rho} \f$ and \f$ \frac{\partial^2 \epsilon}{\partial \nabla_i \rho \partial
 * \nabla_j \rho} \f$ for all \f$ i,j\f$ . For more information, including also the unrestricted expressions, see e.g.
 * M. Boeckers, J. Neugebauer, 2018.
 *
 */
template<Options::SCF_MODES T>
class Kernel {
 public:
  /**
   * @brief Constructor
   *
   * @param act A vector of SystemControllers for which the derivatives shall be calculated.
   * @param env A vector of environment SystemControllers needed for non-additive contributions.
   * @param superSystemGrid If true, use grid of the supersystem to store the derivatives (Default: use grid of active
   * system)
   * @param naddKinFunc A non-additive kinetic energy functional
   * @param naddXCFunc A non-additive exchange-correlation functional
   * @param func The functional used for the kernel evaluation.
   */
  Kernel(std::vector<std::shared_ptr<SystemController>> act, std::vector<std::shared_ptr<SystemController>> env,
         bool superSystemGrid, CompositeFunctionals::KINFUNCTIONALS naddKinFunc,
         CompositeFunctionals::XCFUNCTIONALS naddXCFunc, CompositeFunctionals::XCFUNCTIONALS func);
  virtual ~Kernel() = default;

  ///@brief Returns the PP contribution to the kernel
  const std::shared_ptr<DoublySpinPolarizedData<T, Eigen::VectorXd>> getPP(unsigned int I, unsigned int J,
                                                                           unsigned int blockSize, unsigned int iGridStart);

  ///@brief Returns the GG contribution to the kernel
  const std::shared_ptr<Hessian<DoublySpinPolarizedData<T, Eigen::VectorXd>>>
  getGG(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart);

  ///@brief Returns the PG contribution to the kernel
  const std::shared_ptr<Gradient<DoublySpinPolarizedData<T, Eigen::VectorXd>>>
  getPG(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart);

  ///@brief Returns true if GGA functional is used
  bool isGGA() {
    return _gga;
  };

  /**
   * @brief Returns a shared_ptr of the GridController
   *
   * @todo Let each subsystem have its own GridController!
   */
  std::shared_ptr<GridController> getGridController() {
    return _gridController;
  };

  ///@brief Returns the blocksize.
  unsigned int getBlocksize(unsigned int subsystemNumber);
  ///@brief Returns the Basisfunction Radial Threshold
  unsigned int getbasFuncRadialThreshold(unsigned int subsystemNumber);
  ///@brief Returns the Basisfunction Block Average Value
  double getblockAveThreshold(unsigned int subsystemNumber);

 private:
  // A vector holding all (active and environment) systemControllers
  std::vector<std::shared_ptr<SystemController>> _systems;
  // If true, use super-system grid
  bool _superSystemGrid;

  // Non-additive functionals.
  CompositeFunctionals::KINFUNCTIONALS _naddKinFunc;
  CompositeFunctionals::XCFUNCTIONALS _naddXCFunc;
  // If _func is none, the functional for the active system is taken from settings.dft.
  // Else, all active functionals are evaluated with _func.
  CompositeFunctionals::XCFUNCTIONALS _func;

  // The grid controller controlling the grid on which the kernel is calculated and stored
  std::shared_ptr<GridController> _gridController;

  // True if gga functional is used
  bool _gga;

  // Data objects for each system in _systems
  std::map<std::string, DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>> _pp;
  std::map<std::string, Hessian<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>> _gg;
  std::map<std::string, Gradient<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>> _pg;
  // Data objects needed in the embedding case; evaluated using the total density
  std::shared_ptr<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>> _pptot;
  std::shared_ptr<Hessian<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>> _ggtot;
  std::shared_ptr<Gradient<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>> _pgtot;

  // Calculates the derivatives
  void calculateDerivatives();

  // Stores the derivatives
  void storeDerivatives(FunctionalData<T>& funcData, DensityOnGrid<T>& p, CompositeFunctionals::CLASSES funcClass,
                        DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>& pp,
                        Hessian<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>& gg,
                        Gradient<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>& pg, int pm = 1);

  // Returns the density on grid controller
  std::shared_ptr<DensityMatrixDensityOnGridController<T>> getDensityOnGridController(std::shared_ptr<SystemController> sys);
};

} /* namespace Serenity */

#endif /* LRSCF_KERNEL */
