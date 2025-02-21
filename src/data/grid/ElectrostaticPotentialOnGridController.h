/**
 * @file ElectrostaticPotentialOnGridController.h
 *
 * @author Moritz Bensberg
 * @date May 19, 2020
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

#ifndef DATA_GRID_ELECTROSTATICPOTENTIALONGRIDCONTROLLER_H_
#define DATA_GRID_ELECTROSTATICPOTENTIALONGRIDCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "data/grid/CoulombIntegralsOnGridController.h"
#include "data/grid/GridPotential.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "notification/ObjectSensitiveClass.h" //Notify.
/* Include Std and External Headers */
#include <memory>
#include <string>

namespace H5 {
class H5File;
} // namespace H5

namespace Serenity {
namespace HDF5 {
using H5File = H5::H5File;
} // namespace HDF5

class GridController;
template<Options::SCF_MODES SCFMode>
class DensityMatrixController;
class Geometry;
class BasisController;
class AtomCenteredBasisController;

/**
 * @class
 * @brief A controller for the electrostatic potential on a given grid. This controller uses caching of two center
 *        integrals on disk. The electrostatic potential may depend on a geometry and a density matrix controller.
 */
template<Options::SCF_MODES SCFMode>
class ElectrostaticPotentialOnGridController : public ObjectSensitiveClass<DensityMatrix<SCFMode>>,
                                               public ObjectSensitiveClass<Grid> {
 public:
  /**
   * @brief Constructor.
   * @param gridController               The grid controller.
   * @param densityMatrixController      The density matrix controller.
   * @param geometry                     The geometry.
   * @param fBaseName                    The file base name.
   * @param cacheSize                    The maximum number of two electron center sets to be stored in memory.
   * @param normalVectors                The normal vectors on the molecular surface (only needed for the electric
   * field).
   * @param atomCenteredBasisController  The atom centered basis controller of the density matrix. Only needed for
   * gradients.
   */
  ElectrostaticPotentialOnGridController(std::shared_ptr<GridController> gridController,
                                         std::shared_ptr<DensityMatrixController<SCFMode>> densityMatrixController,
                                         std::shared_ptr<Geometry> geometry, std::string fBaseName,
                                         unsigned int cacheSize = 60, std::shared_ptr<Eigen::Matrix3Xd> normalVectors = nullptr,
                                         std::shared_ptr<AtomCenteredBasisController> atomCenteredBasisController = nullptr);
  /**
   * @brief Destructor. Deletes chache-files on disk.
   */
  ~ElectrostaticPotentialOnGridController();
  /**
   * @brief Getter for the electrostatic potential.
   * @return The electrostatic potential.
   */
  const GridPotential<RESTRICTED>& getPotential();
  /**
   * @brief Getter for the density matrix controller.
   * @return The density matrix controller.
   */
  std::shared_ptr<DensityMatrixController<SCFMode>> getDensityMatrixController();
  /**
   * @brief Setter for the disk mode. If true, the electrostatic potential will be written to disk.
   * @param newMode The new disk mode.
   */
  void setDiskMode(bool newMode);
  /**
   * @brief Force reinit of the electrostatic potential.
   */
  void notify() override final {
    _electrostaticPotential = nullptr;
    _diskUpToDate = false;
  };
  /**
   * @brief Calculate the total Fock matrix rsulting from the interaction with the
   *        given set of charges.
   *
   *        The Fock matrix is expressed in the basis of the
   *        density matrix controller. This makes use of the already cached integrals
   *        which are needed for the potential evaluation at the same points as the
   *        charges are located at:\n
   *          \f$ \pmb{F} = \sum_i -p_i  \pmb{v}_i \f$\n
   *        Here the integrals for point \f$\pmb{r}_i  \f$ are denoted by\n
   *        \f$ (\pmb{v}_i)_{\mu\nu} = \langle \mu | -1/|\pmb{r}-\pmb{r}_i| | \nu \rangle\f$\n
   *        and the charges by \f$ p_i \f$ .
   * @param charges The charges (\f$ p_i \f$).
   * @return The Fock matrix.
   */
  FockMatrix<RESTRICTED> integrateFockMatrix(const GridPotential<RESTRICTED>& charges);
  /**
   * @brief Removes any integrals or potential files from disk and resets the cache.
   */
  void cleanUpDisk();
  /**
   * @brief Getter for the associated basis controller.
   */
  std::shared_ptr<BasisController> getBasisController();
  /**
   * @brief Getter for the associated grid controller.
   */
  std::shared_ptr<GridController> getGridController();
  /**
   * @brief calculate the contribution for the gradient.
   */
  Eigen::MatrixXd calculateGradient(const std::vector<unsigned int>& atomIndicesOfPoints,
                                    const GridPotential<RESTRICTED>& charges);

 private:
  // The controller for the two center integrals.
  std::shared_ptr<CoulombIntegralsOnGridController> _coulombIntegralsOnDiskController;
  // The density matrix controller.
  std::shared_ptr<DensityMatrixController<SCFMode>> _densityMatrixController;
  // The geometry.
  std::shared_ptr<Geometry> _geometry;
  // The current electrostatic potential.
  std::shared_ptr<GridPotential<RESTRICTED>> _electrostaticPotential;
  // The disk mode.
  bool _diskMode = false;
  // If true, the potential saved on disk is up to date. If false, it needs to be recalculated.
  bool _diskUpToDate = false;
  // The file base name.
  std::string _fBaseName;
  // The atom centered basis controller of the density matrix. Only needed for gradients.
  std::shared_ptr<AtomCenteredBasisController> _atomCenteredBasisController;
  // Update the electrostatic potential.
  void update();
  // Make sure that _electrostaticPotential contains the current electrostatic potential.
  void getData();
  // Write the current electrostatic potential to disk.
  void toHDF5();
  // Load the electrostatic potential from disk.
  void fromHDF5();
};

} /* namespace Serenity */

#endif /* DATA_GRID_ELECTROSTATICPOTENTIALONGRIDCONTROLLER_H_ */
