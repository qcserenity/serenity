/**
 * @file DensityMatrixController.h
 *
 * @date May 12, 2016
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

#ifndef BASICS_ELECTRONICSTRUCTUREDATA_MATRICES_DENSITYMATRIXCONTROLLER_H_
#define BASICS_ELECTRONICSTRUCTUREDATA_MATRICES_DENSITYMATRIXCONTROLLER_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "notification/NotifyingClass.h"

namespace Serenity {

/* Forward declarations */
template<Options::SCF_MODES T>
class OrbitalController;
// template<class T>class ObjectSensitiveClass;

/**
 * @class DensityMatrixController DensityMatrixController.h
 * @brief A Controller for the DensityMatrix
 *
 * This object allows to track updates of the DensityMatrix
 *   and also handish alterations of the DensityMatrix that can be tracked
 *   by other parts of the code.
 */
template<Options::SCF_MODES SCFMode>
class DensityMatrixController : public NotifyingClass<DensityMatrix<SCFMode>>,
                                public ObjectSensitiveClass<OrbitalController<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param molecularOrbitals
   * @param nOccupiedOrbitals
   */
  DensityMatrixController(const std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                          const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals);
  /**
   * @brief Constructor.
   * @param molecularOrbitals
   * @param occupations
   */
  DensityMatrixController(const std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                          const SpinPolarizedData<SCFMode, Eigen::VectorXd>& occupations);
  /**
   * @brief Constructor for a single DensityMatrix.
   *        This constructor is meant to be used for a density matrix that has been generated by hand.
   *        Orbitals can be attached later on.
   * @param densityMatrix
   */
  DensityMatrixController(const DensityMatrix<SCFMode>& densityMatrix);
  /**
   * @brief Constructor from HDF5 file.
   * @param fBaseName The basename of the HDF5 files.
   * @param basisController The BasisController of the current System.
   */
  DensityMatrixController(std::string fBaseName, std::shared_ptr<BasisController> basisController, std::string id);

  virtual ~DensityMatrixController() = default;

  /// @returns Returns the DensityMatrix.
  DensityMatrix<SCFMode> getDensityMatrix();

  /// @brief Updates the DensityMatrix based on the linked MOs.
  void updateDensityMatrix();

  ///@brief Notification system implementation.
  void notify();

  /**
   * @brief A function to attach orbitals to the DensityMatrixController after its creation.
   * @param molecularOrbitals
   * @param nOccupiedOrbitals
   * @param update  Switch to disable the automatic DenstyMatrix update that would occur when calling this function.
   *                The idea behind this is that maybe the initial density guess is better than the one generated from
   *                the attached orbitals.
   */
  void attachOrbitals(std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                      const SpinPolarizedData<SCFMode, unsigned int>& nOccupiedOrbitals, bool update = true);

  /**
   * @brief A function to attach orbitals to the DensityMatrixController after its creation.
   * @param molecularOrbitals
   * @param occupations
   * @param update  Switch to disable the automatic DenstyMatrix update that would occur when calling this function.
   *                The idea behind this is that maybe the initial density guess is better than the one generated from
   *                the attached orbitals.
   */
  void attachOrbitals(std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                      const SpinPolarizedData<SCFMode, Eigen::VectorXd>& occupations, bool update = true);

  /**
   * @brief Replaces the current DensityMatrix.
   * @param The new DensityMatrix (will be copied)
   */
  void setDensityMatrix(const DensityMatrix<SCFMode>& densityMatrix);
  /**
   * @param aufbau
   * @return Returns the occupation vectors of this density matrix.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> getOccupations(bool aufbau = false);

  /**
   * @brief Sets the mode for the denity matrix storage.
   *
   * @param diskmode     Iff false keeps data in memory, else writes and reads to disk.
   * @param fBaseName    The base name for HDF5 files. Usually the systems path and systems name.
   *                     (This can be a dummy string is diskmode is set to be off.)
   * @param id           The system id. (This can be a dummy string is diskmode is set to be off.)
   */
  void setDiskMode(bool diskmode, std::string fBaseName, std::string id);

  /**
   * @brief If fractional occupations are to be used, this threshold can be used to fractionally
   * occupy orbitals with energies in this range.
   *
   * @param thresh The degeneracy threshold. 0.1 a.u. is a good suggestion here. Defaults to 0.0
   * via the system settings, leading to no fractional degeneracy.
   */
  void setDegeneracyThreshold(double thresh) {
    _degeneracyThreshold = thresh;
  }

  /**
   * @brief Prints the matrix into a file
   * @param fBaseName  The bas name of the HDF5 file the matrix shall be save to.
   */
  void toHDF5(std::string fBaseName, std::string id);
  /**
   * @brief Reads a DensityMatrix from file
   * @param fBaseName The basename of the HDF5 files.
   */
  void fromHDF5(std::string fBaseName, std::string id);

 private:
  std::unique_ptr<DensityMatrix<SCFMode>> _densityMatrix;
  std::shared_ptr<OrbitalController<SCFMode>> _molecularOrbitals;
  std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>> _occupations;
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _aufbauOccupations;
  const std::shared_ptr<BasisController> _basisController;
  bool _outOfDate;
  std::string _fBaseName;
  std::string _id;
  bool _diskmode = false;
  double _degeneracyThreshold = 0.0;

  void rebuildAufbauOccupations();
};

} /* namespace Serenity */

#endif /* BASICS_ELECTRONICSTRUCTUREDATA_MATRICES_DENSITYMATRIXCONTROLLER_H_ */