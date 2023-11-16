/**
 * @file CoulombIntegralsOnGridController.h
 *
 * @author Moritz Bensberg
 * @date Aug 11, 2020
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

#ifndef DATA_GRID_COULOMBINTEGRALSONGRIDCONTROLLER_H_
#define DATA_GRID_COULOMBINTEGRALSONGRIDCONTROLLER_H_

/* Include Serenity Internal Headers */
#include "basis/Basis.h"               //Loop basis shells.
#include "basis/BasisController.h"     //The basis controller. Shell/basis infos.
#include "grid/GridController.h"       //The grid coordinates.
#include "integrals/wrappers/Libint.h" //Calculate integrals.
#include "io/HDF5.h"                   //Write to disk.
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace Serenity {
/**
 * @class
 * @brief A controller for two center integrals on a grid of type:
 *      \f$(\mu|1/(r-s)|\nu) \f$,  \f$ \nabla (\mu|1/(r-s)|\nu) \f$ and  \f$ (\mu|n_s(r-s)/|r-s|^3|\nu) \f$.
 *      For the integrals \f$(\mu|1/(r-s)|\nu) \f$ extensive caching on disk is used.
 */
class CoulombIntegralsOnGridController : public ObjectSensitiveClass<Grid> {
 public:
  /**
   * @brief Constructor.
   * @param basisController   The basis controller.
   * @param gridController    The grid controller
   * @param fBaseName         Base name for the cache files.
   * @param cacheSize         The maximum number of two center sets hold in memory at the same time.
   * @param normalVectors     The normal vectors at the grid points (optional).
   */
  CoulombIntegralsOnGridController(std::shared_ptr<BasisController> basisController,
                                   std::shared_ptr<GridController> gridController, std::string fBaseName,
                                   unsigned int cacheSize, std::shared_ptr<Eigen::Matrix3Xd> normalVectors = nullptr);
  virtual ~CoulombIntegralsOnGridController() = default;
  /**
   * @brief Remove all cached files.
   */
  void cleanUp();
  /**
   * @brief Save the integrals on file.
   * @param integrals The integrals.
   * @param blockIndex The grid point block index.
   * @param file The HDF5 file.
   */
  void integralsToHDF5(std::vector<std::shared_ptr<Eigen::MatrixXd>> integrals, unsigned int blockIndex, HDF5::H5File& file);
  /**
   * @brief Read a set of integrals from file.
   * @param blockIndex The grid point block index.
   * @param file The HDF5 file.
   * @return The integrals.
   */
  std::vector<std::shared_ptr<Eigen::MatrixXd>> integralsFromHDF5(unsigned int blockIndex, HDF5::H5File& file);
  /**
   * @brief Calculate the integrals for the given set of grid coordinates.
   * @param coordinates The grid point coordinates.
   * @param normalVectors The grid point normal vectors.
   * @return The set of integrals.
   */
  std::vector<std::shared_ptr<Eigen::MatrixXd>> calculateIntegralSet(const Eigen::Matrix3Xd& coordinates,
                                                                     const Eigen::Matrix3Xd& normalVectors);
  /**
   * @brief Getter for the basis controller.
   */
  std::shared_ptr<BasisController> getBasisController();
  /**
   * @brief Getter for the grid controller.
   */
  std::shared_ptr<GridController> getGridController();

  /**
   * @brief Looper function for the two center integrals at the grid points.
   *        Uses chaching.
   * @param distributionFunction The function that is executed for every integral set.
   *
   *  The function has to accept the following parameters:
   *    distributionFunction(
   *    std::shared_ptr<Eigen::MatrixXd> ints,    //All two center integrals for the given grid point i.
   *     unsigned int i,                          //The grid point index.
   *      unsigned int threadId);                 //The thread id.
   */
  template<class Func>
  __attribute__((always_inline)) inline void iterateGridPoints(Func distributionFunction) {
    const unsigned int nTotalSets = _gridController->getNGridPoints();
    if (nTotalSets < _nSetsCached)
      _nSetsCached = nTotalSets;
    if (_cache.size() == 0)
      initializeIntegrals();
    const unsigned int nBlocks = _blockSizes.size();
    unsigned int blockEnd = nTotalSets;
    HDF5::Filepath name(_intFileBaseName);
    HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
    Eigen::setNbThreads(1);
    for (int iBlock = nBlocks - 1; iBlock >= 0; --iBlock) {
      unsigned int n = _blockSizes[iBlock];
      unsigned int start = blockEnd - n;
      if (!_cache[start]) {
        auto ints = integralsFromHDF5(iBlock, file);
        for (unsigned int i = start; i < blockEnd; ++i)
          _cache[i] = ints[i - start];
      }
#pragma omp for schedule(dynamic)
      for (unsigned int i = start; i < blockEnd; ++i) {
        distributionFunction(_cache[i], i, omp_get_thread_num());
        if (nBlocks > 1)
          _cache[i] = nullptr;
      }
      blockEnd = start;
    }
    Eigen::setNbThreads(0);
    file.close();
  }
  // clang-format off
  /**
   * @brief Looper function for the integrals \f$ \nabla (\mu|1/(r-s)|\nu) \f$.
   *        Does not use any caching.
   * @param distributionFunction The function applied to the integral sets.
   *
   * The function accepts the following arguments:
   *   distributionFunction(
   *      Eigen::MatrixXd& ints,             //The integrals for a given shell combination.
   *      				     //Column-wise as: (dx mu|nu),...(mu|dx nu),... (mu|dx 1/(r-x)|nu)...
   *       unsigned int iPoint,              //The point index.
   *       unsigned int offI,                //The basis function index of the first function in shell I.
   *       unsigned int offJ,                //The basis function index of the first function in shell J.
   *       unsigned int nI,                  //The number of basis functions in shell I.
   *       unsigned int nJ,                  //The number of basis functions in shell J.
   *       unsigned int omp_get_thread_num)  //The thread id.
   */
  // clang-format on
  template<class Func>
  __attribute__((always_inline)) inline void iterateGridPoints_gradient(Func distributionFunction) {
    const unsigned int nPoints = _gridController->getNGridPoints();
    const Eigen::Matrix3Xd& coordinates = _gridController->getGridPoints();
    Eigen::setNbThreads(1);
    auto& libint = Libint::getInstance();
    auto basisController = this->getBasisController();
    const auto& bas = basisController->getBasis();
#pragma omp for schedule(dynamic)
    for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint) {
      std::vector<std::pair<double, std::array<double, 3>>> point = {
          {-1.0, {{coordinates(0, iPoint), coordinates(1, iPoint), coordinates(2, iPoint)}}}};
      Eigen::MatrixXd ints;
      libint.initialize(LIBINT_OPERATOR::nuclear, 1, 2, point, 0.0);
      for (unsigned int i = 0; i < bas.size(); ++i) {
        unsigned int offI = basisController->extendedIndex(i);
        const unsigned int nI = bas[i]->getNContracted();
        for (unsigned int j = 0; j <= i; ++j) {
          unsigned int offJ = basisController->extendedIndex(j);
          const unsigned int nJ = bas[j]->getNContracted();
          if (libint.compute(LIBINT_OPERATOR::nuclear, 1, *bas[i], *bas[j], ints)) {
            distributionFunction(ints, iPoint, offI, offJ, nI, nJ, omp_get_thread_num());
          } // if compute
        }   // for j
      }     // for i
      libint.finalize(LIBINT_OPERATOR::nuclear, 1, 2);
    }
    Eigen::setNbThreads(0);
  }
  /**
   * @brief Remove all cached data.
   */
  void notify();

 private:
  // Calculate the two center integrals for the sed of grid coordinates.
  void initializeIntegrals();
  // Calculate integrals of type (mu |1/(r-s)| nu)
  Eigen::MatrixXd calculatePotentialIntegrals(std::vector<std::pair<double, std::array<double, 3>>> point);
  // Calculate integrals of type \f$ (\mu|n_s(r-s)/|r-s|^3|\nu) \f$ numerically.
  void calculateElectricFieldIntegralsNumerically(Eigen::MatrixXd& result,
                                                  std::vector<std::pair<double, std::array<double, 3>>> point,
                                                  const Eigen::Vector3d& normalVector);
  // The basis controller
  std::shared_ptr<BasisController> _basisController;
  // The grid controller
  std::shared_ptr<GridController> _gridController;
  // The base name for the two center integral file.
  std::string _intFileBaseName;
  // Integral cache. For each grid point a set of two center integrals needs to be calculated and saved.
  // This vector will contain the integrals or nullptr, if they have to be loaded from disk.
  std::vector<std::shared_ptr<Eigen::MatrixXd>> _cache;
  // The maximum number of two electron center sets to be stored in memory.
  // Contraction of the integrals with the density matrix is done block wise.
  unsigned int _nSetsCached;
  // The normal vectos n for the electric field integrals.
  std::shared_ptr<Eigen::Matrix3Xd> _normalVectors;
  // The number of two center integral sets in each block.
  // Contraction of the integrals with the density matrix is done block wise.
  std::vector<unsigned int> _blockSizes;
  // Sizes of the cached matrices.
  bool _sizesInitialized = false;
  unsigned int _nRows;
  unsigned int _nCols;
};

} /* namespace Serenity */

#endif /* DATA_GRID_COULOMBINTEGRALSONGRIDCONTROLLER_H_ */
