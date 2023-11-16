/**
 * @file   BasisFunctionOnGridController.h
 *
 * @date   May 7, 2014
 * @author Thomas Dresselhaus
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
#ifndef BASISFUNCTIONONGRIDCONTROLLER_H_
#define BASISFUNCTIONONGRIDCONTROLLER_H_
/* Include Serenity Internal Headers */
#include "notification/ObjectSensitiveClass.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <mutex>
#include <vector>

namespace Serenity {
/* Forward declarations */
class BasisController;
class Grid;
class GridController;
template<class T>
struct Gradient;
template<class T>
struct Hessian;

/**
 * @class BasisFunctionOnGridController BasisFunctionOnGridController.h
 * @brief Lazily calculates and caches the values of basis functions on grid points
 */
class BasisFunctionOnGridController : public ObjectSensitiveClass<Grid> {
 public:
  /**
   * @brief Holds the values of the basis functions on a block of grid points
   *
   * For efficiency reasons grid data is usually processed block-wise. The data of only a single
   * chunk is stored in here.
   *
   */
  struct BasisFunctionBlockOnGridData {
   public:
    /**
     * @param nBasisFunctions The total number of basis functions.
     * @param blockSize The number of grid points per block.
     * @param derivativeLevel The highest order of the derivative of the basis functions
     *                        which should be stored.
     *
     *    MB: This used to be defined in the .h-file. I moved the function body to the .cpp, since I
     *        do not see the necessity of defining it here.
     */
    BasisFunctionBlockOnGridData(const unsigned int nBasisFunctions, const unsigned int blockSize,
                                 const unsigned int derivativeLevel);
    /**
     * @brief Copy constructor for BasisFunctionBlockOnGridData.
     * @param orig The original.
     *
     *    MB: This used to be defined in the .h-file. I moved the function body to the .cpp, since I
     *        do not see the necessity of defining it here.
     */
    explicit BasisFunctionBlockOnGridData(const BasisFunctionBlockOnGridData& orig);
    /** @cond false */
    ///@brief The block center.
    Eigen::Vector3d center;
    ///@brief Quasi bool vector containing a 1 if the basis function is negligible and 0 if not.
    /// TODO: Change to sparse vector (1 for non-negligible) or std::vector containing
    /// only the indices of the non-negligible basis functions.
    Eigen::VectorXi negligible;
    ///@brief The basis function values.
    /// TODO: Reduce the dimension to include only non-negligible values.
    Eigen::MatrixXd functionValues;
    ///@brief Average grid values for each basis function.
    Eigen::VectorXd averageFunctionValues;
    ///@brief The values of the derivatives.
    /// TODO: Reduce the dimension to include only non-negligible values.
    std::unique_ptr<Gradient<Eigen::MatrixXd>> derivativeValues;
    ///@brief The values of the second derivative.
    /// TODO: Reduce the dimension to include only non-negligible values.
    std::unique_ptr<Hessian<Eigen::MatrixXd>> secondDerivativeValues;
    /** @endcond */
  };
  /**
   * @param basisController The basis controller for which the values are calculated.
   * @param gridController The grid controller for the Grid on which the instance will work.
   * @param maxBlockSize The maximum number of point in one block (for parallelization and averages).
   * @param radialThreshold If the radial part of a function is below this threshold the angular
   *                        part is skipped.
   * @param highestDerivative Level of the highest derivative calculated. 0 means basis function
   *                          value only, 1 means also calculate basis function gradients, 2 means
   *                          calculate values, gradients (first derivatives) and Hessian (second
   *                          derivatives).
   */
  BasisFunctionOnGridController(std::shared_ptr<BasisController> basisController,
                                std::shared_ptr<GridController> gridController, const unsigned int maxBlockSize,
                                const double radialThreshold, const unsigned int highestDerivative);
  /**
   * @brief Default destructor.
   */
  virtual ~BasisFunctionOnGridController() = default;

  /**
   * @brief Reinitialize when grid or basis changes.
   */
  void notify();
  /**
   * @brief Getter for the data associated to the given block-index.
   * @param   blockIndex       The index of the block of grid points under consideration.
   *
   * @returns The grid data for this block, which is managed by the controller. Caution: It is
   *          possible that your data become invalid after another function call (from the same
   *          thread)!
   */
  const std::shared_ptr<BasisFunctionBlockOnGridData> getBlockOnGridData(const unsigned int blockIndex);
  /**
   * @brief Getter for the point index of the first point in the given block.
   * @param blockIndex The block number for which the first index is requested.
   * @returns The first index of data objects returned by this class for the block with the given
   *          index.
   */
  unsigned int getFirstIndexOfBlock(const unsigned int blockIndex);
  /**
   * @brief Getter for the number of grid-point blocks.
   * @returns The number of blocks of grid points managed by this controller.
   */
  inline unsigned int getNBlocks() {
    return _nBlocks;
  }
  /**
   * @brief Getter for the number of grid points.
   * @returns (forwarded from GridController)
   */
  unsigned int getNGridPoints() const;
  /**
   * @brief Getter for the number of basis functins.
   * @returns (forwarded from BasisCntroller)
   */
  unsigned int getNBasisFunctions() const;
  /**
   * @brief Calculate the basis function vales/their derivatives etc..
   *
   * Math for spherical basis functions (to generate more, see: sharmonics.py):\n
   * l=0\n
   * \f{eqnarray*}{
   * Y^{ 0}_{ 0}(x,y,z) &=& \sqrt{\frac{  1}{4\pi}} \sqrt{ \frac{1}{1} } (1.0) \\
   * \f}
   * l=1\n
   * \f{eqnarray*}{
   * Y^{-1}_{ 1}(x,y,z) &=& \sqrt{\frac{  3}{4\pi}} \sqrt{ \frac{1}{1} } (1.0y^1) \\
   * Y^{ 0}_{ 1}(x,y,z) &=& \sqrt{\frac{  3}{4\pi}} \sqrt{ \frac{1}{1} } (1.0z^1) \\
   * Y^{ 1}_{ 1}(x,y,z) &=& \sqrt{\frac{  3}{4\pi}} \sqrt{ \frac{1}{1} } (1.0x^1) \\
   * \f}
   * l=2\n
   * \f{eqnarray*}{
   * Y^{-2}_{ 2}(x,y,z) &=& \sqrt{\frac{  5}{4\pi}} \sqrt{ \frac{3}{4} } (2.0x^1y^1) \\
   * Y^{-1}_{ 2}(x,y,z) &=& \sqrt{\frac{  5}{4\pi}} \sqrt{ \frac{3}{1} } (1.0y^1z^1) \\
   * Y^{ 0}_{ 2}(x,y,z) &=& \sqrt{\frac{  5}{4\pi}} \sqrt{ \frac{1}{4} } (2.0z^2) + (-1.0x^2) + (-1.0y^2) \\
   * Y^{ 1}_{ 2}(x,y,z) &=& \sqrt{\frac{  5}{4\pi}} \sqrt{ \frac{3}{1} } (1.0x^1z^1) \\
   * Y^{ 2}_{ 2}(x,y,z) &=& \sqrt{\frac{  5}{4\pi}} \sqrt{ \frac{3}{4} } (1.0x^2) + (-1.0y^2) \\
   * \f}
   * l=3\n
   * \f{eqnarray*}{
   * Y^{-3}_{ 3}(x,y,z) &=& \sqrt{\frac{  7}{4\pi}} \sqrt{ \frac{5}{8} } (3.0x^2y^1) + (-1.0y^3) \\
   * Y^{-2}_{ 3}(x,y,z) &=& \sqrt{\frac{  7}{4\pi}} \sqrt{ \frac{15}{4} } (2.0x^1y^1z^1) \\
   * Y^{-1}_{ 3}(x,y,z) &=& \sqrt{\frac{  7}{4\pi}} \sqrt{ \frac{3}{8} } (4.0y^1z^2) + (-1.0x^2y^1) + (-1.0y^3) \\
   * Y^{ 0}_{ 3}(x,y,z) &=& \sqrt{\frac{  7}{4\pi}} \sqrt{ \frac{1}{4} } (2.0z^3) + (-3.0x^2z^1) + (-3.0y^2z^1) \\
   * Y^{ 1}_{ 3}(x,y,z) &=& \sqrt{\frac{  7}{4\pi}} \sqrt{ \frac{3}{8} } (4.0x^1z^2) + (-1.0x^3) + (-1.0x^1y^2) \\
   * Y^{ 2}_{ 3}(x,y,z) &=& \sqrt{\frac{  7}{4\pi}} \sqrt{ \frac{15}{4} } (1.0x^2z^1) + (-1.0y^2z^1) \\
   * Y^{ 3}_{ 3}(x,y,z) &=& \sqrt{\frac{  7}{4\pi}} \sqrt{ \frac{5}{8} } (1.0x^3) + (-3.0x^1y^2) \\
   * \f}
   * l=4\n
   * \f{eqnarray*}{
   * Y^{-4}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{35}{64} } (4.0x^3y^1) + (-4.0x^1y^3) \\
   * Y^{-3}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{35}{8} } (3.0x^2y^1z^1) + (-1.0y^3z^1) \\
   * Y^{-2}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{5}{16} } (12.0x^1y^1z^2) + (-2.0x^3y^1) + (-2.0x^1y^3)
   * \\
   * Y^{-1}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{5}{8} } (4.0y^1z^3) + (-3.0x^2y^1z^1) + (-3.0y^3z^1) \\
   * Y^{ 0}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{1}{64} } (8.0z^4) + (-24.0x^2z^2) + (-24.0y^2z^2) +
   * (3.0x^4) + (3.0y^4) + (6.0x^2y^2) \\
   * Y^{ 1}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{5}{8} } (4.0x^1z^3) + (-3.0x^3z^1) + (-3.0x^1y^2z^1) \\
   * Y^{ 2}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{5}{16} } (6.0x^2z^2) + (-1.0x^4) + (-6.0y^2z^2) +
   * (1.0y^4) \\
   * Y^{ 3}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{35}{8} } (1.0x^3z^1) + (-3.0x^1y^2z^1) \\
   * Y^{ 4}_{ 4}(x,y,z) &=& \sqrt{\frac{  9}{4\pi}} \sqrt{ \frac{35}{64} } (1.0x^4) + (-6.0x^2y^2) + (1.0y^4) \\
   * \f}
   * l=5\n
   * \f{eqnarray*}{
   * Y^{-5}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{63}{128} } (5.0x^4y^1) + (-10.0x^2y^3) + (1.0y^5) \\
   * Y^{-4}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{315}{64} } (4.0x^3y^1z^1) + (-4.0x^1y^3z^1) \\
   * Y^{-3}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{35}{128} } (24.0x^2y^1z^2) + (-3.0x^4y^1) +
   * (-2.0x^2y^3) + (-8.0y^3z^2) + (1.0y^5) \\
   * Y^{-2}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{105}{16} } (4.0x^1y^1z^3) + (-2.0x^3y^1z^1) +
   * (-2.0x^1y^3z^1) \\
   * Y^{-1}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{15}{64} } (8.0y^1z^4) + (-12.0x^2y^1z^2) +
   * (-12.0y^3z^2) + (1.0x^4y^1) + (1.0y^5) + (2.0x^2y^3) \\
   * Y^{ 0}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{1}{64} } (8.0z^5) + (-40.0x^2z^3) + (-40.0y^2z^3) +
   * (15.0x^4z^1) + (15.0y^4z^1) + (30.0x^2y^2z^1) \\
   * Y^{ 1}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{15}{64} } (8.0x^1z^4) + (-12.0x^3z^2) +
   * (-12.0x^1y^2z^2) + (1.0x^5) + (1.0x^1y^4) + (2.0x^3y^2) \\
   * Y^{ 2}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{105}{16} } (2.0x^2z^3) + (-1.0x^4z^1) + (-2.0y^2z^3) +
   * (1.0y^4z^1) \\
   * Y^{ 3}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{35}{128} } (8.0x^3z^2) + (-1.0x^5) + (2.0x^3y^2) +
   * (-24.0x^1y^2z^2) + (3.0x^1y^4) \\
   * Y^{ 4}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{315}{64} } (1.0x^4z^1) + (-6.0x^2y^2z^1) + (1.0y^4z^1)
   * \\
   * Y^{ 5}_{ 5}(x,y,z) &=& \sqrt{\frac{ 11}{4\pi}} \sqrt{ \frac{63}{128} } (1.0x^5) + (-10.0x^3y^2) + (5.0x^1y^4) \\
   * \f}
   * l=6\n
   * \f{eqnarray*}{
   * Y^{-6}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{231}{512} } (6.0x^5y^1) + (-20.0x^3y^3) + (6.0x^1y^5)
   * \\
   * Y^{-5}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{693}{128} } (5.0x^4y^1z^1) + (-10.0x^2y^3z^1) +
   * (1.0y^5z^1) \\
   * Y^{-4}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{63}{256} } (40.0x^3y^1z^2) + (-4.0x^5y^1) +
   * (-40.0x^1y^3z^2) + (4.0x^1y^5) \\
   * Y^{-3}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{105}{128} } (24.0x^2y^1z^3) + (-9.0x^4y^1z^1) +
   * (-6.0x^2y^3z^1) + (-8.0y^3z^3) + (3.0y^5z^1) \\
   * Y^{-2}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{105}{512} } (32.0x^1y^1z^4) + (-32.0x^3y^1z^2) +
   * (-32.0x^1y^3z^2) + (2.0x^5y^1) + (2.0x^1y^5) + (4.0x^3y^3) \\
   * Y^{-1}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{21}{64} } (8.0y^1z^5) + (-20.0x^2y^1z^3) +
   * (-20.0y^3z^3) + (5.0x^4y^1z^1) + (5.0y^5z^1) + (10.0x^2y^3z^1) \\
   * Y^{ 0}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{1}{256} } (331.0z^6) + (-120.0x^2z^4) + (-120.0y^2z^4)
   * + (-315.0z^6) + (90.0x^4z^2) + (90.0y^4z^2) + (180.0x^2y^2z^2) + (-5.0x^6) + (-5.0y^6) + (-15.0x^4y^2) +
   * (-15.0x^2y^4) \\
   * Y^{ 1}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{21}{64} } (8.0x^1z^5) + (-20.0x^3z^3) +
   * (-20.0x^1y^2z^3) + (5.0x^5z^1) + (5.0x^1y^4z^1) + (10.0x^3y^2z^1) \\
   * Y^{ 2}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{105}{512} } (16.0x^2z^4) + (-16.0x^4z^2) + (1.0x^6) +
   * (-1.0x^2y^4) + (1.0x^4y^2) + (-34.0y^2z^4) + (16.0y^4z^2) + (18.0y^2z^4) + (-1.0y^6) \\
   * Y^{ 3}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{105}{128} } (8.0x^3z^3) + (-3.0x^5z^1) + (6.0x^3y^2z^1)
   * + (-24.0x^1y^2z^3) + (9.0x^1y^4z^1) \\
   * Y^{ 4}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{63}{256} } (10.0x^4z^2) + (-1.0x^6) + (5.0x^4y^2) +
   * (-60.0x^2y^2z^2) + (5.0x^2y^4) + (10.0y^4z^2) + (-1.0y^6) \\
   * Y^{ 5}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{693}{128} } (1.0x^5z^1) + (-10.0x^3y^2z^1) +
   * (5.0x^1y^4z^1) \\
   * Y^{ 6}_{ 6}(x,y,z) &=& \sqrt{\frac{ 13}{4\pi}} \sqrt{ \frac{231}{512} } (1.0x^6) + (-15.0x^4y^2) + (15.0x^2y^4) +
   * (-1.0y^6) \\ \f}
   *
   * @param blockIndex          The index of the block of grid points for which data shall be calculated.
   * @returns the data related to this block of grid points.
   */
  std::shared_ptr<BasisFunctionBlockOnGridData> calculateBasisFunctionData(const unsigned int blockIndex);
  /**
   * @brief Getter for the underlying basis controller.
   * @returns The controller for the basis, the object is working in.
   */
  inline std::shared_ptr<BasisController> getBasisController() const {
    return _basisController;
  }
  /**
   * @brief Getter for the underlying grid controller.
   * @returns The Grid(controller) the object is working on.
   */
  inline std::shared_ptr<GridController> getGridController() const {
    return _gridController;
  }
  /**
   * @brief Getter for the highest derivative of the basis functions to be calculated.
   * @returns A number indicating up to which order derivatives are calculated by this object.
   */
  unsigned int getHighestDerivative() const {
    return _highestDerivative;
  }
  /**
   * @brief sets a new highest derivarive.
   *
   * This indicates up to which order derivatives are calculated by this object.
   * TODO this solution does not look elegant.
   */
  void setHighestDerivative(unsigned int newHighestDerivative);
  /**
   * @brief Getter for the maximum block size.
   * @return The maximum block size.
   */
  unsigned int getMaxBlockSize() {
    return _maxBlockSize;
  }

 private:
  const std::shared_ptr<BasisController> _basisController;
  const std::shared_ptr<GridController> _gridController;
  const unsigned int _maxBlockSize;
  unsigned int _nPoints;
  unsigned int _nBlocks;
  const double _radialThreshold;
  const double _exponentThreshold;
  unsigned int _highestDerivative;
  bool _upToDate = false;

  /**
   * For the blocks of data which are recalculated always the same memory is used. This avoids
   * frequent creation / destruction of data objects (this would slow down the program quite a bit).
   */
  std::vector<std::shared_ptr<BasisFunctionBlockOnGridData>> _workspace;
  std::shared_ptr<BasisFunctionBlockOnGridData> _lastBlock;

  std::mutex _lock;
};

} /* namespace Serenity */

#endif /* BASISFUNCTIONONGRIDCONTROLLER_H_ */
