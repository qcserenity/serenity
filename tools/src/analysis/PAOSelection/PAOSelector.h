/**
 * @file PAOSelector.h
 *
 * @date Jan 25, 2019
 * @author Moritz Bensberg
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

#ifndef ANALYSIS_PAOSELECTION_PAOSELECTOR_H_
#define ANALYSIS_PAOSELECTION_PAOSELECTOR_H_

/* Include Std and External Headers */
#include <Eigen/SparseCore> //Sparse matrices.
#include <memory>           //smrt_ptr.
namespace Serenity {
/**
 * @class PAOSelector PAOSelector.h
 * @brief An interface for different PAO selection algorithms. For general information
 *        on PAOs @see data/PAOController.h
 */
class PAOSelector {
 public:
  /**
   * @brief Default constructor.
   */
  PAOSelector() = default;
  /**
   * @brief Default destructor.
   */
  virtual ~PAOSelector() = default;
  /**
   * @brief Selects the PAOs.
   * @return The PAO selection for each occupied orbital.
   */
  virtual std::shared_ptr<Eigen::SparseMatrix<int>> selectPAOs() = 0;
};

} /* namespace Serenity */

#endif /* ANALYSIS_PAOSELECTION_PAOSELECTOR_H_ */
