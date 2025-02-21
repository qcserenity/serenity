/**
 * @file Besley.h
 *
 * @date Dec 06, 2018
 * @author Johannes Toelle, Niklas Niemeyer
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

#ifndef LRSCF_BESLEY
#define LRSCF_BESLEY

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {
/* Forward declaration */
class SystemController;

/**
 * @class Besley
 * Restricts the occupied and virtual orbital space of response calculations.
 * N.A. Besley (Chemical Physics Letters 390 (2004) 124–129)
 */
template<Options::SCF_MODES SCFMode>
class Besley {
 public:
  /**
   * @brief Constructor
   * @param system the systemController
   * @param nAtoms The number of atoms to be considered by the besley restriction.
   * @param besleyCutoff Both Besley cutoff values.
   */
  Besley(std::shared_ptr<SystemController> system, unsigned int nAtoms, std::vector<double> besleyCutoff);

  ///@brief Default destructor
  virtual ~Besley() = default;

  ///@brief Gets indices of orbitals that are part of the orbital space.
  SpinPolarizedData<SCFMode, std::vector<unsigned int>> getWhiteList();

 private:
  ///@brief The systemController.
  std::shared_ptr<SystemController> _system;

  ///@brief The number of atoms.
  unsigned int _nAtoms;

  ///@brief A vector containing the occupied/virtual cutoffs.
  std::vector<double> _besleyCutoff;

  ///@brief Returns pseudo coefficient matrix.
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> getPseudoCoefficientMatrix();

}; /* class Besley */
} /* namespace Serenity */

#endif /* LRSCF_BESLEY */
