/**
 * @file LRSCFAnalysis.h
 *
 * @date Dec 1, 2016
 * @author M. Boeckers
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

#ifndef ELECTRONICSTRUCTURECALCULATIONS_LRSCF_LRSCFANALYSIS_H_
#define ELECTRONICSTRUCTURECALCULATIONS_LRSCF_LRSCFANALYSIS_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "system/SystemController.h"
/* Include Std and External Headers */
#include <Eigen/Dense>


namespace Serenity {


/**
 * @class LRSCFAnalysis LRSCFAnalysis.h
 */
template<Options::SCF_MODES T>
class LRSCFAnalysis {
public:
  LRSCFAnalysis(
      std::shared_ptr<SystemController> activeSystem,
      std::vector<std::shared_ptr<SystemController> > environmentSystems,
      Eigen::VectorXd& eigenvalues,
      std::vector<Eigen::MatrixXd >& eigenvectors);
  virtual ~LRSCFAnalysis() = default;

  /**
   * @brief Within LRSCF problems, the X and Y coefficients are normalized as
   *        \f[
   *        sum_{ia} X_{ia}^2 - Y_{ia}^2 = 1 \; .
   *        \f]
   *        The squared coefficient (given by ORCA or TURBOMOLE in LRSCF calculations)
   *        can be defined as
   *        \f[
   *        |c_{ia}|^2 = X_{ia}^2 - Y_{ia}^2  \;.
   *        \f]
   *        This function prints |c_{ia}|^2 * 100 for the largest coefficients.
   */
  void printStateInfo(const unsigned int iState);

  /**
   * @brief Prints the excitation vectors in AO basis for analysis. The vectors are written to
   *        .../path/lrscf_type.X and .../path/lrscf_type.Y, respectively. For each element
   *        of the excitation vectors, a unique basis function identifier is given. This function
   *        can be used to e.g. compare excitation vectors of super and subsystem calculation
   *        using the same atomic basis sets.
   */
  void printAOExcitationVectors();


  /**
   * @brief Perform Mulliken population analysis for state iState
   * @param iState
   */
  void mullikenPopulationAnalysis(const unsigned int iState);


private:
  //The system controller of the active subsystem
  std::shared_ptr<SystemController> _activeSystem;

  //System controller of environment subsystems
  std::vector<std::shared_ptr<SystemController> > _environmentSystems;

  //A vector holding the excitation energies
  Eigen::VectorXd _eigenvalues;

  //A vector holding the CI coefficients
  std::vector<Eigen::MatrixXd> _eigenvectors;

  //The number of occupied orbitals
  SpinPolarizedData<T,unsigned int> _nOccupied;

  //Number of virtual orbitals in the active system

  SpinPolarizedData<T, unsigned int> _nVirtual;

  //The number of MOs in the active system
  unsigned int _nMolecularOrbitals;

  Eigen::MatrixXd _x;
  Eigen::MatrixXd _y;

  Eigen::VectorXd Mo2Ao(Eigen::VectorXd moVec);

  void writeExcitationVectors(Eigen::MatrixXd& aoVecs, std::string fname);
;
};

} /* namespace Serenity */

#endif /* ELECTRONICSTRUCTURECALCULATIONS_LRSCF_LRSCFANALYSIS_H_ */
