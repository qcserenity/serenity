/**
 * @file DeltaSpinSquared.h
 *
 * @date Jul. 7, 2021
 * @author Johannes Toelle
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

#ifndef LRSCF_DELTASPINSQUARED
#define LRSCF_DELTASPINSQUARED

/* Include Serenity Internal Headers */
#include "settings/ElectronicStructureOptions.h"
#include "settings/LRSCFOptions.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
class LRSCFController;

/**
 * @class DeltaSpinSquared
 * @brief Computes the Delta S^2 expectation value for a given transition.
 * A TDDFT/TDA calculation needs to be performed beforehand.
 * Formulas are taken from: J. Chem. Phys.134, 134101 (2011) (mistakes are corrected)
 */
template<Options::SCF_MODES SCFMode>
class DeltaSpinSquared {
 public:
  /**
   * @brief Default destructor.
   */
  virtual ~DeltaSpinSquared() = default;

  /**
   * @brief Constructor.
   * @param lrscf The LRSCF Controller.
   * @param nEigen The number of eigenvalues.
   * @param method Linear response method, here only TDA or TDDFT.
   * @param type Isolated, uncoupled, coupled.
   */
  DeltaSpinSquared(std::vector<std::shared_ptr<LRSCFController<SCFMode>>> lrscf, unsigned int nEigen,
                   Options::LR_METHOD method, Options::LRSCF_TYPE type);

  /**
   * @brief Performs the calculation for a specific state.
   * @param iState The state index.
   */
  double calculateSpinSquared(unsigned int iState);

  /**
   * @brief Prints calculated Delta S^2 values to screen.
   */
  void print();

 private:
  std::vector<std::shared_ptr<LRSCFController<SCFMode>>> _lrscf;
  unsigned int _nEigen;
  Eigen::VectorXd _spinSquared;
  Options::LR_METHOD _method;
  Options::LRSCF_TYPE _type;
};

} /* namespace Serenity */
#endif /* LRSCF_DELTASPINSQUARED */
