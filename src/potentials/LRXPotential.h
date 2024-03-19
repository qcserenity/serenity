/**
 * @file LRXPotential.h
 *
 * @date Mar 31, 2017
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

#ifndef POTENTIALS_LRXPOTENTIAL_H_
#define POTENTIALS_LRXPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "potentials/IncrementalFockMatrix.h"
#include "potentials/Potential.h"
#include "settings/Options.h"

namespace Serenity {
/* Forward declaration */
class SystemController;
/**
 * @class LRXPotential LRXPotential.h
 * @brief Calculates the matrix representation of the long range contribution to the
 *        exchange potential in range-separated exchange--correlation functionals.
 */
template<Options::SCF_MODES SCFMode>
class LRXPotential : public Potential<SCFMode>,
                     public ObjectSensitiveClass<Basis>,
                     public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param systemController The system controller.
   * @param dMat The density matrix controller.
   * @param exchangeRatio The exchange ratio.
   * @param prescreeningThreshold The schwartz prescreening threshold.
   * @param prescreeningIncrementStart The start integrals prescreening thresold for the incremental Fock-matrix build
   * @param prescreeningIncrementEnd The end integrals prescreening thresold for the incremental Fock-matrix build
   * @param incrementSteps The number of steps of an incremental Fock-matrix build until it gets rebuild
   * @param mu The range seperation parameter.
   */
  LRXPotential(std::shared_ptr<SystemController> systemController, std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
               const double exchangeRatio, const double prescreeningThreshold, double prescreeningIncrementStart,
               double prescreeningIncrementEnd, unsigned int incrementSteps, const double mu);

  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Adds an increment to an existing potential in matrix form.
   *
   * This function will prescreen based on the delta-density matrix.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP An increment of the density matrix.
   */
  void addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& deltaP);

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Potential is linked to the basis it is defined in.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _outOfDate = true;
  };

  virtual ~LRXPotential() = default;

 private:
  ///@brief The system
  std::weak_ptr<SystemController> _systemController;
  ///@brief The exchange ratio.
  const double _exc;
  ///@brief Density matrix controller for this potential
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief The entire potential
  std::shared_ptr<FockMatrix<SCFMode>> _fullpotential;
  ///@brief The range separation parameter
  const double _mu;
  ///@brief Checks if the data is up to date
  bool _outOfDate;
  ///@brief Screening Threshold for the current iteration
  double _screening;
  // Helper for the incremental fock matrix construction.
  std::shared_ptr<IncrementalFockMatrix<SCFMode>> _incrementHelper;
};

} /* namespace Serenity */

#endif /* POTENTIALS_LRXPOTENTIAL_H_ */
