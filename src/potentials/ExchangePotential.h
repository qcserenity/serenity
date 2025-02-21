/**
 * @file ExchangePotential.h
 *
 * @author Moritz Bensberg
 * @date Sep 23, 2019
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

#ifndef POTENTIALS_EXCHANGEPOTENTIAL_H_
#define POTENTIALS_EXCHANGEPOTENTIAL_H_
/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrixController.h"
#include "potentials/IncrementalFockMatrix.h"
#include "potentials/Potential.h"

namespace Serenity {
/* Forward declaration */
class SystemController;
/**
 * @class ExchangePotential ExchangePotential.h
 * @brief A class that calculates the exchange contribution to the Fock matrix using
 *        four center, two electron integrals.
 */
template<Options::SCF_MODES SCFMode>
class ExchangePotential : public Potential<SCFMode>,
                          public ObjectSensitiveClass<Basis>,
                          public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor.
   * @param systemController The system controller.
   * @param dMat The density matrix controller of the system.
   * @param exchangeRatio The scaling ratio of the exchange operator.
   * @param prescreeningThreshold The Schwarz prescreening threshold.
   * @param prescreeningIncrementStart The start integrals prescreening thresold for the incremental Fock-matrix build
   * @param prescreeningIncrementEnd The end integrals prescreening thresold for the incremental Fock-matrix build
   * @param incrementSteps The number of steps of an incremental Fock-matrix build until it gets rebuild
   * @param clear4CenterCache          If true, the 4-center integral cache is deleted upon destruction of the
   *                                   potential.
   */
  ExchangePotential(std::shared_ptr<SystemController> systemController,
                    std::shared_ptr<DensityMatrixController<SCFMode>> dMat, const double exchangeRatio,
                    const double prescreeningThreshold, double prescreeningIncrementStart, double prescreeningIncrementEnd,
                    unsigned int incrementSteps, bool clear4CenterCache = true, bool transitionDensity = false);
  /**
   * @brief Destructor. Clears the 4-center integral cache.
   */
  virtual ~ExchangePotential();

  /**
   * @brief Getter for the actual potential.
   *
   * @return Returns the potential in matrix form.
   */
  virtual FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @brief Calculates an incremental potential.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP An increment of the density matrix.
   */
  void addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& deltaP);

  /**
   * @brief Calculates an incremental potential for a non-symmetric transition density.
   *
   * @param F      A reference to the potential to be added to.
   * @param deltaP An increment of the density matrix.
   */
  void addToMatrixForTransitionDensity(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& deltaP);

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The active density Matrix.
   * @return The energy of the potential when acting on this density.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
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

 private:
  ///@brief The underlying systemController
  std::weak_ptr<SystemController> _systemController;
  ///@brief The exchange ratio.
  const double _exc;
  ///@brief The basis this potential is defined in.
  std::shared_ptr<DensityMatrixController<SCFMode>> _dMatController;
  ///@brief The entire potential
  std::shared_ptr<FockMatrix<SCFMode>> _fullpotential;
  ///@brief Checks if the data is up to date
  bool _outOfDate;
  ///@brief Screening Threshold for the current iteration
  double _screening;
  ///@brief Internal iteration counter
  unsigned int _counter = 0;
  /// @brief Helper for the incremental fock matrix construction.
  std::shared_ptr<IncrementalFockMatrix<SCFMode>> _incrementHelper;
  /// @brief If true, the 4-center integral cache of the system is deleted upon potential destruction.
  bool _clear4CenterCache;
  bool _transitionDensity;
};

} /* namespace Serenity */

#endif /* POTENTIALS_EXCHANGEPOTENTIAL_H_ */
