/**
 * @file OptEffPotential.h
 *
 * @date Dec 12, 2016
 * @author David Schnieders
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

#ifndef POTENTIALS_OPTEFFPOTENTIAL_H_
#define POTENTIALS_OPTEFFPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/grid/DensityOnGrid.h"
#include "data/grid/GridPotential.h"
#include "data/matrices/FockMatrix.h"
#include "math/Derivatives.h"
#include "potentials/Potential.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>

namespace Serenity {

class BasisFunctionOnGridController;
class OneElectronIntegralController;
template<Options::SCF_MODES SCFMode>
class DensityOnGridCalculator;
template<Options::SCF_MODES SCFMode>
class OrbitalController;

/**
 * @class  OptEffPotential OptEffPotential.h
 * @brief  Performs a potential reconstruction.
 */
template<Options::SCF_MODES SCFMode>
class OptEffPotential : public Potential<SCFMode> {
 public:
  /**
   * @brief Constructor
   * @param basFuncOnGridController Controller for basis functions on a grid
   * @param potBasFuncOnGridController Controller for the basis functions used for the Wu-Yang reconstruction
   * @param oneEIntController Controller for one electron integrals
   * @param densOnGridCalculator Calculator for a real space representation of
   *        the density on a grid
   * @param nOccOrbs Number of occupied orbitals
   * @param smoothFactor Smoothing constraint as used in
   *        Heaton-Burgess, T.; Bulat, F. A.; Yang, W. Phys. Rev. Lett. 98, 256401 (2007)
   *        (Only used in Wu-Yang)
   * @param singValThreshold Threshold for the singular value decomposition in potential reconstruction.
   * @param exc Amount of exact exchange
   */
  OptEffPotential(std::shared_ptr<BasisFunctionOnGridController> basFuncOnGridController,
                  std::shared_ptr<BasisFunctionOnGridController> potBasFuncOnGridController,
                  std::shared_ptr<OneElectronIntegralController> oneEIntController,
                  std::shared_ptr<DensityOnGridCalculator<SCFMode>> densOnGridCalculator,
                  const SpinPolarizedData<SCFMode, unsigned int>& nOccOrbs, const double smoothFactor = 1e-3,
                  const double singValThreshold = 1e-5, const double exc = 0.0);

  virtual ~OptEffPotential() = default;

  /**
   * @brief Getter for the potential matrix.
   * @return The potential in its matrix representation.
   */
  FockMatrix<SCFMode>& getMatrix() override final {
    assert(_potential);
    return *_potential;
  };

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The density matrix for which the energy should be
   * calculated.
   * @return The energy associated with the potential and P.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final {
    assert(_potential);
    auto& pot = *_potential;
    double energy = 0.0;
    for_spin(pot, P) {
      energy += pot_spin.cwiseProduct(P_spin).sum();
    };
    return energy;
  };

  /**
   * @brief Getter for the potential on a grid.
   * @return The reconstructed potential in its grid representation.
   */
  GridPotential<SCFMode>& getPotentialOnGrid() {
    assert(_potentialOnGrid);
    return *_potentialOnGrid;
  };

  /**
   * @brief Geometry gradient contribution from this Potential.
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Performs a Wu-Yang potential reconstruction, yielding a set of orbitals resembling
   *        a ground state electronic density as close as possible to targetDens as well
   *        as the corresponding potential.
   *
   *        For more information on the theory behind this see Ref.:
   *        Quin Wu and Weitao Yang, J. Chem. Phys. 118, 2498 (2003)
   * @param targetDens The target density to be reconstructed.
   * @param resultOrbitals The orbitals of the reconstructed system.
   * @param initialGuess An initial guess for the new potential. It is advisable to at
   *        least choose the nuclear potential here, as it can not be expressed as a
   *        linear combination of gaussian basis functions very well. Additionally,
   *        the Fermi-Amaldi potential is a popular choice.
   */
  void calculateOEP(const DensityOnGrid<SCFMode>& targetDens,
                    std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals, const FockMatrix<SCFMode>& initialGuess);

  /**
   * @brief Overload to take a potential expressed on a grid.
   *        Converts the GridPotential into a FockMatrix in
   *        the correct basis and calls the function above.
   * @param targetDens The target density to be reconstructed.
   * @param resultOrbitals The orbitals of the reconstructed system.
   * @param initialGuess An initial guess for the new potential. It is advisable to at
   *        least choose the nuclear potential here, as it can not be expressed as a
   *        linear combination of gaussian basis functions very well. Additionally,
   *        the Fermi-Amaldi potential is a popular choice.
   */
  void calculateOEP(const DensityOnGrid<SCFMode>& targetDens, std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                    const GridPotential<SCFMode>& initialGuess);

  /**
   * @brief Performs a van Leeuwen--Baerends potential reconstruction
   *
   *        For more information on the theory behind this see Ref.:
   *        van Leeuwen, R.; Baerends, E. J.; Phys. Rev. A 49, 2421 (1994)
   * @param targetDens The target density to be reconstructed.
   * @param resultOrbitals The orbitals of the reconstructed system.
   * @param initialGuess An initial guess for the new potential.
   * @param damping Damping for the density update in the van Leeuwen–Barends potential reconstruction.
   * @param maxCycles Maximum of iterative cycles to be performed.
   */
  SpinPolarizedData<SCFMode, double>
  calculateOEPLB(const DensityOnGrid<SCFMode>& targetDens, std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                 const FockMatrix<SCFMode>& initialGuess, double damping = 0.95, unsigned int maxCycles = 5000);
  /**
   * @brief Overload to take a potential expressed on a grid.
   *        Converts the GridPotential into a FockMatrix in
   *        the correct basis and calls the function above.
   */
  SpinPolarizedData<SCFMode, double>
  calculateOEPLB(const DensityOnGrid<SCFMode>& targetDens, std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                 const GridPotential<SCFMode>& initialGuess, double damping = 0.95, unsigned int maxCycles = 5000);

  /**
   * @brief Performs a Zhang-Carter potential reconstruction.
   *        For more information on the theory behind this see Ref.:
   *        X. Zhang and E. A. Carter, J. Chem. Phys. 148, 034105 (2018)
   * @param targetDens The target density to be reconstructed.
   * @param resultOrbitals The orbitals of the reconstructed system.
   * @param initialGuess An initial guess for the new potential.
   * @param maxCycles Maximum of iterative cycles to be performed.
   */
  void calculateOEPCarter(const DensityMatrix<SCFMode>& targetDens, std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                          const FockMatrix<SCFMode>& initialGuess, const unsigned int maxCycles = 5000);

  /**
   * @brief Calculates the gradient of the density with respect to the basis functions used in the Wu-Yang scheme
   * @param targetDens The target density to be reconstructed.
   * @param densMatOnGrid
   * @param OEPCoeffs Coefficients of the potential expressed in basis functions (Wu-Yang)
   * @return
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> getGradient(const DensityOnGrid<SCFMode>& targetDens,
                                                          const DensityOnGrid<SCFMode>& densMatOnGrid,
                                                          SpinPolarizedData<SCFMode, Eigen::VectorXd>& OEPCoeffs);

 private:
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>>
  getHessian(std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals);
  void updateDensity(DensityOnGrid<SCFMode>& densMatOnGrid, std::shared_ptr<OrbitalController<SCFMode>> resultOrbitals,
                     const FockMatrix<SCFMode>& initialGuess, const SpinPolarizedData<SCFMode, Eigen::VectorXd>& OEPCoeffs);

  FockMatrix<SCFMode> getX(DensityMatrix<SCFMode>& densityMatrix);

  std::shared_ptr<BasisFunctionOnGridController> _basisFuncOnGridController;
  std::shared_ptr<BasisFunctionOnGridController> _potBasisFuncOnGridController;
  std::shared_ptr<OneElectronIntegralController> _oneEIntController;
  std::shared_ptr<DensityOnGridCalculator<SCFMode>> _densOnGridCalculator;
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  std::unique_ptr<GridPotential<SCFMode>> _potentialOnGrid;
  const SpinPolarizedData<SCFMode, unsigned int>& _nOccOrbs;
  const double _smoothFactor;
  const double _singValThreshold;
  const double _exc;
};

} /* namespace Serenity */

#endif /* POTENTIALS_OPTEFFPOTENTIAL_H_ */
