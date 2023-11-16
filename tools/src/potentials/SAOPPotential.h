/**
 * @file   SAOPPotential.h
 *
 * @date   Aug 4, 2017
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

#ifndef POTENTIALS_SAOPPOTENTIAL_H_
#define POTENTIALS_SAOPPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "grid/Grid.h"
#include "notification/ObjectSensitiveClass.h"
#include "potentials/Potential.h"

namespace Serenity {

/* Forward Declarations */
class SystemController;
template<Options::SCF_MODES SCFMode>
class ScalarOperatorToMatrixAdder;
template<Options::SCF_MODES SCFMode>
class DensityOnGridController;
class BasisController;
template<Options::SCF_MODES SCFMode>
class OrbitalController;
class BasisFunctionOnGridController;

/**
 * @class SAOPPotential SAOPPotential.h
 *
 * This potential was created to approximate the exchange-correlation Kohn-Sham potential
 * with a Statistical Average of different Orbital model Potentials (SAOP). The model
 * potential for electrons close to the nucleus is the GLB potential \f$ \nu^{GLB}_{xc~\sigma}(r) \f$ .
 * Electrons in higher orbitals are described by the LB potential \f$ \nu^{LB}_{xc~\sigma}(r) \f$.
 *
 * The SAOP potential is defined as:
 * \f$ \nu_{xc \sigma}^\mathrm{SAOP}(r)=
 * \sum_{i=1}^{N_\sigma} \nu_{xc~i~\sigma}^\mathrm{mod}(r) \frac{|\Psi_{i~\sigma}(r)|^2}{\rho_\sigma (r)} \f$
 *
 * with \f$ \nu_{xc~i~\sigma}^\mathrm{mod}(r)=exp[-2(\epsilon_{N_\sigma}-\epsilon_{i~\sigma})^2]
 * \nu_{xc~\sigma}^\mathrm{LB}
 *                                           + \{
 * 1-exp[-2(\epsilon_{N_\sigma}-\epsilon_{i~\sigma})^2]\}\nu^{GLB}_{xc~\sigma}(r)\f$,\n
 *
 * \f$ \nu^{GLB}_{xc~\sigma}(r) = 2(\epsilon_{x\sigma}^\mathrm{B88}+\epsilon_{c}^\mathrm{PW91})+\sum_{i=1}^{N_\sigma}
 *  \sqrt{\epsilon_{N_\sigma}-\epsilon_{i~\sigma}}\frac{|\Psi_{i~\sigma}(r)|^2}{\rho_\sigma (r)} \f$,\n
 * \f$ \nu^{LB}_{xc~\sigma}(r) = \alpha \nu_{x\sigma}^\mathrm{LDA}(r) + \nu_{c\sigma}^\mathrm{LDA}(r)
 * -\frac{\beta\chi_\sigma^2(r)\rho_\sigma^{1/3}(r)}{1+3\beta\chi_\sigma(r)\ln\{\chi_\sigma(r)+\sqrt{\chi_\sigma^2(r)+1}\}}
 * \f$.\n \n For further information see:\n
 *  /*!<  P. R. T. Schipper,O.V. Gritsenko, S.J. A. van Gisbergen and E.J. Baerends, J. Chem. Phys. 112, 1344 (2000)
 */

template<Options::SCF_MODES SCFMode>
class SAOPPotential : public Potential<SCFMode>,
                      public ObjectSensitiveClass<Grid>,
                      public ObjectSensitiveClass<DensityMatrix<SCFMode>> {
 public:
  /**
   * @brief Constructor for the SAOP potential
   * @param basis The basis controller.
   * @param systemController The system controller. Needed for the settings.
   * @param densityOnGridController The controller for the density in grid representation.
   * @param grid The controller for the used grid.
   */
  SAOPPotential(std::shared_ptr<BasisController> basis, std::shared_ptr<SystemController> systemController,
                std::shared_ptr<DensityOnGridController<SCFMode>> densityOnGridController);
  ///@brief Default destructor.
  virtual ~SAOPPotential() = default;
  /**
   * @brief Getter for the potential in matrix representation.
   * @return The potential as a matrix.
   */
  FockMatrix<SCFMode>& getMatrix() override final;
  /**
   * @brief Getter for the energy associated with the potential.
   * @param P The density associated with the energy. In this case it is only a dummy!
   * @return The energy associated with the potential. Here it equals the sum of
   *         the Becke 88 exchange and Perdew-Wang 91 correlation energy.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;
  /**
   * @brief Dummy-Getter for the gradient.
   * @return Returns 0. No gradients available for this potential!
   */
  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @brief Potential is linked to the grid and density.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  };

 private:
  ///@brief The controller for the total density in grid representation.
  std::shared_ptr<DensityOnGridController<SCFMode>> _densityOnGridController;
  ///@brief Number of occupied orbitals.
  SpinPolarizedData<SCFMode, unsigned int> _nOCC;
  ///@brief The controller for the current orbital set.
  std::shared_ptr<OrbitalController<SCFMode>> _orbitalController;

  std::shared_ptr<SystemController> _systemController;
  ///@brief The controller for the basis function values the grid.
  std::shared_ptr<BasisFunctionOnGridController> _basisFunctionOnGridController;
  ///@brief occupation number of the orbitals.
  unsigned int _occupations;
  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  ///@brief The energy.
  double _energy;
  ///@brief The parameter alpha of the LB potential.
  double _alpha = 1.19;
  ///@brief The parameter beta of the LB potential.
  double _beta = 0.01;
  ///@brief Conversion from grid to matrix.
  std::shared_ptr<ScalarOperatorToMatrixAdder<SCFMode>> _gridToMatrix;
  ///@brief The grid this potential is defined on.
  std::shared_ptr<GridController> _grid;
  ///@brief The parameter K used in the GLB potential.
  double _K = 0.42;
};

} /* namespace Serenity */

#endif /* POTENTIALS_SAOPPOTENTIAL_H_ */
