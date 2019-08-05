/**
 * @file BUReconstructionPotential.h
 *
 * @date Oct 13, 2016
 * @author David Schnieders
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

#ifndef POTENTIALS_BURECONSTRUCTIONPOTENTIAL_H_
#define POTENTIALS_BURECONSTRUCTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "data/grid/GridPotential.h"
#include "potentials/NAddFuncPotential.h"
#include "notification/ObjectSensitiveClass.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>
namespace Serenity {

class SystemController;
template<Options::SCF_MODES SCFMode>class DensityMatrixController;

template<Options::SCF_MODES SCFMode>
class BUReconstructionPotential : public NAddFuncPotential<SCFMode>{
public:
  /**
   * @brief Constructor
   * @param actSys The active system.
   * @param supSys The supersystem.
   */
  BUReconstructionPotential(std::shared_ptr<SystemController> actSys,
      std::shared_ptr<DensityMatrixController<SCFMode> > activeDMat,
      std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > envDMats,
      std::shared_ptr<GridController> superSystemGrid,
      Functional functional,
      std::vector<std::shared_ptr<SystemController>> envSystems,
      double smoothFactor=0.0,
      std::string potBasisLabel="",
      const double singValThreshold=0.0,
      double lbDamping = 0.995,
      unsigned int lbCycles=0,
      unsigned int carterCycles=0);
  virtual ~BUReconstructionPotential() = default;

  /**
   * @brief Getter for the potential matrix.
   * @return The potential in its matrix representation.
   */
  FockMatrix<SCFMode>& getMatrix() override final{
    if(!this->_potential){
      calculatePotential();
    };
    return *this->_potential;
  };

  /**
   * @brief Getter for the energy associated with this potential.
   * @param P The density matrix for which the energy should be
   * calculated.
   * @return The energy associated with the potential and P.
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final{
    if(!this->_potential){
      calculatePotential();
    };
    double energy=_supEnergy;
    auto& libint = Libint::getInstance();
    auto kin=libint.compute1eInts(libint2::Operator::kinetic, this->_system->getBasisController());
    for_spin(P){
      double kinE=P_spin.cwiseProduct(kin).sum();
      energy-=kinE;
    };

    for(auto sys : this->_envSystems){
      auto kin=libint.compute1eInts(libint2::Operator::kinetic, sys->getBasisController());
      auto mat=sys->template getElectronicStructure<SCFMode>()->getDensityMatrix();
      for_spin(mat){
        double kinE=mat_spin.cwiseProduct(kin).sum();
        energy-=kinE;
      };
    }

    return energy;
  };

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */

  Eigen::MatrixXd getGeomGradients() override;

  /**
   * @brief Overrides the notifying function of NAddFuncPotential
   *        to do nothing (we do not want to reperform the potential
   *        reconstructions on every SCF step).
   */
  void notify() override final{
  };

private:

  /**
   * @brief Calculates the non-additive kinetic potential as described
   * in J. Chem. Phys. 133, 084103 (2010)
   */
  void calculatePotential();

  std::shared_ptr<BasisController> _basis;
  std::vector<std::shared_ptr<SystemController>> _envSystems;
  double _smoothFactor;
  std::string _potBasisLabel;
  const double _singValThreshold;
  double _lbDamping;
  unsigned int _lbCycles;
  unsigned int _carterCycles;
  std::unique_ptr<GridPotential<SCFMode>> _potentialOnGrid;
  double _supEnergy;


};

} /* namespace Serenity */

#endif /* BASICS_GRID_DATAONGRID_EXACTNADDKINPOTCALCULATOR_H_ */
