/**
 * @file TDReconstructionPotential.h
 *
 * @date Apr 03, 2018
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

#ifndef POTENTIALS_TDRECONSTRUCTIONPOTENTIAL_H_
#define POTENTIALS_TDRECONSTRUCTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "potentials/Potential.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <memory>
#include <vector>
namespace Serenity {

// Todo add tests!

class SystemController;

/**
 * @class TDReconstructionPotential TDReconstructionPotential.h
 *
 * @brief A class to perform a top-down potential reconstruction.
 *        Details concerning the procedure/implementation: J. Chem. Phys. 149, 054103 (2018)
 */
template<Options::SCF_MODES SCFMode>
class TDReconstructionPotential : public Potential<SCFMode> {
 public:
  /**
   * @brief Constructor.
   * @param actSys The active system.
   * @param supSys The supersystem.
   * @param envSystems The environment systems.
   * @param smoothFactor The smooth factor.
   * @param potBasisLabel The label of the basis used for the Wu-Yang potential reconstruction.
   * @param singValThreshold Threshold for the singular value decomposition.
   * @param lbDamping The damping factor for the van Leeuwen-Baerends reconstruction.
   * @param lbCycles Number of cycles for the van Leeuwen-Baerends reconstruction.
   * @param carterCycles Number of cycles for the Huang--Carter reconstruction.
   * @param noSupRec Do not reconstruct the supersystem potential.
   */
  TDReconstructionPotential(std::shared_ptr<SystemController> actSys, std::shared_ptr<SystemController> supSys,
                            std::vector<std::shared_ptr<SystemController>> envSystems, double smoothFactor = 0.0,
                            std::string potBasisLabel = "", const double singValThreshold = 0.0, double lbDamping = 0.995,
                            unsigned int lbCycles = 0, unsigned int carterCycles = 0, bool noSupRec = false);

  /// @brief Default destructor.
  virtual ~TDReconstructionPotential() = default;

  /**
   * @brief Getter for the potential matrix.
   * @return The potential in its matrix representation.
   */
  FockMatrix<SCFMode>& getMatrix() override {
    if (!this->_potential) {
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
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;

  /**
   * @brief Geometry gradient contribution from this Potential.
   *
   * @return The geometry gradient contribution resulting from this Potential.
   */
  Eigen::MatrixXd getGeomGradients() override;

 private:
  void calculatePotential();

  ///@brief The potential.
  std::unique_ptr<FockMatrix<SCFMode>> _potential;
  std::weak_ptr<SystemController> _actSys;
  std::weak_ptr<SystemController> _supSys;
  std::vector<std::weak_ptr<SystemController>> _envSystems;
  double _smoothFactor;
  std::string _potBasisLabel;
  const double _singValThreshold;
  double _lbDamping;
  unsigned int _lbCycles;
  unsigned int _carterCycles;
  bool _noSupRec;
  double _supEnergy;
  double _supXEnergy;
};

} /* namespace Serenity */

#endif /* POTENTIALS_TDRECONSTRUCTIONPOTENTIAL_H_ */
