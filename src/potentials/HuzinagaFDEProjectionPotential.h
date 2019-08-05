/**
 * @file HuzinagaFDEProjectionPotential.h
 *
 * @date Nov 23, 2017
 * @author Moritz Bensberg
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

#ifndef POTENTIALS_HUZINAGAFDEPROJECTIONPOTENTIAL_H_
#define POTENTIALS_HUZINAGAFDEPROJECTIONPOTENTIAL_H_

/* Include Serenity Internal Headers */
#include "data/grid/ExternalDensityOnGridController.h"
#include "data/ElectronicStructure.h"
#include "data/matrices/DensityMatrixController.h"
#include "dft/Functional.h"
#include "grid/GridController.h"
#include "potentials/bundles/FDEPotentials.h"
#include "potentials/Potential.h"
#include "settings/Options.h"
#include "potentials/ABFockMatrixConstruction/ABPotential.h"

namespace Serenity {

/* Forward declarations */
class SystemController;
/**
 * @class HuzinagaFDEProjectionPotential HuzinagaFDEProjectionPotential.h
 * @brief A class for handling the Huzinaga operator and Hoffmann's external orthogonality
 *        approach in any FDE like calculation.
 *
 * The Huzinaga operator:\n
 * According to:\n
 *  J. Chem. Theory Comput. 13, 1503-1508 (2017)\n
 *  J. Chem. Phys. 55, 5543 (1971)
 *
 * Hoffmann's external orthogonality approach:\n
 * According to:\n
 *  Theory. Annu. Rep. Comput. Chem. 8,53-70 (2012) and \n
 *  J. Phys. Chem. A 118, 9182–9200 (2014).
 */
template <Options::SCF_MODES SCFMode>
class HuzinagaFDEProjectionPotential: public Potential<SCFMode>,
public ObjectSensitiveClass<Basis>,
public ObjectSensitiveClass<DensityMatrix<SCFMode> >  {
public:
  /**
   * @brief Constructor.
   *
   * This constructor prepares the calculations by building pairs of environment systems with the active system.
   * These pairs are needed in order to calculate the Fock matrix and overlap matrix of the joint system.
   * The constructor expects a guess electronic structure of the active system with SPIN=SCFMode!
   *
   * @param activeSystem The active system.
   * @param environmentSystems The environment systems.
   * @param truncatedProjection A flag in order to truncate the projection operator.
   * @param projectionOverlapThreshold The overlap threshold for the truncation.
   * @param buildHoffmannsOperator A flag whether Hoffmann's operator should be constructed.
   * @param distantKinFunc A flag whether not projected subsystems are treated with a kin. energy func. .(optional)
   * @param supersystemgrid The supersystem grid for the kin. energy func. evaluation. (optional)
   * @param naddKinFunc The non--additive kin. energy functional. (optional)
   * @param gridCutOff The grid cut off. (optional)
   * @param allEConts The energy component controllers. (optional)
   */
  HuzinagaFDEProjectionPotential(
      std::shared_ptr<SystemController> activeSystem,
      std::vector<std::shared_ptr<SystemController> > environmentSystems,
      Functional naddXCFunc,
      bool truncatedProjection = false,
      double projectionOverlapThreshold = 0.0,
      bool buildHoffmannsOperator = false,
      bool topDown = false,
      bool distantKinFunc = false,
      double basisFunctionRatio = 0.0,
      double borderAtomThreshold = 0.02,
      std::shared_ptr<GridController> supersystemgrid = nullptr,
      Options::KINFUNCTIONALS naddKinFunc = Options::KINFUNCTIONALS::NONE,
      double gridCutOff = 0.0,
      std::vector<std::shared_ptr<EnergyComponentController> > allEConts = {nullptr}
      );
  /**
   * @brief Default destructor.
   */
  virtual ~HuzinagaFDEProjectionPotential()=default;

  /**
   * @return The fock matrix for the embedded/active system.
   */
  FockMatrix<SCFMode>& getMatrix() override final;

  /**
   * @return Dummy(0) matrix.
   */
  Eigen::MatrixXd getGeomGradients() override final;

  /**
   * @param P The density matrix.
   * @return 0. There is no energy associated with this potential.
   *
   *TODO
   * In general the orbitals may not be orthogonal after applying this projection operator due to
   * the lack of basis functions in the space of the projected orbitals.
   * This would result in the so called overlap energy. For more information see:
   *   J. Chem. Phys. 97, 9 (1992) 6504-6508
   */
  double getEnergy(const DensityMatrix<SCFMode>& P) override final;
  /**
   * @brief Potential is linked to the grid and density.
   *        This is used for lazy evaluation.
   *        (see ObjectSensitiveClass and NotifyingClass)
   */
  void notify() override final {
    _potential = nullptr;
  }

private:
  /// @brief The potential in matrix representation.
  std::unique_ptr<FockMatrix<SCFMode> >_potential;
  /// @brief The active system.
  std::shared_ptr<SystemController> _activeSystem;
  /// @brief The environment systems.
  std::vector<std::shared_ptr<SystemController> > _environmentSystems;
  /// @brief The density matrix controllers of the environment systems.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > _envDensityCont;
  /// @brief the supersystem geometry (without duplicated atoms!)
  std::shared_ptr<Geometry> _supersystemGeometry;
  /// @brief The transformer for the exchange--correlation part for all AB grid combinations.
  std::vector<std::shared_ptr<ScalarOperatorToMatrixAdder<SCFMode> > > _gridToMatrix_AB;
  ///@brief The non-additive exchange--correlation functional.
  Functional _naddXCFunc;
  /*
   * Projection truncation
   */
  /// @brief The overlap matrix of the active system basis set with all environment basis sets.
  std::vector<std::shared_ptr<Eigen::MatrixXd> >_s_ABs;
  /// @brief A flag whether the projection operator should be truncated.
  bool _truncatedProjection;
  /// @brief A flag for each environment system, which is true when the system is projected into the
  /// active system
  Eigen::VectorXi _projectedEnvSystems;
  /// @brief The overlap threshold for the determination of the projected environment systems.
  double _projectionOverlapThreshold;
  /*
   * Hoffmann only
   */
  /// @brief A flag whether Hoffmann's projection operator should be constructed.
  bool _buildHoffmannsOperator;
  /// @brief The FDE potentials for the environment systems (f_BB) (Hoffmann only).
  std::vector<std::shared_ptr<FDEPotentials<SCFMode> > > _env_i_fdePot;

  /// @brief Special case of a top-down calculation.
  bool _topDown;
  /*
   * Distant kin. functional only
   */
  /// @brief A flag whether not projected subsystems should be treated with a non additve kin. energy functional.
  bool _distantKinFunc;
  /// @brief The supersystem grid (optional).
  std::shared_ptr<GridController> _supersystemgrid;
  /// @brief The non additive kinetic energy potential for not projected subsystems (optional).
  std::shared_ptr<NAddFuncPotential<SCFMode> > _naddKinPot;

  unsigned int _printLevel;
  /// @brief The density controllers of the not-projected densities.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode> > > _notProjectedEnvDensities;
  /*
   * AB-Potentials
   */
  /// @brief The coulomb contributions in the basis A x B for all environment systems.
  std::vector<std::shared_ptr<ABPotential<SCFMode> > > _coulombAB_Potentials;
  /// @brief The coulomb contribution in the basis A x B for the active systems.
  std::vector<std::shared_ptr<ABPotential<SCFMode> > > _activeCoulombAB_Potentials;
  /// @brief The exchange--correlation contribution in the basis A x B for the active system on all AB grids.
  std::vector<std::shared_ptr<ABPotential<SCFMode> > > _exchangeAB_Potentials;
  /// @brief The core hamiltonian in the basis A x B.
  std::vector<std::shared_ptr<ABPotential<SCFMode> > > _coreAB_Potentials;
  /// @brief The non-additive exchange--correlation potential for all environment systems.
  std::vector<std::shared_ptr<ABPotential<SCFMode> > > _naddExchangeAB_Potentials;
  /// @brief The non-additive exchange--correlation potential for all environment systems.
  std::vector<std::shared_ptr<ABPotential<SCFMode> > > _naddKinAB_Potentials; /// only for hybrid methods
  /// @brief The combined AB auxillary basis.
  std::vector<std::shared_ptr<BasisController> > _abAuxBasisController;

  /* Helper Functions */
  /**
   * @brief builds the outer diagonal Fock matrix of act. + env[iEnv]
   * @param iEnv The index of the environment system.
   */
  SPMatrix<SCFMode> buildOuterDiagonalFockMatrix(
      unsigned int iEnv);
  /**
   * @brief Writes the current average overlap of all occupied environmental orbitals with the occ. active system orbitals.
   */
  void writeInterSubsystemOccOverlap();
};

} /* namespace Serenity */

#endif /* POTENTIALS_HUZINAGAFDEPROJECTIONPOTENTIAL_H_ */
