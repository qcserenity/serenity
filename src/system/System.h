/**
 * @file   System.h
 *
 * @date   Mar 25, 2013
 * @author Thomas Dresselhaus
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
#ifndef SYSTEM_H_
#define SYSTEM_H_
/* Include Serenity Internal Headers */
#include "settings/BasisOptions.h"
#include "settings/GridOptions.h"
#include "settings/Settings.h"
/* Include Std and External Headers */
#include <cassert>
#include <map>
#include <memory>
#include <string>

namespace Serenity {
/* Forward declarations */
class AtomCenteredBasisController;
class BasisController;
template<Options::SCF_MODES SCFMode>
class ElectronicStructure;
class Geometry;
class GridController;
class OneElectronIntegralController;
template<Options::SCF_MODES SCFMode>
class OrbitalController;
class SystemController;
enum class MOLECULAR_SURFACE_TYPES;
template<Options::SCF_MODES SCFMode>
class ElectrostaticPotentialOnGridController;
class MolecularSurfaceController;
class ExternalChargeController;
/**
 * @class System System.h
 * @brief Some kind of system on which can be worked, e.g. a molecule
 *
 * or groups of molecules, or a single atom... It is defined by a geometry, charge
 * and spin, and possibly also a name. Additional information like basis set(s),
 * integrals, orbitals and so on can be stored alongside with it. The SystemController
 * manages objects of this class.
 */
class System {
  friend SystemController;

 public:
  /**
   * @param charge
   * @param spin
   * @param geometry
   * @param name
   */
  System(std::shared_ptr<Geometry> geometry, Settings settings)
    : _geometry(geometry),
      _settings(settings),
      _lastSCFMode(_settings.spin == 0 ? Options::SCF_MODES::RESTRICTED : Options::SCF_MODES::UNRESTRICTED) {
    assert(_geometry);
  };
  virtual ~System() = default;

 private:
  std::shared_ptr<Geometry> _geometry;
  std::unique_ptr<Eigen::MatrixXd> _pointChargeGradients;
  std::shared_ptr<ExternalChargeController> _externalChargeController;
  Settings _settings;
  std::map<Options::BASIS_PURPOSES, std::shared_ptr<AtomCenteredBasisController>> _basisControllers;
  std::map<Options::GRID_PURPOSES, std::shared_ptr<GridController>> _gridControllers;

  std::map<MOLECULAR_SURFACE_TYPES, std::shared_ptr<MolecularSurfaceController>> _molecularSurfaces;
  std::map<MOLECULAR_SURFACE_TYPES, std::shared_ptr<ElectrostaticPotentialOnGridController<Options::SCF_MODES::RESTRICTED>>> _restrictedElectrostaticPotentialsOnGridControllers;
  std::map<MOLECULAR_SURFACE_TYPES, std::shared_ptr<ElectrostaticPotentialOnGridController<Options::SCF_MODES::UNRESTRICTED>>>
      _unrestrictedElectrostaticPotentialsOnGridControllers;

  enum Options::SCF_MODES _lastSCFMode;

  std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>> _restrictedElectronicStructure;
  std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> _unrestrictedElectronicStructure;
  unsigned int _nElectrons;
};

} /* namespace Serenity */

#endif /* SYSTEM_H_ */
