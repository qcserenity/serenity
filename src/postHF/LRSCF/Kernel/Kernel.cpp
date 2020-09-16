/**
 * @file   Kernel.cpp
 *
 * @date   Mar 19, 2017
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

/* Include Class Header*/
#include "postHF/LRSCF/Kernel/Kernel.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "data/grid/BasisFunctionOnGridControllerFactory.h"
#include "data/grid/DensityMatrixDensityOnGridController.h"
#include "data/grid/DensityOnGrid.h"
#include "data/grid/DensityOnGridCalculator.h"
#include "data/grid/ExternalDensityOnGridController.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/DensityMatrixController.h"
#include "dft/functionals/wrappers/XCFun.h"
#include "geometry/Geometry.h"
#include "grid/GridControllerFactory.h"
#include "misc/SerenityError.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"

namespace Serenity {

template<Options::SCF_MODES T>
Kernel<T>::Kernel(std::vector<std::shared_ptr<SystemController>> act, std::vector<std::shared_ptr<SystemController>> env,
                  bool superSystemGrid, CompositeFunctionals::KINFUNCTIONALS naddKinFunc,
                  CompositeFunctionals::XCFUNCTIONALS naddXCFunc, CompositeFunctionals::XCFUNCTIONALS func)
  : _superSystemGrid(superSystemGrid),
    _naddKinFunc(naddKinFunc),
    _naddXCFunc(naddXCFunc),
    _func(func),
    _gridController(act[0]->getGridController()),
    _gga(false) {
  // Build systems
  for (auto sys : act) {
    _systems.push_back(sys);
  }
  for (auto sys : env) {
    _systems.push_back(sys);
  }
  // Build supersystem grid if desired
  if (_superSystemGrid) {
    // Build supersystem geometry
    std::shared_ptr<Geometry> geometry = std::make_shared<Geometry>();
    for (unsigned int I = 0; I < _systems.size(); ++I) {
      (*geometry) += (*_systems[I]->getGeometry());
    }
    geometry->deleteIdenticalAtoms();
    // produce grid controller for supersystem geometry
    _gridController = GridControllerFactory::produce(geometry, _systems[0]->getSettings(), Options::GRID_PURPOSES::DEFAULT);
  }
  assert(_gridController);
  // Get and check data from input functionals
  if (_systems.size() == 1) {
    // If only one system is present, non-additive contribution is zero and need not be calculated.
    // Though NONE is default, we reassure that non-additive functionals are not used.
    _naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
    _naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::NONE;
  }
  Functional EnaddXC = resolveFunctional(naddXCFunc);
  Functional EnaddKIN = resolveFunctional(naddKinFunc);
  if (EnaddKIN.isHybrid() || EnaddKIN.isRSHybrid())
    throw SerenityError("ERROR: Found hybrid functional for non-additive kinetic kernel.");
  if (EnaddXC.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA ||
      EnaddKIN.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA)
    _gga = true;
  for (auto sys : _systems) {
    Functional EXC = resolveFunctional(sys->getSettings().dft.functional);
    if (EXC.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA)
      _gga = true;
  }
  if (resolveFunctional(_func).getFunctionalClass() == CompositeFunctionals::CLASSES::GGA)
    _gga = true;
  // Initialize storage objects (this ensures that key is already present when working with it)
  for (auto sys : _systems) {
    _pp.insert(std::make_pair(sys->getSettings().name,
                              DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>(_gridController)));
    if (_gga)
      _pg.insert(std::make_pair(
          sys->getSettings().name,
          makeGradient<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController)));
    if (_gga)
      _gg.insert(std::make_pair(
          sys->getSettings().name,
          makeHessian<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController)));
  }
  if (_systems.size() > 1) {
    _pptot = std::make_shared<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
    _ggtot = makeHessianPtr<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
    _pgtot = makeGradientPtr<DoublySpinPolarizedData<T, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
  }
  // Calculate and store derivatives
  calculateDerivatives();
}

template<>
const std::shared_ptr<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>
Kernel<Options::SCF_MODES::RESTRICTED>::getPP(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart) {
  auto pp_Block = std::make_shared<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>(blockSize);
  (*pp_Block).setZero();
  // If more than one system for the non-additive contributions
  if (_systems.size() > 1)
    *pp_Block += (*_pptot).segment(iGridStart, blockSize);
  // For non-additive contributions in case of only more than one subsystem
  // Or evaluation of the entire Kernel for only one subsystem (Supersystem)
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& pp_system = _pp.find(name)->second;
    *pp_Block += pp_system.segment(iGridStart, blockSize);
  }
  return pp_Block;
}

template<>
const std::shared_ptr<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>
Kernel<Options::SCF_MODES::UNRESTRICTED>::getPP(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart) {
  auto pp_Block = std::make_shared<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>(blockSize);
  (*pp_Block).aa.setZero();
  (*pp_Block).ba.setZero();
  (*pp_Block).ab.setZero();
  (*pp_Block).bb.setZero();
  if (_systems.size() > 1) {
    (*pp_Block).aa += (*_pptot).aa.segment(iGridStart, blockSize);
    (*pp_Block).ab += (*_pptot).ab.segment(iGridStart, blockSize);
    ;
    (*pp_Block).ba += (*_pptot).ba.segment(iGridStart, blockSize);
    ;
    (*pp_Block).bb += (*_pptot).bb.segment(iGridStart, blockSize);
    ;
  }
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& ppsub = _pp.find(name)->second;
    (*pp_Block).aa += ppsub.aa.segment(iGridStart, blockSize);
    (*pp_Block).ab += ppsub.ab.segment(iGridStart, blockSize);
    (*pp_Block).ba += ppsub.ba.segment(iGridStart, blockSize);
    (*pp_Block).bb += ppsub.bb.segment(iGridStart, blockSize);
  }
  return pp_Block;
}

template<>
const std::shared_ptr<Gradient<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>>
Kernel<Options::SCF_MODES::RESTRICTED>::getPG(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart) {
  auto pg_Block = std::make_shared<Gradient<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>>(
      makeGradient<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>(blockSize));
  (*pg_Block).x.setZero();
  (*pg_Block).y.setZero();
  (*pg_Block).z.setZero();
  if (_systems.size() > 1) {
    (*pg_Block).x += (*_pgtot).x.segment(iGridStart, blockSize);
    (*pg_Block).y += (*_pgtot).y.segment(iGridStart, blockSize);
    (*pg_Block).z += (*_pgtot).z.segment(iGridStart, blockSize);
  }
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& pg_system = _pg.find(name)->second;
    (*pg_Block).x += pg_system.x.segment(iGridStart, blockSize);
    (*pg_Block).y += pg_system.y.segment(iGridStart, blockSize);
    (*pg_Block).z += pg_system.z.segment(iGridStart, blockSize);
  }
  return pg_Block;
}

template<>
const std::shared_ptr<Gradient<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>>
Kernel<Options::SCF_MODES::UNRESTRICTED>::getPG(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart) {
  auto pg_Block = std::make_shared<Gradient<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>>(
      makeGradient<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>(blockSize));

  (*pg_Block).x.aa.setZero();
  (*pg_Block).x.ab.setZero();
  (*pg_Block).x.ba.setZero();
  (*pg_Block).x.bb.setZero();
  (*pg_Block).y.aa.setZero();
  (*pg_Block).y.ab.setZero();
  (*pg_Block).y.ba.setZero();
  (*pg_Block).y.bb.setZero();
  (*pg_Block).z.aa.setZero();
  (*pg_Block).z.ab.setZero();
  (*pg_Block).z.ba.setZero();
  (*pg_Block).z.bb.setZero();

  if (_systems.size() > 1) {
    (*pg_Block).x.aa += (*_pgtot).x.aa.segment(iGridStart, blockSize);
    (*pg_Block).x.ab += (*_pgtot).x.ab.segment(iGridStart, blockSize);
    (*pg_Block).x.ba += (*_pgtot).x.ba.segment(iGridStart, blockSize);
    (*pg_Block).x.bb += (*_pgtot).x.bb.segment(iGridStart, blockSize);
    (*pg_Block).y.aa += (*_pgtot).y.aa.segment(iGridStart, blockSize);
    (*pg_Block).y.ab += (*_pgtot).y.ab.segment(iGridStart, blockSize);
    (*pg_Block).y.ba += (*_pgtot).y.ba.segment(iGridStart, blockSize);
    (*pg_Block).y.bb += (*_pgtot).y.bb.segment(iGridStart, blockSize);
    (*pg_Block).z.aa += (*_pgtot).z.aa.segment(iGridStart, blockSize);
    (*pg_Block).z.ab += (*_pgtot).z.ab.segment(iGridStart, blockSize);
    (*pg_Block).z.ba += (*_pgtot).z.ba.segment(iGridStart, blockSize);
    (*pg_Block).z.bb += (*_pgtot).z.bb.segment(iGridStart, blockSize);
  }
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& pgsub = _pg.find(name)->second;
    (*pg_Block).x.aa += pgsub.x.aa.segment(iGridStart, blockSize);
    (*pg_Block).x.ab += pgsub.x.ab.segment(iGridStart, blockSize);
    (*pg_Block).x.ba += pgsub.x.ba.segment(iGridStart, blockSize);
    (*pg_Block).x.bb += pgsub.x.bb.segment(iGridStart, blockSize);
    (*pg_Block).y.aa += pgsub.y.aa.segment(iGridStart, blockSize);
    (*pg_Block).y.ab += pgsub.y.ab.segment(iGridStart, blockSize);
    (*pg_Block).y.ba += pgsub.y.ba.segment(iGridStart, blockSize);
    (*pg_Block).y.bb += pgsub.y.bb.segment(iGridStart, blockSize);
    (*pg_Block).z.aa += pgsub.z.aa.segment(iGridStart, blockSize);
    (*pg_Block).z.ab += pgsub.z.ab.segment(iGridStart, blockSize);
    (*pg_Block).z.ba += pgsub.z.ba.segment(iGridStart, blockSize);
    (*pg_Block).z.bb += pgsub.z.bb.segment(iGridStart, blockSize);
  }
  return pg_Block;
}

template<>
const std::shared_ptr<Hessian<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>>
Kernel<Options::SCF_MODES::RESTRICTED>::getGG(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart) {
  auto gg_Block = std::make_shared<Hessian<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>>(
      makeHessian<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>(blockSize));

  (*gg_Block).xx.setZero();
  (*gg_Block).xy.setZero();
  (*gg_Block).xz.setZero();
  (*gg_Block).yy.setZero();
  (*gg_Block).yz.setZero();
  (*gg_Block).zz.setZero();

  if (_systems.size() > 1) {
    (*gg_Block).xx += (*_ggtot).xx.segment(iGridStart, blockSize);
    (*gg_Block).xy += (*_ggtot).xy.segment(iGridStart, blockSize);
    (*gg_Block).xz += (*_ggtot).xz.segment(iGridStart, blockSize);
    (*gg_Block).yy += (*_ggtot).yy.segment(iGridStart, blockSize);
    (*gg_Block).yz += (*_ggtot).yz.segment(iGridStart, blockSize);
    (*gg_Block).zz += (*_ggtot).zz.segment(iGridStart, blockSize);
  }
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& gg_system = _gg.find(name)->second;
    (*gg_Block).xx += gg_system.xx.segment(iGridStart, blockSize);
    (*gg_Block).xy += gg_system.xy.segment(iGridStart, blockSize);
    (*gg_Block).xz += gg_system.xz.segment(iGridStart, blockSize);
    (*gg_Block).yy += gg_system.yy.segment(iGridStart, blockSize);
    (*gg_Block).yz += gg_system.yz.segment(iGridStart, blockSize);
    (*gg_Block).zz += gg_system.zz.segment(iGridStart, blockSize);
  }
  return gg_Block;
}

template<>
const std::shared_ptr<Hessian<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>>
Kernel<Options::SCF_MODES::UNRESTRICTED>::getGG(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart) {
  auto gg_Block = std::make_shared<Hessian<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>>(
      makeHessian<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>(blockSize));

  (*gg_Block).xx.aa.setZero();
  (*gg_Block).xx.ab.setZero();
  (*gg_Block).xx.ba.setZero();
  (*gg_Block).xx.bb.setZero();
  (*gg_Block).xy.aa.setZero();
  (*gg_Block).xy.ab.setZero();
  (*gg_Block).xy.ba.setZero();
  (*gg_Block).xy.bb.setZero();
  (*gg_Block).xz.aa.setZero();
  (*gg_Block).xz.ab.setZero();
  (*gg_Block).xz.ba.setZero();
  (*gg_Block).xz.bb.setZero();
  (*gg_Block).yy.aa.setZero();
  (*gg_Block).yy.ab.setZero();
  (*gg_Block).yy.ba.setZero();
  (*gg_Block).yy.bb.setZero();
  (*gg_Block).yz.aa.setZero();
  (*gg_Block).yz.ab.setZero();
  (*gg_Block).yz.ba.setZero();
  (*gg_Block).yz.bb.setZero();
  (*gg_Block).zz.aa.setZero();
  (*gg_Block).zz.ab.setZero();
  (*gg_Block).zz.ba.setZero();
  (*gg_Block).zz.bb.setZero();

  if (_systems.size() > 1) {
    (*gg_Block).xx.aa += (*_ggtot).xx.aa.segment(iGridStart, blockSize);
    (*gg_Block).xx.ab += (*_ggtot).xx.ab.segment(iGridStart, blockSize);
    (*gg_Block).xx.ba += (*_ggtot).xx.ba.segment(iGridStart, blockSize);
    (*gg_Block).xx.bb += (*_ggtot).xx.bb.segment(iGridStart, blockSize);
    (*gg_Block).xy.aa += (*_ggtot).xy.aa.segment(iGridStart, blockSize);
    (*gg_Block).xy.ab += (*_ggtot).xy.ab.segment(iGridStart, blockSize);
    (*gg_Block).xy.ba += (*_ggtot).xy.ba.segment(iGridStart, blockSize);
    (*gg_Block).xy.bb += (*_ggtot).xy.bb.segment(iGridStart, blockSize);
    (*gg_Block).xz.aa += (*_ggtot).xz.aa.segment(iGridStart, blockSize);
    (*gg_Block).xz.ab += (*_ggtot).xz.ab.segment(iGridStart, blockSize);
    (*gg_Block).xz.ba += (*_ggtot).xz.ba.segment(iGridStart, blockSize);
    (*gg_Block).xz.bb += (*_ggtot).xz.bb.segment(iGridStart, blockSize);
    (*gg_Block).yy.aa += (*_ggtot).yy.aa.segment(iGridStart, blockSize);
    (*gg_Block).yy.ab += (*_ggtot).yy.ab.segment(iGridStart, blockSize);
    (*gg_Block).yy.ba += (*_ggtot).yy.ba.segment(iGridStart, blockSize);
    (*gg_Block).yy.bb += (*_ggtot).yy.bb.segment(iGridStart, blockSize);
    (*gg_Block).yz.aa += (*_ggtot).yz.aa.segment(iGridStart, blockSize);
    (*gg_Block).yz.ab += (*_ggtot).yz.ab.segment(iGridStart, blockSize);
    (*gg_Block).yz.ba += (*_ggtot).yz.ba.segment(iGridStart, blockSize);
    (*gg_Block).yz.bb += (*_ggtot).yz.bb.segment(iGridStart, blockSize);
    (*gg_Block).zz.aa += (*_ggtot).zz.aa.segment(iGridStart, blockSize);
    (*gg_Block).zz.ab += (*_ggtot).zz.ab.segment(iGridStart, blockSize);
    (*gg_Block).zz.ba += (*_ggtot).zz.ba.segment(iGridStart, blockSize);
    (*gg_Block).zz.bb += (*_ggtot).zz.bb.segment(iGridStart, blockSize);
  }
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& ggsub = _gg.find(name)->second;
    (*gg_Block).xx.aa += ggsub.xx.aa.segment(iGridStart, blockSize);
    (*gg_Block).xx.ab += ggsub.xx.ab.segment(iGridStart, blockSize);
    (*gg_Block).xx.ba += ggsub.xx.ba.segment(iGridStart, blockSize);
    (*gg_Block).xx.bb += ggsub.xx.bb.segment(iGridStart, blockSize);
    (*gg_Block).xy.aa += ggsub.xy.aa.segment(iGridStart, blockSize);
    (*gg_Block).xy.ab += ggsub.xy.ab.segment(iGridStart, blockSize);
    (*gg_Block).xy.ba += ggsub.xy.ba.segment(iGridStart, blockSize);
    (*gg_Block).xy.bb += ggsub.xy.bb.segment(iGridStart, blockSize);
    (*gg_Block).xz.aa += ggsub.xz.aa.segment(iGridStart, blockSize);
    (*gg_Block).xz.ab += ggsub.xz.ab.segment(iGridStart, blockSize);
    (*gg_Block).xz.ba += ggsub.xz.ba.segment(iGridStart, blockSize);
    (*gg_Block).xz.bb += ggsub.xz.bb.segment(iGridStart, blockSize);
    (*gg_Block).yy.aa += ggsub.yy.aa.segment(iGridStart, blockSize);
    (*gg_Block).yy.ab += ggsub.yy.ab.segment(iGridStart, blockSize);
    (*gg_Block).yy.ba += ggsub.yy.ba.segment(iGridStart, blockSize);
    (*gg_Block).yy.bb += ggsub.yy.bb.segment(iGridStart, blockSize);
    (*gg_Block).yz.aa += ggsub.yz.aa.segment(iGridStart, blockSize);
    (*gg_Block).yz.ab += ggsub.yz.ab.segment(iGridStart, blockSize);
    (*gg_Block).yz.ba += ggsub.yz.ba.segment(iGridStart, blockSize);
    (*gg_Block).yz.bb += ggsub.yz.bb.segment(iGridStart, blockSize);
    (*gg_Block).zz.aa += ggsub.zz.aa.segment(iGridStart, blockSize);
    (*gg_Block).zz.ab += ggsub.zz.ab.segment(iGridStart, blockSize);
    (*gg_Block).zz.ba += ggsub.zz.ba.segment(iGridStart, blockSize);
    (*gg_Block).zz.bb += ggsub.zz.bb.segment(iGridStart, blockSize);
  }
  return gg_Block;
}
template<>
void Kernel<Options::SCF_MODES::RESTRICTED>::storeDerivatives(
    FunctionalData<Options::SCF_MODES::RESTRICTED>& funcData, DensityOnGrid<Options::SCF_MODES::RESTRICTED>& p,
    CompositeFunctionals::CLASSES funcClass,
    DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, GridData<Options::SCF_MODES::RESTRICTED>>& pp,
    Hessian<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, GridData<Options::SCF_MODES::RESTRICTED>>>& gg,
    Gradient<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, GridData<Options::SCF_MODES::RESTRICTED>>>& pg, int pm) {
  // Note: Keys are already present upon initialization. No need to check whether key is present or not.
  pp += Eigen::VectorXd(pm * (*funcData.d2FdRho2));
  if (funcClass == CompositeFunctionals::CLASSES::GGA) {
    pg.x += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).x);
    pg.y += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).y);
    pg.z += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).z);
    gg.xx += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xx);
    gg.xy += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xy);
    gg.xz += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xz);
    gg.yy += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yy);
    gg.yz += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yz);
    gg.zz += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).zz);
  }
  // Screen Kernel:
  // If the density is close to zero, p^-x goes to infinity and will cause numerical
  // problems. However, if the density at this grid points is low, also the molecular
  // orbitals at this point will be small so that the contribution to the XC kernel
  // integrals can be neglected. Since the kernel at this points goes to infinity,
  // while the molecular orbitals go to zero, it is not guaranteed that the integrals
  // will numerically be zero. Thus, we set the kernel to zero at positions
  // where the density falls below a threshold.
  // ToDo: Make threshold input parameter
  double screeningThreshold = 1.0e-8;
  for (unsigned int iPoint = 0; iPoint < _gridController->getNGridPoints(); ++iPoint) {
    if (p[iPoint] < screeningThreshold) {
      pp(iPoint) = 0.0;
      if (_gga) {
        pg.x(iPoint) = 0.0;
        pg.y(iPoint) = 0.0;
        pg.z(iPoint) = 0.0;
        gg.xx(iPoint) = 0.0;
        gg.xy(iPoint) = 0.0;
        gg.xz(iPoint) = 0.0;
        gg.yy(iPoint) = 0.0;
        gg.yz(iPoint) = 0.0;
        gg.zz(iPoint) = 0.0;
      }
    }
  }
}

template<>
void Kernel<Options::SCF_MODES::UNRESTRICTED>::storeDerivatives(
    FunctionalData<Options::SCF_MODES::UNRESTRICTED>& funcData, DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>& p,
    CompositeFunctionals::CLASSES funcClass,
    DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<Options::SCF_MODES::RESTRICTED>>& pp,
    Hessian<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<Options::SCF_MODES::RESTRICTED>>>& gg,
    Gradient<DoublySpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, GridData<Options::SCF_MODES::RESTRICTED>>>& pg,
    int pm) {
  // Note: Keys are already present upon initialization. No need to check whether key is present or not.
  pp.aa += Eigen::VectorXd(pm * (*funcData.d2FdRho2).aa);
  pp.ab += Eigen::VectorXd(pm * (*funcData.d2FdRho2).ab);
  pp.ba += Eigen::VectorXd(pm * (*funcData.d2FdRho2).ab);
  pp.bb += Eigen::VectorXd(pm * (*funcData.d2FdRho2).bb);
  if (funcClass == CompositeFunctionals::CLASSES::GGA) {
    pg.x.aa += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).x.aa);
    pg.x.ab += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).x.ab);
    pg.x.ba += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).x.ba);
    pg.x.bb += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).x.bb);
    pg.y.aa += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).y.aa);
    pg.y.ab += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).y.ab);
    pg.y.ba += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).y.ba);
    pg.y.bb += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).y.bb);
    pg.z.aa += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).z.aa);
    pg.z.ab += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).z.ab);
    pg.z.ba += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).z.ba);
    pg.z.bb += Eigen::VectorXd(pm * (*funcData.d2FdRhodGradRho).z.bb);
    gg.xx.aa += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xx.aa);
    gg.xx.ab += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xx.ab);
    gg.xx.ba += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xx.ab);
    gg.xx.bb += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xx.bb);
    gg.xy.aa += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xy.aa);
    gg.xy.ab += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xy.ab);
    gg.xy.ba += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xy.ab);
    gg.xy.bb += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xy.bb);
    gg.xz.aa += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xz.aa);
    gg.xz.ab += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xz.ab);
    gg.xz.ba += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xz.ab);
    gg.xz.bb += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).xz.bb);
    gg.yy.aa += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yy.aa);
    gg.yy.ab += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yy.ab);
    gg.yy.ba += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yy.ab);
    gg.yy.bb += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yy.bb);
    gg.yz.aa += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yz.aa);
    gg.yz.ab += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yz.ab);
    gg.yz.ba += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yz.ab);
    gg.yz.bb += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).yz.bb);
    gg.zz.aa += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).zz.aa);
    gg.zz.ab += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).zz.ab);
    gg.zz.ba += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).zz.ab);
    gg.zz.bb += Eigen::VectorXd(pm * (*funcData.d2FdGradRho2).zz.bb);
  }
  // Screen Kernel: For comment, see restricted function
  double screeningThreshold = 1.0e-8;
  for (unsigned int iPoint = 0; iPoint < _gridController->getNGridPoints(); ++iPoint) {
    if (p.alpha(iPoint) < screeningThreshold) {
      pp.aa(iPoint) = 0.0;
      pp.ab(iPoint) = 0.0;
      pp.ba(iPoint) = 0.0;
      if (_gga) {
        pg.x.aa(iPoint) = 0.0;
        pg.x.ab(iPoint) = 0.0;
        pg.x.ba(iPoint) = 0.0;
        pg.y.aa(iPoint) = 0.0;
        pg.y.ab(iPoint) = 0.0;
        pg.y.ba(iPoint) = 0.0;
        pg.z.aa(iPoint) = 0.0;
        pg.z.ab(iPoint) = 0.0;
        pg.z.ba(iPoint) = 0.0;
        gg.xx.aa(iPoint) = 0.0;
        gg.xx.ab(iPoint) = 0.0;
        gg.xx.ba(iPoint) = 0.0;
        gg.xy.aa(iPoint) = 0.0;
        gg.xy.ab(iPoint) = 0.0;
        gg.xy.ba(iPoint) = 0.0;
        gg.xz.aa(iPoint) = 0.0;
        gg.xz.ab(iPoint) = 0.0;
        gg.xz.ba(iPoint) = 0.0;
        gg.yy.aa(iPoint) = 0.0;
        gg.yy.ab(iPoint) = 0.0;
        gg.yy.ba(iPoint) = 0.0;
        gg.yz.aa(iPoint) = 0.0;
        gg.yz.ab(iPoint) = 0.0;
        gg.yz.ba(iPoint) = 0.0;
        gg.zz.aa(iPoint) = 0.0;
        gg.zz.ab(iPoint) = 0.0;
        gg.zz.ba(iPoint) = 0.0;
      }
    }
    if (p.beta(iPoint) < screeningThreshold) {
      pp.ab(iPoint) = 0.0;
      pp.ba(iPoint) = 0.0;
      pp.bb(iPoint) = 0.0;
      if (_gga) {
        pg.x.ab(iPoint) = 0.0;
        pg.x.ba(iPoint) = 0.0;
        pg.x.bb(iPoint) = 0.0;
        pg.y.ab(iPoint) = 0.0;
        pg.y.ba(iPoint) = 0.0;
        pg.y.bb(iPoint) = 0.0;
        pg.z.ab(iPoint) = 0.0;
        pg.z.ba(iPoint) = 0.0;
        pg.z.bb(iPoint) = 0.0;
        gg.xx.ab(iPoint) = 0.0;
        gg.xx.ba(iPoint) = 0.0;
        gg.xx.bb(iPoint) = 0.0;
        gg.xy.aa(iPoint) = 0.0;
        gg.xy.ab(iPoint) = 0.0;
        gg.xy.ba(iPoint) = 0.0;
        gg.xy.bb(iPoint) = 0.0;
        gg.xz.ab(iPoint) = 0.0;
        gg.xz.ba(iPoint) = 0.0;
        gg.xz.bb(iPoint) = 0.0;
        gg.yy.ab(iPoint) = 0.0;
        gg.yy.ba(iPoint) = 0.0;
        gg.yy.bb(iPoint) = 0.0;
        gg.yz.ab(iPoint) = 0.0;
        gg.yz.ba(iPoint) = 0.0;
        gg.yz.bb(iPoint) = 0.0;
        gg.zz.ab(iPoint) = 0.0;
        gg.zz.ba(iPoint) = 0.0;
        gg.zz.bb(iPoint) = 0.0;
      }
    }
  }
}

template<Options::SCF_MODES T>
void Kernel<T>::calculateDerivatives() {
  // Combine functionals and build new one for total density contributions
  Functional EnaddXC = resolveFunctional(_naddXCFunc);
  Functional EnaddKIN = resolveFunctional(_naddKinFunc);
  std::vector<BasicFunctionals::BASIC_FUNCTIONALS> naddBasicFunctionals;
  std::vector<double> naddMixingFactors;
  // XCFun gets segmentation fault if functionals are provided in the order GGA NONE?!
  // If the functionals are provided in the order NONE GGA, everything is fine; the order
  // LDA NONE is also fine.
  // To prevent this segmentation fault, none functionals are excluded if a functional from
  // a different functional class is used. If both EnaddXC and EnaddKin are of functional class
  // none, continue with none.
  //
  // As Sept 7, 2020 this should be outdated and could be resolved much cleaner
  // This is to be checked - JU
  if (EnaddXC.getFunctionalClass() == CompositeFunctionals::CLASSES::NONE &&
      EnaddKIN.getFunctionalClass() == CompositeFunctionals::CLASSES::NONE) {
    naddBasicFunctionals.push_back(BasicFunctionals::BASIC_FUNCTIONALS::NONE);
    naddMixingFactors.push_back(0.0);
  }
  else {
    for (unsigned int i = 0; i < EnaddKIN.getNBasicFunctionals(); ++i) {
      if (BasicFunctionals::getClass[(int)EnaddKIN.getBasicFunctionals()[i]] == BasicFunctionals::CLASSES::NONE)
        continue;
      naddBasicFunctionals.push_back(EnaddKIN.getBasicFunctionals()[i]);
      naddMixingFactors.push_back(EnaddKIN.getMixingFactors()[i]);
    }
    if (!EnaddXC.isRSHybrid() && !EnaddXC.isHybrid()) {
      // Special treatment for hybrid functionals
      for (unsigned int i = 0; i < EnaddXC.getNBasicFunctionals(); ++i) {
        if (BasicFunctionals::getClass[(int)EnaddXC.getBasicFunctionals()[i]] == BasicFunctionals::CLASSES::NONE)
          continue;
        naddBasicFunctionals.push_back(EnaddXC.getBasicFunctionals()[i]);
        naddMixingFactors.push_back(EnaddXC.getMixingFactors()[i]);
      }
    }
    if (naddBasicFunctionals.empty())
      naddBasicFunctionals.push_back(BasicFunctionals::BASIC_FUNCTIONALS::NONE);
    if (naddMixingFactors.empty())
      naddMixingFactors.push_back(0.0);
  }
  Functional totFunc(EnaddXC.implementation(), naddBasicFunctionals, naddMixingFactors);
  // Prepare total density and gradient
  std::unique_ptr<DensityOnGrid<T>> totalDensity = std::unique_ptr<DensityOnGrid<T>>(new DensityOnGrid<T>(_gridController));
  std::shared_ptr<Gradient<DensityOnGrid<T>>> totalDensityGradient = makeGradientPtr<DensityOnGrid<T>>(_gridController);
  for (auto sys : _systems) {
    // Combine functionals and build new one for active density contributions
    std::vector<BasicFunctionals::BASIC_FUNCTIONALS> sysBasicFunctionals;
    std::vector<double> sysMixingFacotrs;
    // Subtract non-additive functionals (if any)
    for (unsigned int i = 0; i < totFunc.getNBasicFunctionals(); ++i) {
      sysBasicFunctionals.push_back(totFunc.getBasicFunctionals()[i]);
      sysMixingFacotrs.push_back(-1.0 * totFunc.getMixingFactors()[i]);
    }
    Functional EXC = resolveFunctional(sys->getSettings().dft.functional);
    if (_func != CompositeFunctionals::XCFUNCTIONALS::NONE)
      EXC = resolveFunctional(_func);
    if (!EXC.isRSHybrid() && !EXC.isHybrid()) {
      // hybrids need separate treatment
      for (unsigned int i = 0; i < EXC.getNBasicFunctionals(); ++i) {
        sysBasicFunctionals.push_back(EXC.getBasicFunctionals()[i]);
        sysMixingFacotrs.push_back(EXC.getMixingFactors()[i]);
      }
    }
    Functional sysFunc(totFunc.implementation(), sysBasicFunctionals, sysMixingFacotrs);
    // Calculate density of sys and add to total density
    auto densityOnGridController = getDensityOnGridController(sys);
    auto density = std::unique_ptr<DensityOnGrid<T>>(new DensityOnGrid<T>(_gridController));
    (*density) = densityOnGridController->getDensityOnGrid();
    (*totalDensity) += (*density);
    // Calculate density gradient of sys and add to total density gradient
    std::shared_ptr<Gradient<DensityOnGrid<T>>> densityGradient = nullptr;
    densityGradient = makeGradientPtr<DensityOnGrid<T>>(_gridController);
    (*densityGradient) = densityOnGridController->getDensityGradientOnGrid();
    (*totalDensityGradient) += (*densityGradient);
    // Calculate derivatives of sys functional and add to storage objects
    XCFun<T> xcfun(sys->getSettings().grid.blocksize);
    auto funcData = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, sysFunc, densityOnGridController, 2);
    auto name = sys->getSettings().name;
    storeDerivatives(funcData, (*density), sysFunc.getFunctionalClass(), _pp.find(name)->second, _gg.find(name)->second,
                     _pg.find(name)->second);
    if (EXC.isRSHybrid() || EXC.isHybrid()) {
      // Separate treatment for hybrids
      auto funcDataHybrid = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, EXC, densityOnGridController, 2);
      storeDerivatives(funcDataHybrid, (*density), EXC.getFunctionalClass(), _pp.find(name)->second,
                       _gg.find(name)->second, _pg.find(name)->second);
    }

    if (EnaddXC.isRSHybrid() || EnaddXC.isHybrid()) {
      // Separate treatment for hybrids
      auto funcDataHybrid = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, EnaddXC, densityOnGridController, 2);
      storeDerivatives(funcDataHybrid, (*density), EnaddXC.getFunctionalClass(), _pp.find(name)->second,
                       _gg.find(name)->second, _pg.find(name)->second, -1);
    }
  }
  // All done for conventional TDDFT
  if (_systems.size() == 1)
    return;
  // If none functionals are used, all is done
  if (EnaddKIN.getFunctionalClass() == CompositeFunctionals::CLASSES::NONE &&
      EnaddXC.getFunctionalClass() == CompositeFunctionals::CLASSES::NONE)
    return;
  // Build DensityOnGridController for total density
  auto totDens = makeGradientPtr<DensityOnGrid<T>>(_gridController);
  (*totDens) = (*totalDensityGradient);
  std::shared_ptr<ExternalDensityOnGridController<T>> totalDensityOnGridController =
      std::make_shared<ExternalDensityOnGridController<T>>(totalDensity, totDens);
  // Calculate derivatives of totFunc and add to storage objects
  XCFun<T> xcfun(_systems[0]->getSettings().grid.blocksize);
  auto funcData = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, totFunc, totalDensityOnGridController, 2);
  DensityOnGrid<T> tmp(_gridController);
  tmp = totalDensityOnGridController->getDensityOnGrid();
  storeDerivatives(funcData, tmp, totFunc.getFunctionalClass(), (*_pptot), (*_ggtot), (*_pgtot));

  if (EnaddXC.isRSHybrid() || EnaddXC.isHybrid()) {
    // Separate treatment for hybrids
    auto funcDataHybrid = xcfun.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, EnaddXC, totalDensityOnGridController, 2);
    storeDerivatives(funcDataHybrid, tmp, EnaddXC.getFunctionalClass(), (*_pptot), (*_ggtot), (*_pgtot));
  }
}

template<>
std::shared_ptr<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>>
Kernel<Options::SCF_MODES::RESTRICTED>::getDensityOnGridController(std::shared_ptr<SystemController> sys) {
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(sys->getSettings(), sys->getBasisController(), _gridController);
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>>(
      basisFunctionOnGridController, sys->getSettings().grid.blockAveThreshold);
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::RESTRICTED>> densityMatrixController = nullptr;
  if (sys->getLastSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
    // ToDo
    throw SerenityError("UNRESTRICTED ENVIRONMENT FOR RESTRICTED FDE-TDDFT CALCULATIONS NOT IMPLEMENTED");
  }
  else {
    densityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::RESTRICTED>>(
        sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>(),
        sys->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>());
  }
  return std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>>(densityOnGridCalculator,
                                                                                                densityMatrixController);
}

template<>
std::shared_ptr<DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>
Kernel<Options::SCF_MODES::UNRESTRICTED>::getDensityOnGridController(std::shared_ptr<SystemController> sys) {
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(sys->getSettings(), sys->getBasisController(), _gridController);
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::UNRESTRICTED>>(
      basisFunctionOnGridController, sys->getSettings().grid.blockAveThreshold);
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> densityMatrixController = nullptr;
  // If restricted environment system is used, build unrestricted DensityMatrixController
  if (sys->getLastSCFMode() == Options::SCF_MODES::RESTRICTED) {
    WarningTracker::printWarning(
        (std::string) "WARNING: USE OF RESTRICTED ENVIRONMENT SYSTEM IN UNRESTRICTED EMBEDDING RESPONSE CALCULATION", true);
    // Build unrestricted CoefficientMatrix from restricted one
    std::unique_ptr<CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>> c =
        std::unique_ptr<CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>>(
            new CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>(sys->getBasisController()));
    (*c).alpha = sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
    (*c).beta = sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
    // Build unrestricted eigenvalues from restricted ones
    std::unique_ptr<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>> e =
        std::unique_ptr<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>(
            new SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>);
    (*e).alpha = sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getEigenvalues();
    (*e).beta = sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getEigenvalues();
    // Build unrestricted OrbitalController from restricted data
    std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED>> orbitalController =
        std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(std::move(c), sys->getBasisController(),
                                                                              std::move(e));
    // Build DensityMatrixController
    densityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(
        orbitalController, sys->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>());
  }
  else {
    densityMatrixController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(
        sys->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>(),
        sys->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>());
  }
  return std::make_shared<DensityMatrixDensityOnGridController<Options::SCF_MODES::UNRESTRICTED>>(densityOnGridCalculator,
                                                                                                  densityMatrixController);
}
template<Options::SCF_MODES SCFMode>
unsigned int Kernel<SCFMode>::getBlocksize(unsigned int subsystemNumber) {
  return _systems[subsystemNumber]->getSettings().grid.blocksize;
};
template<Options::SCF_MODES SCFMode>
unsigned int Kernel<SCFMode>::getbasFuncRadialThreshold(unsigned int subsystemNumber) {
  return _systems[subsystemNumber]->getSettings().grid.basFuncRadialThreshold;
};
template<Options::SCF_MODES SCFMode>
double Kernel<SCFMode>::getblockAveThreshold(unsigned int subsystemNumber) {
  return _systems[subsystemNumber]->getSettings().grid.blockAveThreshold;
}

template class Kernel<Options::SCF_MODES::RESTRICTED>;
template class Kernel<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
