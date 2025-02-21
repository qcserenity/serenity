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
#include "dft/Functional.h"
#include "dft/functionals/BasicFunctionals.h"
#include "dft/functionals/CompositeFunctionals.h"
#include "dft/functionals/FunctionalLibrary.h"
#include "dft/functionals/wrappers/XCFun.h"
#include "geometry/Geometry.h"
#include "geometry/GeometryAdderFactory.h"
#include "grid/AtomCenteredGridControllerFactory.h"
#include "io/FormattedOutputStream.h" //Filtered output streams
#include "misc/SerenityError.h"
#include "misc/WarningTracker.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
Kernel<SCFMode>::Kernel(std::vector<std::shared_ptr<SystemController>> act, std::vector<std::shared_ptr<SystemController>> env,
                        const LRSCFTaskSettings& settings, unsigned nLRSCF, Options::GRID_PURPOSES gridAccuracy)
  : _settings(settings),
    _naddKinFunc(settings.embedding.naddKinFunc),
    _naddXCFunc(settings.embedding.naddXCFunc),
    _gga(false),
    _nLRSCF(nLRSCF) {
  // Build system
  for (auto sys : act) {
    _systems.push_back(sys);
  }
  for (auto sys : env) {
    _systems.push_back(sys);
  }

  //  Build reference geometry.
  std::shared_ptr<Geometry> geometry = GeometryAdderFactory::produce(_systems, _settings.subsystemgrid);
  _gridController = AtomCenteredGridControllerFactory::produce(geometry, _settings.grid, gridAccuracy);

  // Get and check data from input functionals
  for (auto sys : _systems) {
    // HF when system uses HF and no other Kernel functional is specified
    if (sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF &&
        settings.func == CompositeFunctionals::XCFUNCTIONALS::NONE && !(settings.customFunc.basicFunctionals.size())) {
      _func.push_back(resolveFunctional(CompositeFunctionals::XCFUNCTIONALS::HF));
    }
    // LRSCFTask before system settings and custom before composite
    // use the custom LRSCF functional if there is one
    else if (settings.customFunc.basicFunctionals.size()) {
      _func.push_back(Functional(settings.customFunc));
    }
    // use the composite LRSCF functionl if there is one
    else if (settings.func != CompositeFunctionals::XCFUNCTIONALS::NONE) {
      _func.push_back(resolveFunctional(settings.func));
    }
    // no LRSCF functional -> use the system settings
    else if (sys->getSettings().customFunc.basicFunctionals.size()) {
      _func.push_back(Functional(sys->getSettings().customFunc));
    }
    // system settings composite functional at last
    else {
      _func.push_back(resolveFunctional(sys->getSettings().dft.functional));
    }
  }

  if (_systems.size() == 1 || settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::NONE) {
    if (settings.embedding.embeddingModeList.empty()) {
      _naddKinFunc = CompositeFunctionals::KINFUNCTIONALS::NONE;
      _naddXCFunc = CompositeFunctionals::XCFUNCTIONALS::NONE;
    }
  }

  Functional EnaddXC = settings.embedding.customNaddXCFunc.basicFunctionals.size()
                           ? Functional(settings.embedding.customNaddXCFunc)
                           : resolveFunctional(_naddXCFunc);
  Functional EnaddKIN = settings.embedding.customNaddKinFunc.basicFunctionals.size()
                            ? Functional(settings.embedding.customNaddKinFunc)
                            : resolveFunctional(_naddKinFunc);
  if (EnaddKIN.isHybrid() || EnaddKIN.isRSHybrid()) {
    throw SerenityError("ERROR: Found hybrid functional for non-additive kinetic kernel.");
  }
  if (EnaddXC.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
    _gga = true;
  }
  if (EnaddKIN.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
    _gga = true;
  }

  for (auto funcXC : _func) {
    if (funcXC.getFunctionalClass() == CompositeFunctionals::CLASSES::GGA) {
      _gga = true;
    }
  }

  // Initialize storage objects (this ensures that key is already present when working with it)
  for (auto sys : _systems) {
    auto sysname = sys->getSystemName();
    _pp.insert(std::make_pair(
        sysname, DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>(_gridController)));
    if (_gga) {
      _pg.insert(std::make_pair(
          sysname, makeGradient<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController)));
      _gg.insert(std::make_pair(
          sysname, makeHessian<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController)));
    }
  }
  if (_systems.size() > 1) {
    _pptot = std::make_shared<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
    _ggtot = makeHessianPtr<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
    _pgtot = makeGradientPtr<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
  }

  // Check if same density keyword is used correctly
  if (_settings.samedensity.size() > 0) {
    if (_settings.samedensity.size() != _systems.size()) {
      throw SerenityError("You need to specify the same amount of systems for the samedensity keyword as systems in "
                          "your system (act + env)!");
    }
  }
  // Check for mixed exact approx embedding
  if (settings.embedding.embeddingModeList.size() == 0) {
    this->calculateDerivatives();
  }
  else {
    auto embeddingmodes = settings.embedding.embeddingModeList;
    if (embeddingmodes.size() != _systems.size()) {
      throw SerenityError("You need to specify the same amount of systems as embedding modes!");
    }
    _mixedEmbeddingUsed = true;
    // Additonal contributions for exact Embedding
    _ppExact = std::make_shared<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
    _ggExact = makeHessianPtr<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
    _pgExact = makeGradientPtr<DoublySpinPolarizedData<SCFMode, GridData<Options::SCF_MODES::RESTRICTED>>>(_gridController);
    // Calculate and store mixed derivatives
    this->calculateDerivativesMixedEmbedding();
  }
}

template<>
const std::shared_ptr<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>
Kernel<Options::SCF_MODES::RESTRICTED>::getPP(unsigned int I, unsigned int J, unsigned int blockSize, unsigned int iGridStart) {
  auto pp_Block = std::make_shared<DoublySpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>>(blockSize);
  (*pp_Block).setZero();
  if (_systems.size() > 1 &&
      (I != J || _nLRSCF == 1 || _mixedEmbeddingUsed || _settings.partialResponseConstruction || _settings.noCoupling)) {
    *pp_Block += (*_pptot).segment(iGridStart, blockSize);
  }
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& pp_system = _pp.find(name)->second;
    *pp_Block += pp_system.segment(iGridStart, blockSize);
  }
  if (_mixedEmbeddingUsed) {
    if ((_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) ||
        (_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA)) {
      *pp_Block += (*_ppExact).segment(iGridStart, blockSize);
    }
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
  if (_systems.size() > 1 &&
      (I != J || _nLRSCF == 1 || _mixedEmbeddingUsed || _settings.partialResponseConstruction || _settings.noCoupling)) {
    (*pp_Block).aa += (*_pptot).aa.segment(iGridStart, blockSize);
    (*pp_Block).ab += (*_pptot).ab.segment(iGridStart, blockSize);
    (*pp_Block).ba += (*_pptot).ba.segment(iGridStart, blockSize);
    (*pp_Block).bb += (*_pptot).bb.segment(iGridStart, blockSize);
  }
  if (I == J) {
    auto name = _systems[I]->getSettings().name;
    auto& ppsub = _pp.find(name)->second;
    (*pp_Block).aa += ppsub.aa.segment(iGridStart, blockSize);
    (*pp_Block).ab += ppsub.ab.segment(iGridStart, blockSize);
    (*pp_Block).ba += ppsub.ba.segment(iGridStart, blockSize);
    (*pp_Block).bb += ppsub.bb.segment(iGridStart, blockSize);
  }
  if (_mixedEmbeddingUsed) {
    if ((_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) ||
        (_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA)) {
      (*pp_Block).aa += (*_ppExact).aa.segment(iGridStart, blockSize);
      (*pp_Block).ab += (*_ppExact).ab.segment(iGridStart, blockSize);
      (*pp_Block).ba += (*_ppExact).ba.segment(iGridStart, blockSize);
      (*pp_Block).bb += (*_ppExact).bb.segment(iGridStart, blockSize);
    }
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
  if (_systems.size() > 1 &&
      (I != J || _nLRSCF == 1 || _mixedEmbeddingUsed || _settings.partialResponseConstruction || _settings.noCoupling)) {
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
  if (_mixedEmbeddingUsed) {
    if ((_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) ||
        (_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA)) {
      (*pg_Block).x += (*_pgExact).x.segment(iGridStart, blockSize);
      (*pg_Block).y += (*_pgExact).y.segment(iGridStart, blockSize);
      (*pg_Block).z += (*_pgExact).z.segment(iGridStart, blockSize);
    }
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

  if (_systems.size() > 1 &&
      (I != J || _nLRSCF == 1 || _mixedEmbeddingUsed || _settings.partialResponseConstruction || _settings.noCoupling)) {
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
  if (_mixedEmbeddingUsed) {
    if ((_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) ||
        (_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA)) {
      (*pg_Block).x.aa += (*_pgExact).x.aa.segment(iGridStart, blockSize);
      (*pg_Block).x.ab += (*_pgExact).x.ab.segment(iGridStart, blockSize);
      (*pg_Block).x.ba += (*_pgExact).x.ba.segment(iGridStart, blockSize);
      (*pg_Block).x.bb += (*_pgExact).x.bb.segment(iGridStart, blockSize);
      (*pg_Block).y.aa += (*_pgExact).y.aa.segment(iGridStart, blockSize);
      (*pg_Block).y.ab += (*_pgExact).y.ab.segment(iGridStart, blockSize);
      (*pg_Block).y.ba += (*_pgExact).y.ba.segment(iGridStart, blockSize);
      (*pg_Block).y.bb += (*_pgExact).y.bb.segment(iGridStart, blockSize);
      (*pg_Block).z.aa += (*_pgExact).z.aa.segment(iGridStart, blockSize);
      (*pg_Block).z.ab += (*_pgExact).z.ab.segment(iGridStart, blockSize);
      (*pg_Block).z.ba += (*_pgExact).z.ba.segment(iGridStart, blockSize);
      (*pg_Block).z.bb += (*_pgExact).z.bb.segment(iGridStart, blockSize);
    }
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

  if (_systems.size() > 1 &&
      (I != J || _nLRSCF == 1 || _mixedEmbeddingUsed || _settings.partialResponseConstruction || _settings.noCoupling)) {
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
  if (_mixedEmbeddingUsed) {
    if ((_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) ||
        (_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA)) {
      (*gg_Block).xx += (*_ggExact).xx.segment(iGridStart, blockSize);
      (*gg_Block).xy += (*_ggExact).xy.segment(iGridStart, blockSize);
      (*gg_Block).xz += (*_ggExact).xz.segment(iGridStart, blockSize);
      (*gg_Block).yy += (*_ggExact).yy.segment(iGridStart, blockSize);
      (*gg_Block).yz += (*_ggExact).yz.segment(iGridStart, blockSize);
      (*gg_Block).zz += (*_ggExact).zz.segment(iGridStart, blockSize);
    }
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

  if (_systems.size() > 1 &&
      (I != J || _nLRSCF == 1 || _mixedEmbeddingUsed || _settings.partialResponseConstruction || _settings.noCoupling)) {
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
  if (_mixedEmbeddingUsed) {
    if ((_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) ||
        (_settings.embedding.embeddingModeList[I] == _settings.embedding.embeddingModeList[J] &&
         _settings.embedding.embeddingModeList[I] == Options::KIN_EMBEDDING_MODES::HUZINAGA)) {
      (*gg_Block).xx.aa += (*_ggExact).xx.aa.segment(iGridStart, blockSize);
      (*gg_Block).xx.ab += (*_ggExact).xx.ab.segment(iGridStart, blockSize);
      (*gg_Block).xx.ba += (*_ggExact).xx.ba.segment(iGridStart, blockSize);
      (*gg_Block).xx.bb += (*_ggExact).xx.bb.segment(iGridStart, blockSize);
      (*gg_Block).xy.aa += (*_ggExact).xy.aa.segment(iGridStart, blockSize);
      (*gg_Block).xy.ab += (*_ggExact).xy.ab.segment(iGridStart, blockSize);
      (*gg_Block).xy.ba += (*_ggExact).xy.ba.segment(iGridStart, blockSize);
      (*gg_Block).xy.bb += (*_ggExact).xy.bb.segment(iGridStart, blockSize);
      (*gg_Block).xz.aa += (*_ggExact).xz.aa.segment(iGridStart, blockSize);
      (*gg_Block).xz.ab += (*_ggExact).xz.ab.segment(iGridStart, blockSize);
      (*gg_Block).xz.ba += (*_ggExact).xz.ba.segment(iGridStart, blockSize);
      (*gg_Block).xz.bb += (*_ggExact).xz.bb.segment(iGridStart, blockSize);
      (*gg_Block).yy.aa += (*_ggExact).yy.aa.segment(iGridStart, blockSize);
      (*gg_Block).yy.ab += (*_ggExact).yy.ab.segment(iGridStart, blockSize);
      (*gg_Block).yy.ba += (*_ggExact).yy.ba.segment(iGridStart, blockSize);
      (*gg_Block).yy.bb += (*_ggExact).yy.bb.segment(iGridStart, blockSize);
      (*gg_Block).yz.aa += (*_ggExact).yz.aa.segment(iGridStart, blockSize);
      (*gg_Block).yz.ab += (*_ggExact).yz.ab.segment(iGridStart, blockSize);
      (*gg_Block).yz.ba += (*_ggExact).yz.ba.segment(iGridStart, blockSize);
      (*gg_Block).yz.bb += (*_ggExact).yz.bb.segment(iGridStart, blockSize);
      (*gg_Block).zz.aa += (*_ggExact).zz.aa.segment(iGridStart, blockSize);
      (*gg_Block).zz.ab += (*_ggExact).zz.ab.segment(iGridStart, blockSize);
      (*gg_Block).zz.ba += (*_ggExact).zz.ba.segment(iGridStart, blockSize);
      (*gg_Block).zz.bb += (*_ggExact).zz.bb.segment(iGridStart, blockSize);
    }
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

template<Options::SCF_MODES SCFMode>
void Kernel<SCFMode>::calculateDerivatives() {
  // total density objects
  auto totalDensity = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_gridController));
  auto totalDensityGradient = makeGradientPtr<DensityOnGrid<SCFMode>>(_gridController);
  Functional exc_nadd = _settings.embedding.customNaddXCFunc.basicFunctionals.size()
                            ? Functional(_settings.embedding.customNaddXCFunc)
                            : resolveFunctional(_naddXCFunc);
  Functional kin_nadd = _settings.embedding.customNaddKinFunc.basicFunctionals.size()
                            ? Functional(_settings.embedding.customNaddKinFunc)
                            : resolveFunctional(_naddKinFunc);
  // Loop over systems
  // Subsystem contributions
  unsigned int counter = 0;
  for (auto sys : _systems) {
    // Returns the density/density gradients on the custom grid
    // Calculate density of sys and add to total density
    auto densityOnGridController = getDensityOnGridController(sys);
    auto density = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_gridController));
    (*density) = densityOnGridController->getDensityOnGrid();
    // Calculate density gradient of sys and add to total density gradient
    std::shared_ptr<Gradient<DensityOnGrid<SCFMode>>> densityGradient = nullptr;
    densityGradient = makeGradientPtr<DensityOnGrid<SCFMode>>(_gridController);
    (*densityGradient) = densityOnGridController->getDensityGradientOnGrid();
    // Summed up density and gradient
    if (_settings.samedensity.size() > 0) {
      if (_settings.samedensity[counter] == (counter + 1)) {
        (*totalDensity) += (*density);
        (*totalDensityGradient) += (*densityGradient);
      }
    }
    else {
      (*totalDensity) += (*density);
      (*totalDensityGradient) += (*densityGradient);
    }
    // Calculate derivatives of sys functional and add to storage objects
    auto name = sys->getSettings().name;
    FunctionalLibrary<SCFMode> flib(sys->getSettings().grid.blocksize);
    auto funcData = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _func[counter], densityOnGridController, 2);
    // Complete kernel contributions
    storeDerivatives(funcData, (*density), _func[counter].getFunctionalClass(), _pp.find(name)->second,
                     _gg.find(name)->second, _pg.find(name)->second);
    // Subsystem specific kernel contributions
    if (exc_nadd.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
      auto funcData_naddxc = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, exc_nadd, densityOnGridController, 2);
      storeDerivatives(funcData_naddxc, (*density), exc_nadd.getFunctionalClass(), _pp.find(name)->second,
                       _gg.find(name)->second, _pg.find(name)->second, -1.0);
    }
    if (kin_nadd.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
      auto funcData_naddkin = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, kin_nadd, densityOnGridController, 2);
      storeDerivatives(funcData_naddkin, (*density), kin_nadd.getFunctionalClass(), _pp.find(name)->second,
                       _gg.find(name)->second, _pg.find(name)->second, -1.0);
    }
    counter++;
  }
  // Naddkin evaluation for complete density
  std::shared_ptr<ExternalDensityOnGridController<SCFMode>> totalDensityOnGridController =
      std::make_shared<ExternalDensityOnGridController<SCFMode>>(totalDensity, totalDensityGradient);
  FunctionalLibrary<SCFMode> flib(_systems[0]->getSettings().grid.blocksize);
  if (exc_nadd.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
    auto funcData_naddxc = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, exc_nadd, totalDensityOnGridController, 2);
    DensityOnGrid<SCFMode> totaldens(_gridController);
    totaldens = totalDensityOnGridController->getDensityOnGrid();
    storeDerivatives(funcData_naddxc, totaldens, exc_nadd.getFunctionalClass(), (*_pptot), (*_ggtot), (*_pgtot));
  }
  if (kin_nadd.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
    auto funcData_naddkin = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, kin_nadd, totalDensityOnGridController, 2);
    DensityOnGrid<SCFMode> totaldens(_gridController);
    totaldens = totalDensityOnGridController->getDensityOnGrid();
    storeDerivatives(funcData_naddkin, totaldens, kin_nadd.getFunctionalClass(), (*_pptot), (*_ggtot), (*_pgtot));
  }
}

template<Options::SCF_MODES SCFMode>
void Kernel<SCFMode>::calculateDerivativesMixedEmbedding() {
  // Embedding modes
  auto embeddingmodes = _settings.embedding.embeddingModeList;
  // Exact and approximate treated systems
  std::vector<std::shared_ptr<SystemController>> treatedExact;
  std::vector<std::shared_ptr<SystemController>> treatedApprox;
  // Total density/gradient for exact and approximated treated subsystems
  std::unique_ptr<DensityOnGrid<SCFMode>> totalDensityForExact =
      std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_gridController));
  std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>> totalDensityForExactGradient =
      makeGradientPtr<DensityOnGrid<SCFMode>>(_gridController);
  std::unique_ptr<DensityOnGrid<SCFMode>> totalDensity =
      std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_gridController));
  std::unique_ptr<Gradient<DensityOnGrid<SCFMode>>> totalDensityGradient =
      makeGradientPtr<DensityOnGrid<SCFMode>>(_gridController);
  // Mixing factors and basic functionals
  std::vector<BasicFunctionals::BASIC_FUNCTIONALS> naddBasicFunctionals_approx;
  std::vector<double> naddMixingFactors_approx;
  std::vector<BasicFunctionals::BASIC_FUNCTIONALS> naddBasicFunctionals_exact;
  std::vector<double> naddMixingFactors_exact;
  // NaddXC functionals
  Functional exc_exact = resolveFunctional(_settings.embedding.naddXCFuncList[0]);
  Functional exc_approx = resolveFunctional(_settings.embedding.naddXCFuncList[1]);
  // Resolve NaddKinfunctional
  Functional kin_nadd = resolveFunctional(_naddKinFunc);
  // Subsystem contributions
  for (unsigned int iSys = 0; iSys < _systems.size(); iSys++) {
    // Density on custom grid
    auto densityOnGridController = getDensityOnGridController(_systems[iSys]);
    auto density = std::unique_ptr<DensityOnGrid<SCFMode>>(new DensityOnGrid<SCFMode>(_gridController));
    (*density) = densityOnGridController->getDensityOnGrid();
    // Calculate density gradient of sys and add to total density gradient
    std::shared_ptr<Gradient<DensityOnGrid<SCFMode>>> densityGradient = nullptr;
    densityGradient = makeGradientPtr<DensityOnGrid<SCFMode>>(_gridController);
    (*densityGradient) = densityOnGridController->getDensityGradientOnGrid();
    // Calculate derivatives of sys functional and add to storage objects
    auto name = _systems[iSys]->getSettings().name;
    FunctionalLibrary<SCFMode> flib(_systems[iSys]->getSettings().grid.blocksize);
    auto funcData = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, _func[iSys], densityOnGridController, 2);
    // Subsystem kernel contributions
    storeDerivatives(funcData, (*density), _func[iSys].getFunctionalClass(), _pp.find(name)->second,
                     _gg.find(name)->second, _pg.find(name)->second);
    if (embeddingmodes[iSys] == Options::KIN_EMBEDDING_MODES::HUZINAGA ||
        embeddingmodes[iSys] == Options::KIN_EMBEDDING_MODES::LEVELSHIFT) {
      // Systems from exact embedding
      treatedExact.push_back(_systems[iSys]);
      // Evaulate contributions for systems with the same density
      if (_settings.samedensity.size() > 0) {
        if (_settings.samedensity[iSys] == (iSys + 1)) {
          (*totalDensity) += (*density);
          (*totalDensityForExact) += (*density);
          (*totalDensityGradient) += (*densityGradient);
          (*totalDensityForExactGradient) += (*densityGradient);
        }
      }
      else {
        (*totalDensity) += (*density);
        (*totalDensityForExact) += (*density);
        (*totalDensityGradient) += (*densityGradient);
        (*totalDensityForExactGradient) += (*densityGradient);
      }
      // Subsystem nadd kernel
      if (exc_exact.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
        auto funcData_naddxc_exact = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, exc_exact, densityOnGridController, 2);
        storeDerivatives(funcData_naddxc_exact, (*density), exc_exact.getFunctionalClass(), _pp.find(name)->second,
                         _gg.find(name)->second, _pg.find(name)->second, -1.0);
      }
    }
    else if (embeddingmodes[iSys] == Options::KIN_EMBEDDING_MODES::FERMI_SHIFTED_HUZINAGA ||
             embeddingmodes[iSys] == Options::KIN_EMBEDDING_MODES::HOFFMANN ||
             embeddingmodes[iSys] == Options::KIN_EMBEDDING_MODES::RECONSTRUCTION) {
      throw SerenityError("Exact embedding mode not supported in list input yet!");
    }
    else {
      // Systems from approx embedding
      treatedApprox.push_back(_systems[iSys]);
      // Evaulate contributions for systems with the same density
      if (_settings.samedensity.size() > 0) {
        if (_settings.samedensity[iSys] == (iSys + 1)) {
          (*totalDensity) += (*density);
          (*totalDensityGradient) += (*densityGradient);
        }
      }
      else {
        (*totalDensity) += (*density);
        (*totalDensityGradient) += (*densityGradient);
      }
      // Subsystem specific kernel contributions
      if (exc_approx.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
        auto funcData_naddxc_approx = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, exc_approx, densityOnGridController, 2);
        storeDerivatives(funcData_naddxc_approx, (*density), exc_approx.getFunctionalClass(), _pp.find(name)->second,
                         _gg.find(name)->second, _pg.find(name)->second, -1.0);
      }
      if (kin_nadd.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
        auto funcData_naddkin = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, kin_nadd, densityOnGridController, 2);
        storeDerivatives(funcData_naddkin, (*density), kin_nadd.getFunctionalClass(), _pp.find(name)->second,
                         _gg.find(name)->second, _pg.find(name)->second, -1.0);
      }
    }
  }
  // System sepcific debug output
  OutputControl::dOut << " -- Mixed Exact-Approx Embedding: -- " << std::endl;
  OutputControl::dOut << "  Systems treated exact - " << treatedExact.size() << std::endl;
  OutputControl::dOut << "  Systems treated approx. - " << treatedApprox.size() << std::endl;

  // Naddkin evaluation for complete density
  std::shared_ptr<ExternalDensityOnGridController<SCFMode>> totalDensityOnGridController =
      std::make_shared<ExternalDensityOnGridController<SCFMode>>(totalDensity, totalDensityGradient);
  FunctionalLibrary<SCFMode> flib(_systems[0]->getSettings().grid.blocksize);
  if (exc_approx.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
    auto funcData_naddxc_approx = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, exc_approx, totalDensityOnGridController, 2);
    DensityOnGrid<SCFMode> totaldens(_gridController);
    totaldens = totalDensityOnGridController->getDensityOnGrid();
    storeDerivatives(funcData_naddxc_approx, totaldens, exc_approx.getFunctionalClass(), (*_pptot), (*_ggtot), (*_pgtot));
  }
  if (kin_nadd.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
    auto funcData_naddkin = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, kin_nadd, totalDensityOnGridController, 2);
    DensityOnGrid<SCFMode> totaldens(_gridController);
    totaldens = totalDensityOnGridController->getDensityOnGrid();
    storeDerivatives(funcData_naddkin, totaldens, kin_nadd.getFunctionalClass(), (*_pptot), (*_ggtot), (*_pgtot));
  }
  // Naddkin evaluation for density of exactly treates systems
  std::shared_ptr<ExternalDensityOnGridController<SCFMode>> totalExactDensityOnGridController =
      std::make_shared<ExternalDensityOnGridController<SCFMode>>(totalDensityForExact, totalDensityForExactGradient);
  if (exc_exact.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
    auto funcData_naddxc_exact =
        flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, exc_exact, totalExactDensityOnGridController, 2);
    DensityOnGrid<SCFMode> totaldens(_gridController);
    totaldens = totalExactDensityOnGridController->getDensityOnGrid();
    storeDerivatives(funcData_naddxc_exact, totaldens, exc_exact.getFunctionalClass(), (*_ppExact), (*_ggExact), (*_pgExact));
  }
  if (exc_approx.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
    auto funcData_naddxc_approx =
        flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, exc_approx, totalExactDensityOnGridController, 2);
    DensityOnGrid<SCFMode> totaldens(_gridController);
    totaldens = totalExactDensityOnGridController->getDensityOnGrid();
    storeDerivatives(funcData_naddxc_approx, totaldens, exc_approx.getFunctionalClass(), (*_ppExact), (*_ggExact),
                     (*_pgExact), -1.0);
  }
  if (kin_nadd.getFunctionalClass() != CompositeFunctionals::CLASSES::NONE) {
    auto funcData_naddkin = flib.calcData(FUNCTIONAL_DATA_TYPE::GRADIENTS, kin_nadd, totalExactDensityOnGridController, 2);
    DensityOnGrid<SCFMode> totaldens(_gridController);
    totaldens = totalExactDensityOnGridController->getDensityOnGrid();
    storeDerivatives(funcData_naddkin, totaldens, kin_nadd.getFunctionalClass(), (*_ppExact), (*_ggExact), (*_pgExact), -1.0);
  }
}

template<>
std::shared_ptr<DensityMatrixDensityOnGridController<Options::SCF_MODES::RESTRICTED>>
Kernel<Options::SCF_MODES::RESTRICTED>::getDensityOnGridController(std::shared_ptr<SystemController> sys) {
  auto basisFunctionOnGridController =
      BasisFunctionOnGridControllerFactory::produce(sys->getSettings(), sys->getBasisController(), _gridController);
  auto densityOnGridCalculator = std::make_shared<DensityOnGridCalculator<Options::SCF_MODES::RESTRICTED>>(
      basisFunctionOnGridController, _settings.grid.blockAveThreshold);
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
      basisFunctionOnGridController, _settings.grid.blockAveThreshold);
  std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> densityMatrixController = nullptr;
  // If restricted environment system is used, build unrestricted DensityMatrixController
  if (sys->getLastSCFMode() == Options::SCF_MODES::RESTRICTED) {
    // Build unrestricted CoefficientMatrix from restricted one
    auto c = std::make_unique<CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>>(
        sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients());
    // Build unrestricted eigenvalues from restricted ones
    auto e = std::make_unique<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>>(
        sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getEigenvalues());
    // Build unrestricted core orbital vector from restricted ones
    auto core = std::make_unique<SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXi>>(
        sys->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getOrbitalFlags());
    // Build unrestricted OrbitalController from restricted data
    std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED>> orbitalController =
        std::make_shared<OrbitalController<Options::SCF_MODES::UNRESTRICTED>>(std::move(c), sys->getBasisController(),
                                                                              std::move(e), std::move(core));
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
unsigned Kernel<SCFMode>::getBlocksize(unsigned I) {
  return _systems[I]->getSettings().grid.blocksize;
};

template<Options::SCF_MODES SCFMode>
unsigned Kernel<SCFMode>::getbasFuncRadialThreshold(unsigned I) {
  return _systems[I]->getSettings().grid.basFuncRadialThreshold;
};

template<Options::SCF_MODES SCFMode>
double Kernel<SCFMode>::getblockAveThreshold(unsigned I) {
  return _systems[I]->getSettings().grid.blockAveThreshold;
}

template<Options::SCF_MODES SCFMode>
unsigned int Kernel<SCFMode>::getNSystems() {
  return _systems.size();
}

template<Options::SCF_MODES SCFMode>
bool Kernel<SCFMode>::usesMixedEmbedding() {
  return _mixedEmbeddingUsed;
}

template class Kernel<Options::SCF_MODES::RESTRICTED>;
template class Kernel<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
