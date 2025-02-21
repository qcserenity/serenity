/**
 * @file DensityMatrixController.cpp
 *
 * @date May 12, 2016
 * @author Jan Unsleber
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
#include "data/matrices/DensityMatrixController.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/HDF5.h"
#include "math/linearAlgebra/MatrixFunctions.h"

namespace Serenity {

template<>
DensityMatrixController<Options::SCF_MODES::RESTRICTED>::DensityMatrixController(
    const std::shared_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>> molecularOrbitals,
    const SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>& nOccupiedOrbitals)
  : _densityMatrix(nullptr),
    _molecularOrbitals(molecularOrbitals),
    _occupations(nullptr),
    _aufbauOccupations(Eigen::VectorXd(molecularOrbitals->getBasisController()->getNBasisFunctions())),
    _basisController(molecularOrbitals->getBasisController()),
    _outOfDate(true) {
  _molecularOrbitals->addSensitiveObject(this->_self);
  _aufbauOccupations.setZero();
  for (unsigned int i = 0; i < nOccupiedOrbitals; ++i) {
    _aufbauOccupations[i] = 2.0;
  }
}

template<>
DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>::DensityMatrixController(
    const std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED>> molecularOrbitals,
    const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int>& nOccupiedOrbitals)
  : _densityMatrix(nullptr),
    _molecularOrbitals(molecularOrbitals),
    _occupations(nullptr),
    _aufbauOccupations(makeUnrestrictedFromPieces<Eigen::VectorXd>(
        Eigen::VectorXd(molecularOrbitals->getBasisController()->getNBasisFunctions()),
        Eigen::VectorXd(molecularOrbitals->getBasisController()->getNBasisFunctions()))),
    _basisController(molecularOrbitals->getBasisController()),
    _outOfDate(true) {
  _molecularOrbitals->addSensitiveObject(this->_self);
  for_spin(nOccupiedOrbitals, _aufbauOccupations) {
    _aufbauOccupations_spin.setZero();
    for (unsigned int i = 0; i < nOccupiedOrbitals_spin; ++i) {
      _aufbauOccupations_spin[i] = 1.0;
    }
  };
}

template<Options::SCF_MODES SCFMode>
DensityMatrixController<SCFMode>::DensityMatrixController(const std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                                                          const SpinPolarizedData<SCFMode, Eigen::VectorXd>& occupations)
  : _densityMatrix(nullptr),
    _molecularOrbitals(molecularOrbitals),
    _occupations(nullptr),
    _aufbauOccupations(occupations),
    _basisController(molecularOrbitals->getBasisController()),
    _outOfDate(true) {
  _molecularOrbitals->addSensitiveObject(this->_self);
}

template<>
DensityMatrixController<RESTRICTED>::DensityMatrixController(const DensityMatrix<RESTRICTED>& densityMatrix)
  : _densityMatrix(new DensityMatrix<RESTRICTED>(densityMatrix)),
    _molecularOrbitals(nullptr),
    _occupations(nullptr),
    _aufbauOccupations(Eigen::VectorXd(densityMatrix.getBasisController()->getNBasisFunctions())),
    _basisController(densityMatrix.getBasisController()),
    _outOfDate(false) {
}

template<>
DensityMatrixController<UNRESTRICTED>::DensityMatrixController(const DensityMatrix<UNRESTRICTED>& densityMatrix)
  : _densityMatrix(new DensityMatrix<UNRESTRICTED>(densityMatrix)),
    _molecularOrbitals(nullptr),
    _occupations(nullptr),
    _aufbauOccupations(makeUnrestrictedFromPieces<Eigen::VectorXd>(
        Eigen::VectorXd(densityMatrix.getBasisController()->getNBasisFunctions()),
        Eigen::VectorXd(densityMatrix.getBasisController()->getNBasisFunctions()))),
    _basisController(densityMatrix.getBasisController()),
    _outOfDate(false) {
}
template<Options::SCF_MODES SCFMode>
DensityMatrixController<SCFMode>::DensityMatrixController(std::string fBaseName,
                                                          std::shared_ptr<BasisController> basisController, std::string id)
  : _densityMatrix(nullptr),
    _molecularOrbitals(new OrbitalController<SCFMode>(fBaseName, basisController, id)),
    _occupations(nullptr),
    _aufbauOccupations(basisController->getNBasisFunctions()),
    _basisController(basisController),
    _outOfDate(false) {
  fromHDF5(fBaseName, id);
  _molecularOrbitals->addSensitiveObject(this->_self);
}

template<>
void DensityMatrixController<Options::SCF_MODES::RESTRICTED>::attachOrbitals(
    std::shared_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>> molecularOrbitals,
    const SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int>& nOccupiedOrbitals, bool update) {
  _molecularOrbitals = molecularOrbitals;
  _molecularOrbitals->addSensitiveObject(this->_self);
  _aufbauOccupations.setZero();
  for (unsigned int i = 0; i < nOccupiedOrbitals; ++i) {
    _aufbauOccupations[i] = 2.0;
  }
  if (update)
    _outOfDate = true;
}

template<>
void DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>::attachOrbitals(
    std::shared_ptr<OrbitalController<Options::SCF_MODES::UNRESTRICTED>> molecularOrbitals,
    const SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int>& nOccupiedOrbitals, bool update) {
  _molecularOrbitals = molecularOrbitals;
  _molecularOrbitals->addSensitiveObject(this->_self);
  for_spin(nOccupiedOrbitals, _aufbauOccupations) {
    _aufbauOccupations_spin.setZero();
    for (unsigned int i = 0; i < nOccupiedOrbitals_spin; ++i) {
      _aufbauOccupations_spin[i] = 1.0;
    }
  };
  if (update)
    _outOfDate = true;
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixController<SCFMode>::attachOrbitals(std::shared_ptr<OrbitalController<SCFMode>> molecularOrbitals,
                                                      const SpinPolarizedData<SCFMode, Eigen::VectorXd>& occupations,
                                                      bool update) {
  _molecularOrbitals = molecularOrbitals;
  _molecularOrbitals->addSensitiveObject(this->_self);
  _aufbauOccupations = occupations;
  if (update)
    _outOfDate = true;
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixController<SCFMode>::setDiskMode(bool diskmode, std::string fBaseName, std::string id) {
  if (diskmode) {
    _fBaseName = fBaseName;
    _id = id;
    assert(!_fBaseName.empty() && "Need to set file path when setting DensityMatrixController to disk mode.");
    assert(!_id.empty() && "Need to set file ID when setting DensityMatrixController to disk mode.");
  }
  if (diskmode and !_diskmode) {
    this->toHDF5(_fBaseName, _id);
    _densityMatrix.reset();
    _occupations.reset();
  }
  else if (!diskmode and _diskmode) {
    fromHDF5(_fBaseName, _id);
  }
  _diskmode = diskmode;
}

template<Options::SCF_MODES SCFMode>
DensityMatrix<SCFMode> DensityMatrixController<SCFMode>::getDensityMatrix() {
  if (_outOfDate) {
    bool tmp = _diskmode;
    this->setDiskMode(false, _fBaseName, _id);
    this->updateDensityMatrix();
    DensityMatrix<SCFMode> dmat(*_densityMatrix);
    this->setDiskMode(tmp, _fBaseName, _id);
    return dmat;
  }
  else if (_diskmode) {
    assert(!_fBaseName.empty());
    assert(!_id.empty());
    this->fromHDF5(_fBaseName, _id);
    DensityMatrix<SCFMode> dmat(*_densityMatrix);
    _densityMatrix.reset();
    _occupations.reset();
    return dmat;
  }
  else {
    return (DensityMatrix<SCFMode>)*_densityMatrix;
  }
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, Eigen::VectorXd> DensityMatrixController<SCFMode>::getOccupations(bool aufbau) {
  if (aufbau)
    return _aufbauOccupations;

  if (_outOfDate) {
    bool tmp = _diskmode;
    this->setDiskMode(false, _fBaseName, _id);
    this->updateDensityMatrix();
    SpinPolarizedData<SCFMode, Eigen::VectorXd> occ(*_occupations);
    this->setDiskMode(tmp, _fBaseName, _id);
    return occ;
  }
  else if (_diskmode) {
    assert(!_fBaseName.empty());
    assert(!_id.empty());
    this->fromHDF5(_fBaseName, _id);
    SpinPolarizedData<SCFMode, Eigen::VectorXd> occ(*_occupations);
    _densityMatrix.reset();
    _occupations.reset();
    return occ;
  }
  else {
    if (_occupations == nullptr)
      return _aufbauOccupations;
    SpinPolarizedData<SCFMode, Eigen::VectorXd> occ(*_occupations);
    return *_occupations;
  }
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixController<SCFMode>::setDensityMatrix(const DensityMatrix<SCFMode>& densityMatrix) {
  assert(densityMatrix.getBasisController() == _basisController);
  bool tmp = _diskmode;
  this->setDiskMode(false, _fBaseName, _id);
  _densityMatrix.reset(new DensityMatrix<SCFMode>(densityMatrix));
  _outOfDate = false;
  this->notifyObjects();
  this->setDiskMode(tmp, _fBaseName, _id);
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixController<SCFMode>::updateDensityMatrix() {
  assert(_molecularOrbitals);
  // Ensure that the aufbau occupations have the dimension of the basis.
  rebuildAufbauOccupations();
  SpinPolarizedData<SCFMode, Eigen::VectorXd> energies = _molecularOrbitals->getEigenvalues();
  _occupations.reset(new SpinPolarizedData<SCFMode, Eigen::VectorXd>(_aufbauOccupations));
  auto& occupations = *_occupations;
  if (_degeneracyThreshold != 0.0) {
    for_spin(energies, occupations) {
      for (unsigned int k = 0; k < occupations_spin.size(); ++k) {
        double sum = occupations_spin[k];
        unsigned int t = 0;
        for (t = k + 1; t < occupations_spin.size(); ++t) {
          if (fabs(energies_spin[k] - energies_spin[t]) < _degeneracyThreshold) {
            sum += occupations_spin[t];
          }
          else {
            break;
          }
        }
        sum /= (double)(t - k);
        for (; k < t; k++)
          occupations_spin[k] = sum;
        k -= 1;
      }
    };
  }
  CoefficientMatrix<SCFMode> coefficients = _molecularOrbitals->getCoefficients();
  DensityMatrix<SCFMode> dmat(_molecularOrbitals->getBasisController());
  if (_oneEIntController) {
    const double occ = (SCFMode == RESTRICTED) ? 2.0 : 1.0;
    Eigen::MatrixXd S = _oneEIntController->getOverlapIntegrals();
    for_spin(dmat, coefficients, occupations) {
      const unsigned int nOcc = occupations_spin.sum() / occ;
      Eigen::MatrixXd S_MO = coefficients_spin.leftCols(nOcc).transpose() * S * coefficients_spin.leftCols(nOcc);
      if (nOcc > 1) {
        // Jacobi SVD for pseudo inverse of overlapMO
        Eigen::MatrixXd singMatrix = Eigen::MatrixXd::Zero(nOcc, nOcc);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S_MO, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::VectorXd sing = svd.singularValues();
        // Invert singular values that are above a certain threshold and add to diagonal matrix
        double invThreshold = _basisController->getPrescreeningThreshold() * nOcc * sing.maxCoeff();
        for (unsigned i = 0; i < nOcc; ++i) {
          singMatrix(i, i) = fabs(sing(i)) >= invThreshold ? 1 / sing(i) : 0.0;
        }
        Eigen::MatrixXd inverseS_MO = svd.matrixU() * singMatrix * svd.matrixV().transpose();
        dmat_spin =
            occ * symmetrize(coefficients_spin.leftCols(nOcc) * inverseS_MO * coefficients_spin.leftCols(nOcc).transpose());
      }
      else if (nOcc == 1) {
        dmat_spin = occ * symmetrize(coefficients_spin.leftCols(nOcc) * 1 / S_MO(0, 0) *
                                     coefficients_spin.leftCols(nOcc).transpose());
      }
    };
  }
  else {
    for_spin(dmat, occupations, coefficients) {
      dmat_spin = symmetrize(coefficients_spin * occupations_spin.asDiagonal() * coefficients_spin.transpose());
    };
  }
  bool tmp = _diskmode;
  this->setDiskMode(false, _fBaseName, _id);
  _densityMatrix.reset(new DensityMatrix<SCFMode>(dmat));
  _outOfDate = false;
  this->setDiskMode(tmp, _fBaseName, _id);
}

template<>
void DensityMatrixController<Options::SCF_MODES::RESTRICTED>::toHDF5(std::string fBaseName, std::string id) {
  if (!_occupations or !_densityMatrix) {
    bool tmp = _diskmode;
    this->setDiskMode(false, _fBaseName, _id);
    this->getDensityMatrix();
    this->setDiskMode(tmp, _fBaseName, _id);
  }
  auto occupations = this->getOccupations();
  std::string name = fBaseName + ".dmat.res.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "densityMatrix", *_densityMatrix);
  HDF5::save(file, "occupations", occupations);
  HDF5::save_scalar_attribute(file, "ID", id);
  file.close();
}

template<>
void DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>::toHDF5(std::string fBaseName, std::string id) {
  if (!_occupations or !_densityMatrix) {
    bool tmp = _diskmode;
    this->setDiskMode(false, _fBaseName, _id);
    this->getDensityMatrix();
    this->setDiskMode(tmp, _fBaseName, _id);
  }
  auto occupations = this->getOccupations();
  std::string name = fBaseName + ".dmat.unres.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "densityMatrix_alpha", _densityMatrix->alpha);
  HDF5::save(file, "densityMatrix_beta", _densityMatrix->beta);
  HDF5::save(file, "occupations_alpha", occupations.alpha);
  HDF5::save(file, "occupations_beta", occupations.beta);
  HDF5::save_scalar_attribute(file, "ID", id);
  file.close();

  // This code generates a total alpha and beta density matrix, for unrestricted SCF calculation.
  // std::string name2 = fBaseName + ".dmat.res.h5";
  // HDF5::H5File file2(name2.c_str(), H5F_ACC_TRUNC);
  // HDF5::save(file2, "densityMatrix", _densityMatrix->total());
  // HDF5::save_scalar_attribute(file2, "ID", id);
  // file2.close();
}
template<>
void DensityMatrixController<Options::SCF_MODES::RESTRICTED>::fromHDF5(std::string fBaseName, std::string id) {
  _densityMatrix.reset(new DensityMatrix<RESTRICTED>(_basisController));
  _occupations.reset(new SpinPolarizedData<RESTRICTED, Eigen::VectorXd>(_basisController->getNBasisFunctions()));
  HDF5::Filepath name(fBaseName + ".dmat.res.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "densityMatrix");
  HDF5::dataset_exists(file, "occupations");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  HDF5::load(file, "densityMatrix", *_densityMatrix);
  HDF5::load(file, "occupations", *_occupations);
  *_densityMatrix = symmetrize(*_densityMatrix);
  unsigned int nBasisFunctions = _basisController->getNBasisFunctions();
  if (nBasisFunctions != _densityMatrix->cols() || nBasisFunctions != _densityMatrix->rows() ||
      nBasisFunctions != _occupations->size())
    throw SerenityError("ERROR: Failed loading a density matrix from disk. The dimensions are incorrect!");
  file.close();
}
template<>
void DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>::fromHDF5(std::string fBaseName, std::string id) {
  _densityMatrix.reset(new DensityMatrix<UNRESTRICTED>(_basisController));
  _occupations.reset(new SpinPolarizedData<UNRESTRICTED, Eigen::VectorXd>(_basisController->getNBasisFunctions()));
  HDF5::Filepath name(fBaseName + ".dmat.unres.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
  HDF5::dataset_exists(file, "densityMatrix_alpha");
  HDF5::dataset_exists(file, "densityMatrix_beta");
  HDF5::dataset_exists(file, "occupations_alpha");
  HDF5::dataset_exists(file, "occupations_beta");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  HDF5::load(file, "densityMatrix_alpha", _densityMatrix->alpha);
  HDF5::load(file, "densityMatrix_beta", _densityMatrix->beta);
  HDF5::load(file, "occupations_alpha", _occupations->alpha);
  HDF5::load(file, "occupations_beta", _occupations->beta);
  unsigned int nBasisFunctions = _basisController->getNBasisFunctions();
  if (nBasisFunctions != _densityMatrix->alpha.cols() || nBasisFunctions != _densityMatrix->alpha.rows() ||
      nBasisFunctions != _densityMatrix->beta.cols() || nBasisFunctions != _densityMatrix->beta.rows() ||
      nBasisFunctions != _occupations->alpha.size() || nBasisFunctions != _occupations->beta.size())
    throw SerenityError("ERROR: Failed loading a density matrix from disk. The dimensions are incorrect!");
  file.close();
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixController<SCFMode>::notify() {
  _outOfDate = true;
  _densityMatrix.reset();
  _occupations.reset();
  this->notifyObjects();
}

template<Options::SCF_MODES SCFMode>
void DensityMatrixController<SCFMode>::rebuildAufbauOccupations() {
  // Reset the aufbau occupations,since it is not clear
  // whether a notify was triggered by the basis controller
  //(via the OrbitalController)
  const double occ = (SCFMode == RESTRICTED) ? 2.0 : 1.0;
  const unsigned int nBasisFunctions = _basisController->getNBasisFunctions();
  for_spin(_aufbauOccupations) {
    const unsigned int nOcc = _aufbauOccupations_spin.sum() / occ;
    _aufbauOccupations_spin = Eigen::VectorXd::Zero(nBasisFunctions);
    _aufbauOccupations_spin.head(nOcc) = Eigen::VectorXd::Constant(nOcc, occ);
  };
}

template class DensityMatrixController<Options::SCF_MODES::RESTRICTED>;
template class DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
