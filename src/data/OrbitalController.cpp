/**
 * @file   OrbitalController.cpp
 *
 * @date   last rework Jan 16. 2017
 * @author Thomas Dresselhaus, Jan Unsleber
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
#include "data/OrbitalController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "data/SpinPolarizedData.h"
#include "data/matrices/CoefficientMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "data/matrices/MatrixInBasis.h"
#include "integrals/OneElectronIntegralController.h"
#include "io/FormattedOutput.h"
#include "io/HDF5.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "misc/Timing.h"
#include "misc/WarningTracker.h"

namespace Serenity {
using namespace std;

template<Options::SCF_MODES T>
OrbitalController<T>::OrbitalController(std::unique_ptr<CoefficientMatrix<T>> coefficients,
                                        std::shared_ptr<BasisController> basisController,
                                        std::unique_ptr<SpinPolarizedData<T, Eigen::VectorXd>> eigenvalues)
  : NotifyingClass<OrbitalController<T>>(),
    _coefficients(move(coefficients)),
    _basisController(basisController),
    _eigenvalues(move(eigenvalues)),
    _canOrthTh(1.0e-7),
    _linearDependent(false),
    _nZero(0) {
  auto& tmp = *_eigenvalues;
  for_spin(tmp) {
    if (!isDefinedInSameBasis(*_coefficients, *this))
      throw SerenityError("OrbitalController: Coefficients are not defined in the correct basis");
    if (tmp_spin.size() != _basisController->getNBasisFunctions())
      throw SerenityError("OrbitalController: The number of eigenvalues does not match the number of orbitals.");
  };
  _basisController->addSensitiveObject(this->_self);
}

template<Options::SCF_MODES T>
OrbitalController<T>::OrbitalController(std::shared_ptr<BasisController> basisController)
  : NotifyingClass<OrbitalController<T>>(),
    _coefficients(new CoefficientMatrix<T>(basisController)),
    _basisController(basisController),
    _eigenvalues(new SpinPolarizedData<T, Eigen::VectorXd>(basisController->getNBasisFunctions())),
    _canOrthTh(1.0e-7),
    _linearDependent(false),
    _nZero(0) {
  assert(_basisController);
  _basisController->addSensitiveObject(this->_self);
}

template<Options::SCF_MODES T>
OrbitalController<T>::OrbitalController(const OrbitalController<T>& orig)
  : NotifyingClass<OrbitalController<T>>(),
    _coefficients(new CoefficientMatrix<T>(*orig._coefficients)),
    _basisController(orig._basisController),
    _eigenvalues(new SpinPolarizedData<T, Eigen::VectorXd>(*orig._eigenvalues)),
    _canOrthTh(orig._canOrthTh),
    _nZero(0) {
  assert(isDefinedInSameBasis(orig, *this));
  _basisController->addSensitiveObject(this->_self);
}

template<Options::SCF_MODES T>
OrbitalController<T>::OrbitalController(std::string filePath, std::shared_ptr<BasisController> basisController, std::string id)
  : _coefficients(new CoefficientMatrix<T>(basisController)),
    _basisController(basisController),
    _eigenvalues(new SpinPolarizedData<T, Eigen::VectorXd>(basisController->getNBasisFunctions())),
    _canOrthTh(1.0e-7),
    _linearDependent(false),
    _nZero(0),
    _fBaseName(filePath),
    _id(id) {
  fromHDF5(filePath, id);
  _basisController->addSensitiveObject(this->_self);
}
template<Options::SCF_MODES T>
OrbitalController<T>::~OrbitalController() {
}

template<Options::SCF_MODES T>
unsigned int OrbitalController<T>::getNOrbitals() const {
  return _basisController->getNBasisFunctions();
}

template<Options::SCF_MODES T>
void OrbitalController<T>::calculateTransformationX(std::shared_ptr<OneElectronIntegralController> oneIntController) {
  if (_X.cols() > 0)
    return;
  // Calculate canonical orthogonalization matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(oneIntController->getOverlapIntegrals());
  auto U = es.eigenvectors();
  auto s = es.eigenvalues();
  int n = s.rows();
  _nZero = 0;
  Eigen::VectorXd sigma(n);
  Eigen::VectorXd sigmaInvers(n);
  sigma.setZero();
  sigmaInvers.setZero();
  for (int i = n - 1; i >= 0; --i) {
    if (s(i) > _canOrthTh) {
      sigma[i] = s[i];
      sigmaInvers[i] = 1 / s[i];
      ++_nZero;
    }
    else {
      i = 0;
    }
  }
  _X = Eigen::MatrixXd(U * sigmaInvers.array().sqrt().matrix().asDiagonal()).eval();
  _X = _X.rightCols(_nZero).eval();
  _Xinv = Eigen::MatrixXd(U * sigma.array().sqrt().matrix().asDiagonal()).eval();
  _Xinv = _Xinv.rightCols(_nZero).eval();
  _nZero = n - _nZero;
  _linearDependent = _nZero > 0;
  if (_linearDependent) {
    WarningTracker::printWarning((std::string) "Warning: Basis-Set (near) linear dependent. Will try to use canonical "
                                               "orthogonalization. Removed " +
                                     _nZero + " columns from transformation matrix.",
                                 true);
  }
}

template<Options::SCF_MODES T>
void OrbitalController<T>::updateOrbitals(const FockMatrix<T>& fockMatrix,
                                          std::shared_ptr<OneElectronIntegralController> oneIntController) {
  auto shift = std::make_pair<Eigen::VectorXd, SpinPolarizedData<T, Eigen::VectorXd>>(
      Eigen::VectorXd(Eigen::VectorXd::Zero(2)), SpinPolarizedData<T, Eigen::VectorXd>(Eigen::VectorXd(0)));
  this->updateOrbitals(shift, fockMatrix, oneIntController);
}

template<Options::SCF_MODES T>
void OrbitalController<T>::updateOrbitals(const std::pair<Eigen::VectorXd, SpinPolarizedData<T, Eigen::VectorXd>> levelshift,
                                          const FockMatrix<T>& fockMatrix,
                                          std::shared_ptr<OneElectronIntegralController> oneIntController) {
  Timings::takeTime("Tech. -    Fock Matrix Solving");

  auto& occupation = levelshift.second;

  // Calculate transformation matrix
  if (!_fIsInOthoBasis and !_customS)
    calculateTransformationX(oneIntController);

  auto eps(std::move(this->getEigenvalues()));
  auto c(std::move(this->getCoefficients()));

  if (_customS) {
    const auto& S = *_customS;
    for_spin(fockMatrix, eps, c, S) {
      Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(fockMatrix_spin, S_spin);
      c_spin = es.eigenvectors();
      eps_spin = es.eigenvalues();
    };
  }
  else if (_fIsInOthoBasis) {
    for_spin(fockMatrix, eps, c) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(fockMatrix_spin);
      c_spin = es.eigenvectors();
      eps_spin = es.eigenvalues();
    };
  }
  else {
    for_spin(fockMatrix, eps, c, occupation) {
      // transform F into orthonormal basis
      Eigen::MatrixXd tmpFockMatrix(_X.transpose() * fockMatrix_spin * _X);

      // Level-shift
      if (levelshift.first[0] > 0.0 and !_firstIteration) {
        auto C_ortho = _Xinv.transpose() * c_spin;
        // transform into old MO space
        tmpFockMatrix = C_ortho.transpose() * tmpFockMatrix * C_ortho;
        const double diagonalShift = levelshift.first[0];
        for (unsigned int i = 0; i < occupation_spin.size(); i++) {
          if (occupation_spin[i] < 1e-9) {
            tmpFockMatrix(i, i) += diagonalShift;
            for (unsigned int j = 0; j < occupation_spin.size(); j++) {
              if (occupation_spin[j] > 1e-9) {
                tmpFockMatrix(i, j) *= levelshift.first[1];
                tmpFockMatrix(j, i) *= levelshift.first[1];
              }
            }
          }
        }
        // transform back from old MO space
        tmpFockMatrix = C_ortho * tmpFockMatrix * C_ortho.transpose();
      }

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(tmpFockMatrix);
      c_spin.setZero();
      Eigen::MatrixXd newCoeffs = _X * es.eigenvectors();
      for (unsigned int i = 0; i < newCoeffs.cols(); ++i) {
        if (newCoeffs(0, i) < 0.0)
          newCoeffs.col(i) *= -1.0;
      }
      c_spin.leftCols(newCoeffs.cols()) = newCoeffs;
      eps_spin = Eigen::VectorXd::Constant(_basisController->getNBasisFunctions(), numeric_limits<double>::infinity());
      Eigen::VectorXd newEigenvalues = es.eigenvalues();
      eps_spin.segment(0, newEigenvalues.size()) = newEigenvalues;
    };
  }
  _firstIteration = false;
  this->updateOrbitals(c, eps);
  Timings::timeTaken("Tech. -    Fock Matrix Solving");
}

template<Options::SCF_MODES T>
void OrbitalController<T>::setDiskMode(bool diskmode, std::string fBaseName, std::string id) {
  if (diskmode) {
    _fBaseName = fBaseName;
    _id = id;
    assert(!_fBaseName.empty() && "Need to set file path when setting OrbitalController to disk mode.");
    assert(!_id.empty() && "Need to set file ID when setting OrbitalController to disk mode.");
  }
  if (diskmode and _keepInMemory) {
    this->toHDF5(_fBaseName, _id);
    _coefficients.reset();
    _eigenvalues.reset();
  }
  else if (!diskmode and !_keepInMemory) {
    fromHDF5(_fBaseName, _id);
  }
  _keepInMemory = !diskmode;
}

/// @returns the coefficients determining the orbitals in connection with the basis.
template<Options::SCF_MODES T>
CoefficientMatrix<T> OrbitalController<T>::getCoefficients() {
  if (!_keepInMemory && (_coefficients == nullptr)) {
    assert(!_fBaseName.empty());
    assert(!_id.empty());
    this->coefficientsfromHDF5(_fBaseName, _id);
    CoefficientMatrix<T> coefficients(*_coefficients);
    _coefficients.reset();
    return coefficients;
  }
  else {
    return (CoefficientMatrix<T>)*_coefficients;
  }
}

/// @returns the orbital energies.
template<Options::SCF_MODES T>
SpinPolarizedData<T, Eigen::VectorXd> OrbitalController<T>::getEigenvalues() {
  if (!_keepInMemory && (_eigenvalues == nullptr)) {
    assert(!_fBaseName.empty());
    assert(!_id.empty());
    this->eigenvaluesfromHDF5(_fBaseName, _id);
    SpinPolarizedData<T, Eigen::VectorXd> eigenvalues(*_eigenvalues);
    _eigenvalues.reset();
    return eigenvalues;
  }
  else {
    return (SpinPolarizedData<T, Eigen::VectorXd>)*_eigenvalues;
  }
}

template<Options::SCF_MODES T>
void OrbitalController<T>::updateOrbitals(const CoefficientMatrix<T>& updatedCoefficients,
                                          const SpinPolarizedData<T, Eigen::VectorXd>& updatedEigenvalues) {
  assert(updatedCoefficients.getBasisController() == _basisController);
  bool keeptmp = _keepInMemory;
  if (!_keepInMemory) {
    assert(!_fBaseName.empty());
    assert(!_id.empty());
    _eigenvalues.reset(new SpinPolarizedData<T, Eigen::VectorXd>(_basisController->getNBasisFunctions()));
    _coefficients.reset(new CoefficientMatrix<T>(_basisController));
    _keepInMemory = true;
  }
  *_eigenvalues = updatedEigenvalues;
  *_coefficients = updatedCoefficients;
  this->notifyObjects();

  _keepInMemory = keeptmp;
  if (!_keepInMemory) {
    this->toHDF5(_fBaseName, _id);
    _eigenvalues.reset();
  }
}

template<Options::SCF_MODES T>
void OrbitalController<T>::fromHDF5(std::string fBaseName, std::string id) {
  this->coefficientsfromHDF5(fBaseName, id);
  this->eigenvaluesfromHDF5(fBaseName, id);
  _firstIteration = false;
}

template<>
void OrbitalController<Options::SCF_MODES::RESTRICTED>::coefficientsfromHDF5(std::string fBaseName, std::string id) {
  HDF5::Filepath name(fBaseName + ".orbs.res.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5::dataset_exists(file, "coefficients");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  _coefficients.reset(new CoefficientMatrix<Options::SCF_MODES::RESTRICTED>(_basisController));
  HDF5::load(file, "coefficients", *_coefficients);
  file.close();
}
template<>
void OrbitalController<Options::SCF_MODES::UNRESTRICTED>::coefficientsfromHDF5(std::string fBaseName, std::string id) {
  HDF5::Filepath name(fBaseName + ".orbs.unres.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5::dataset_exists(file, "coefficients_alpha");
  HDF5::dataset_exists(file, "coefficients_beta");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  _coefficients.reset(new CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED>(_basisController));
  HDF5::load(file, "coefficients_alpha", _coefficients->alpha);
  HDF5::load(file, "coefficients_beta", _coefficients->beta);
  file.close();
}

template<>
void OrbitalController<Options::SCF_MODES::RESTRICTED>::toHDF5(std::string fBaseName, std::string id) {
  std::string name = fBaseName + ".orbs.res.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "eigenvalues", *_eigenvalues);
  HDF5::save(file, "coefficients", *_coefficients);
  HDF5::save_scalar_attribute(file, "ID", id);
  file.close();
}
template<>
void OrbitalController<Options::SCF_MODES::UNRESTRICTED>::toHDF5(std::string fBaseName, std::string id) {
  std::string name = fBaseName + ".orbs.unres.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "eigenvalues_alpha", _eigenvalues->alpha);
  HDF5::save(file, "eigenvalues_beta", _eigenvalues->beta);
  HDF5::save(file, "coefficients_alpha", _coefficients->alpha);
  HDF5::save(file, "coefficients_beta", _coefficients->beta);
  HDF5::save_scalar_attribute(file, "ID", id);
  file.close();
}
template<>
void OrbitalController<Options::SCF_MODES::RESTRICTED>::eigenvaluesfromHDF5(std::string fBaseName, std::string id) {
  HDF5::Filepath name(fBaseName + ".orbs.res.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5::dataset_exists(file, "eigenvalues");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  _eigenvalues.reset(
      new SpinPolarizedData<Options::SCF_MODES::RESTRICTED, Eigen::VectorXd>(_basisController->getNBasisFunctions()));
  HDF5::load(file, "eigenvalues", *_eigenvalues);
  file.close();
}
template<>
void OrbitalController<Options::SCF_MODES::UNRESTRICTED>::eigenvaluesfromHDF5(std::string fBaseName, std::string id) {
  HDF5::Filepath name(fBaseName + ".orbs.unres.h5");
  HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5::dataset_exists(file, "eigenvalues_alpha");
  HDF5::dataset_exists(file, "eigenvalues_beta");
  HDF5::attribute_exists(file, "ID");
  HDF5::check_attribute(file, "ID", id);
  _eigenvalues.reset(
      new SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::VectorXd>(_basisController->getNBasisFunctions()));
  HDF5::load(file, "eigenvalues_alpha", _eigenvalues->alpha);
  HDF5::load(file, "eigenvalues_beta", _eigenvalues->beta);
  file.close();
}

template class OrbitalController<Options::SCF_MODES::RESTRICTED>;
template class OrbitalController<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
