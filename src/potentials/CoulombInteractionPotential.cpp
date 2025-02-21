/**
 * @file CoulombInteractionPotential.cpp
 * @author: Kevin Klahr
 *
 * @date 29. November 2016
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
#include "potentials/CoulombInteractionPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"            //Atom to basis function mapping for gradients.
#include "basis/AtomCenteredBasisControllerFactory.h"     //Construction of joined supersystem basis in gradients
#include "basis/BasisController.h"                        //Basis controller definition.
#include "basis/BasisFunctionMapper.h"                    //Construct joined fitting basis.
#include "data/ElectronicStructure.h"                     //Access to density matrices for gradients.
#include "data/matrices/DensityMatrixController.h"        //Access to environment density matrices.
#include "data/matrices/FockMatrix.h"                     //Construct joined fitting basis.
#include "geometry/Geometry.h"                            //Needed for gradients.
#include "integrals/RI_J_IntegralControllerFactory.h"     //Used as dummy for metric inversion/llT
#include "integrals/looper/CoulombInteractionIntLooper.h" //Loop integrals.
#include "integrals/looper/TwoElecThreeCenterIntLooper.h" //Loop integrals.
#include "integrals/wrappers/Libint.h"                    //getSharedPtr()
#include "io/FormattedOutputStream.h"                     //Filtered output streams.
#include "io/HDF5.h"                                      //HDF5 IO
#include "misc/Timing.h"                                  //Timings.
#include "potentials/RICoulombPotential.h"                //Gradient construction. Substract intra-subsystem part.
#include "settings/Settings.h"                            //Settings definition (gradients only).
#include "system/SystemController.h"                      //Settings and getElectronicStructure
/* Include Std and External Headers */
#include <sys/stat.h> //Check for an already existing file on disk.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
CoulombInteractionPotential<SCFMode>::CoulombInteractionPotential(
    std::shared_ptr<SystemController> actSystem, std::vector<std::shared_ptr<SystemController>> envSystems,
    const std::shared_ptr<BasisController> actBasis,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensityMatrixController,
    const std::shared_ptr<BasisController> actAuxBasis, std::vector<std::shared_ptr<BasisController>> envAuxBasis, bool isPassive)
  : Potential<SCFMode>(actBasis),
    _actSystem(actSystem),
    _envDMatController(envDensityMatrixController),
    _actAuxBasis(actAuxBasis),
    _envAuxBasis(envAuxBasis),
    _isPassive(isPassive) {
  Timings::takeTime("FDE -        Coulomb Pot.");
  for (auto e : envSystems)
    _envSystems.push_back(e);
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_actAuxBasis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  for (auto& mat : envDensityMatrixController) {
    mat->getDensityMatrix().getBasisController()->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  }
  for (auto& bas : envAuxBasis) {
    bas->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  }
  _mode = Options::DENS_FITS::RI;
  _fBaseName = actSystem->getSystemPath() + actSystem->getSystemName();
  Timings::timeTaken("FDE -        Coulomb Pot.");
};

template<Options::SCF_MODES SCFMode>
CoulombInteractionPotential<SCFMode>::CoulombInteractionPotential(
    std::shared_ptr<SystemController> actSystem, std::vector<std::shared_ptr<SystemController>> envSystems,
    const std::shared_ptr<BasisController> actBasis,
    std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> envDensityMatrixController, bool isPassive)
  : Potential<SCFMode>(actBasis),
    _actSystem(actSystem),
    _envDMatController(envDensityMatrixController),
    _actAuxBasis(nullptr),
    _envAuxBasis(0, nullptr),
    _isPassive(isPassive) {
  for (auto e : envSystems)
    _envSystems.push_back(e);
  Timings::takeTime("FDE -        Coulomb Pot.");
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  for (auto& mat : envDensityMatrixController) {
    mat->getDensityMatrix().getBasisController()->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  }
  _mode = Options::DENS_FITS::NONE;
  _fBaseName = actSystem->getSystemPath() + actSystem->getSystemName();
  Timings::timeTaken("FDE -        Coulomb Pot.");
};

template<Options::SCF_MODES SCFMode>
void CoulombInteractionPotential<SCFMode>::getFromDisk() {
  // Check if the file already exists.
  if (not fromHDF5()) {
    calculateFockMatrix();
    OutputControl::dOut << "Writing passive Coulomb contributions to file." << std::endl;
    toHDF5();
  }
}
template<Options::SCF_MODES SCFMode>
bool CoulombInteractionPotential<SCFMode>::fromHDF5() {
  std::string fileName = _fBaseName + ".pasCoulomb.h5";
  struct stat buffer;
  auto actSystem = _actSystem.lock();
  if (stat(fileName.c_str(), &buffer) == 0) {
    HDF5::Filepath name(fileName);
    HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
    try {
      HDF5::attribute_exists(file, "ID");
      HDF5::dataset_exists(file, "passiveCoulombContribution");
      // Note: This lines will enforce a recalculation for old/wrong id
      // files. If you want to reuse old Fock matrices. You could try to remove
      // this check.
      HDF5::check_attribute(file, "ID", actSystem->getSystemIdentifier());
    }
    catch (...) {
      return false;
    }
    OutputControl::dOut << "Loading passive Coulomb contributions from file." << std::endl;
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    FockMatrix<RESTRICTED> f(this->_basis);
    HDF5::load(file, "passiveCoulombContribution", f);
    auto& pot = *_potential;
    for_spin(pot) {
      pot_spin = f;
    };
    file.close();
  }
  else {
    return false;
  }
  return true;
}
template<>
void CoulombInteractionPotential<RESTRICTED>::toHDF5() {
  auto actSystem = _actSystem.lock();
  std::string fileName = _fBaseName + ".pasCoulomb.h5";
  HDF5::H5File file(fileName.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "passiveCoulombContribution", *_potential);
  HDF5::save_scalar_attribute(file, "ID", actSystem->getSystemIdentifier());
  file.close();
}
template<>
void CoulombInteractionPotential<UNRESTRICTED>::toHDF5() {
  auto actSystem = _actSystem.lock();
  std::string fileName = _fBaseName + ".pasCoulomb.h5";
  HDF5::H5File file(fileName.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "passiveCoulombContribution", _potential->alpha);
  HDF5::save_scalar_attribute(file, "ID", actSystem->getSystemIdentifier());
  file.close();
}

template<Options::SCF_MODES SCFMode>
double CoulombInteractionPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (!_potential)
    this->getMatrix();
  Timings::takeTime("FDE -        Coulomb Pot.");
  auto& pot = *_potential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("FDE -        Coulomb Pot.");
  return energy;
};
template<Options::SCF_MODES SCFMode>
void CoulombInteractionPotential<SCFMode>::calculateFockMatrix() {
  _potential.reset(new FockMatrix<SCFMode>(this->_basis));
  auto& F = *_potential;
  for_spin(F) {
    F_spin.setZero();
  };

  unsigned nb_act = this->_basis->getNBasisFunctions();

  BasisFunctionMapper basisMapper(this->_basis);
  BasisFunctionMapper auxBasisMapper(_actAuxBasis);
  unsigned int nThreads = omp_get_max_threads();

  double screening = _actSystem.lock()->getSettings().basis.integralThreshold;
  double screeningAct = (screening != 0) ? screening : this->_basis->getPrescreeningThreshold();

  if (_mode == Options::DENS_FITS::RI) {
    for (unsigned int i = 0; i < _envDMatController.size(); ++i) {
      std::vector<MatrixInBasis<RESTRICTED>> fock_threads(nThreads, MatrixInBasis<RESTRICTED>(this->_basis));

      auto envBasis = _envDMatController[i]->getDensityMatrix().getBasisController();
      auto envAuxBasis = _envAuxBasis[i];

      MatrixInBasis<RESTRICTED> envDensityMatrix = _envDMatController[i]->getDensityMatrix().total();
      auto densptr = envDensityMatrix.data();
      unsigned nb_env = envDensityMatrix.rows();

      auto maxDens = envDensityMatrix.shellWiseAbsMax();
      auto maxDensPtr = maxDens.data();
      unsigned ns = maxDens.rows();

      // Environment prescreening threshold.
      double screeningEnv = (screening != 0) ? screening : envBasis->getPrescreeningThreshold();

      // Supersystem basis controller.
      auto superSystemBasis = basisMapper.getCombinedBasis(envBasis);
      auto superSystemAuxBasis = auxBasisMapper.getCombinedBasis(envAuxBasis);
      unsigned nx = superSystemAuxBasis->getNBasisFunctions();
      Eigen::MatrixXd sumMat = Eigen::MatrixXd::Zero(nx, nThreads);
      auto sumptr = sumMat.data();

      auto ri_j_IntController = RI_J_IntegralControllerFactory::getInstance().produce(superSystemBasis, superSystemAuxBasis);

      // Distribute function for first contraction.
      auto distribute1 = [&](unsigned i, unsigned j, unsigned K, double integral, unsigned threadId) {
        double perm = (i == j ? 1.0 : 2.0);
        sumptr[nx * threadId + K] += perm * integral * densptr[i * nb_env + j];
      };

      // Prescreening function for first contraction.
      auto prescreen1 = [&](unsigned i, unsigned j, unsigned, double schwarz) {
        unsigned long ij = i * ns + j;
        return (maxDensPtr[ij] * schwarz < screeningEnv);
      };

      // Perform first contraction.
      TwoElecThreeCenterIntLooper looper1(LIBINT_OPERATOR::coulomb, 0, envBasis, superSystemAuxBasis, screeningEnv);
      looper1.loopNoDerivative(distribute1, prescreen1);

      // Solve linear system.
      Eigen::VectorXd sumPMuNu_DMuNu = sumMat.rowwise().sum();
      sumMat.resize(1, 1);
      Eigen::VectorXd coefficients = ri_j_IntController->getLLTMetric().solve(sumPMuNu_DMuNu).eval();

      // Distribute function for second contraction.
      auto coeff = coefficients.data();
      auto coeffMax = superSystemAuxBasis->shellWiseAbsMax(coefficients);
      auto distribute2 = [&](unsigned i, unsigned j, unsigned K, double integral, unsigned threadId) {
        // Only perform half the needed contractions and symmetrize afterwards.
        fock_threads[threadId].data()[i * nb_act + j] += integral * coeff[K];
      };

      // Prescreening function for second contraction.
      auto prescreen2 = [&](unsigned, unsigned, unsigned K, double schwarz) {
        return (coeffMax(K) * schwarz < screening);
      };

      // Perform second contraction.
      TwoElecThreeCenterIntLooper looper2(LIBINT_OPERATOR::coulomb, 0, this->_basis, superSystemAuxBasis, screeningAct);
      looper2.loopNoDerivative(distribute2, prescreen2);

      // Add to Fock matrix.
      Eigen::Ref<Eigen::MatrixXd> Fc = fock_threads[0];
      for (unsigned i = 1; i < nThreads; i++) {
        Fc += fock_threads[i];
      }
      Eigen::MatrixXd Fc_sym = Fc + Fc.transpose();
      // The diagonal was fine, so it needs to be halved.
      Fc_sym.diagonal() *= 0.5;
      for_spin(F) {
        F_spin += Fc_sym;
      };
    }
  }
  else {
    for (unsigned int i = 0; i < _envDMatController.size(); ++i) {
      std::vector<MatrixInBasis<RESTRICTED>> eriContr(nThreads, MatrixInBasis<RESTRICTED>(this->_basis));
      auto envMat = _envDMatController[i]->getDensityMatrix().total();
      double screeningEnv = (screening != 0) ? screening : envMat.getBasisController()->getPrescreeningThreshold();
      double minPrescreeningActEnv = std::min(screeningAct, screeningEnv);

      CoulombInteractionIntLooper coulLooper(LIBINT_OPERATOR::coulomb, 0, this->_basis, envMat.getBasisController(),
                                             minPrescreeningActEnv);
      const auto maxDens = envMat.shellWiseAbsMax();

      auto prescreen = [&](const unsigned, const unsigned, unsigned int a, unsigned int b, const double schwartz) {
        return (std::abs(maxDens(a, b) * schwartz) < minPrescreeningActEnv);
      };
      auto const coulLooperFunction = [&eriContr, &envMat](const unsigned int& i, const unsigned int& j,
                                                           const unsigned int& a, const unsigned int& b,
                                                           double intValues, const unsigned int& threadID) {
        double perm = 1.0;
        if (a != b) {
          perm *= 2.0;
        }

        const double coul = perm * envMat(a, b) * intValues;
        eriContr[threadID](i, j) += coul;
        if (i != j) {
          eriContr[threadID](j, i) += coul;
        }
      };
      coulLooper.loopNoDerivative(coulLooperFunction, prescreen);

      // sum over all threads
      for (const auto& eriCont : eriContr) {
        for_spin(F) {
          F_spin += eriCont;
        };
      }
    }
  }
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& CoulombInteractionPotential<SCFMode>::getMatrix() {
  Timings::takeTime("FDE -        Coulomb Pot.");
  if (_isPassive and !_potential)
    getFromDisk();
  if (!_potential)
    calculateFockMatrix();
  Timings::timeTaken("FDE -        Coulomb Pot.");
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd CoulombInteractionPotential<SCFMode>::getGeomGradients() {
  auto actSystem = _actSystem.lock();
  auto atomsAct = actSystem->getAtoms();
  unsigned int nAtomsAct = atomsAct.size();
  Matrix<double> activeSystemGradientContr(atomsAct.size(), 3);
  activeSystemGradientContr.setZero();
  DensityMatrix<RESTRICTED> densityMatrix(actSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix().total());

  for (unsigned int i = 0; i < _envSystems.size(); ++i) {
    // Creating supersystem Geometry and Basis, necessary still for inverseM
    auto envSystem = _envSystems[i].lock();
    auto atomsEnv = envSystem->getAtoms();
    std::vector<std::shared_ptr<Atom>> ActiveFrozenAtoms;
    for (auto& atom : atomsAct) {
      ActiveFrozenAtoms.push_back(atom);
    }
    for (auto& atom : atomsEnv) {
      ActiveFrozenAtoms.push_back(atom);
    }
    unsigned int nAtoms = ActiveFrozenAtoms.size();
    Matrix<double> eriContr(nAtoms, 3);
    eriContr.setZero();

    auto superSystemGeometry = std::make_shared<Geometry>(ActiveFrozenAtoms);
    auto superSystemBasis = AtomCenteredBasisControllerFactory::produce(
        superSystemGeometry, actSystem->getSettings().basis.basisLibPath, actSystem->getSettings().basis.makeSphericalBasis,
        true, actSystem->getSettings().basis.firstECP, actSystem->getSettings().basis.label);
    auto superSystemAuxBasis = AtomCenteredBasisControllerFactory::produce(
        superSystemGeometry, actSystem->getSettings().basis.basisLibPath,
        actSystem->getSettings().basis.makeSphericalBasis, false, 999999999, actSystem->getSettings().basis.auxJLabel);

    // Getting necessary information from the active and environment system
    // Active
    auto basisController = actSystem->getAtomCenteredBasisController();
    auto auxBasisController = actSystem->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
    auto& auxBasisFunctions = auxBasisController->getBasis();
    auto basisIndicesRed = actSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
    auto mapping = basisController->getAtomIndicesOfBasisShells();
    auto nAuxBasFunc = auxBasisController->getNBasisFunctions();
    auto auxBasisIndicesRed = auxBasisController->getBasisIndicesRed();
    auto auxMapping = auxBasisController->getAtomIndicesOfBasisShells();

    // Environment
    auto basisControllerEnv = envSystem->getAtomCenteredBasisController();
    auto auxBasisControllerEnv = envSystem->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
    auto& auxBasisFunctionsEnv = auxBasisControllerEnv->getBasis();
    auto basisIndicesRedEnv = envSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
    auto mappingEnv = basisControllerEnv->getAtomIndicesOfBasisShells();
    auto nAuxBasFuncEnv = auxBasisControllerEnv->getNBasisFunctions();
    auto auxBasisIndicesRedEnv = auxBasisControllerEnv->getBasisIndicesRed();
    auto auxMappingEnv = auxBasisControllerEnv->getAtomIndicesOfBasisShells();
    DensityMatrix<RESTRICTED> densityMatrixEnv =
        envSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix().total();

    auto ri_j_IntController = RI_J_IntegralControllerFactory::getInstance().produce(superSystemBasis, superSystemAuxBasis);
    /*
     * Constructing the different parts of eq. ([1].[6]).
     *
     * Starting with the C vectors. The first step combines the
     * Density Matrix elements ij with the respecting 3c Integral <J|ij>.
     * This is technically a Matrix-Vector multiplication.
     * Has to be done 4 times in order to get all combinations of
     * active and environment systems basis and auxBasis.
     *
     */
#ifdef _OPENMP
    // create a vector of matrices for each thread
    std::vector<Eigen::VectorXd> CvecActPriv(omp_get_max_threads(), Eigen::VectorXd::Zero(nAuxBasFunc));
    std::vector<Eigen::VectorXd> CvecEnvPriv(omp_get_max_threads(), Eigen::VectorXd::Zero(nAuxBasFuncEnv));
#else
    // or just one
    std::vector<Eigen::VectorXd> CvecActPriv(1, Eigen::VectorXd::Zero(nAuxBasFunc));
    std::vector<Eigen::VectorXd> CvecEnvPriv(1, Eigen::VectorXd::Zero(nAuxBasFuncEnv));
#endif
    double screening = _actSystem.lock()->getSettings().basis.integralThreshold;
    double screeningEnv = (screening != 0) ? screening : basisControllerEnv->getPrescreeningThreshold();
    TwoElecThreeCenterIntLooper looper1(LIBINT_OPERATOR::coulomb, 0, basisControllerEnv, auxBasisController, screeningEnv);
    const auto maxDens = densityMatrixEnv.shellWiseAbsMax();
    auto const calcCVector1 = [&CvecActPriv, &densityMatrixEnv](const unsigned int& i, const unsigned int& j,
                                                                const unsigned int& K, Eigen::VectorXd& intValues,
                                                                const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecActPriv[threadId][K] += perm * intValues(0) * densityMatrixEnv(i, j);
    };
    looper1.loop(calcCVector1, maxDens.maxCoeff());
    double screeningAct = (screening != 0) ? screening : basisController->getPrescreeningThreshold();
    TwoElecThreeCenterIntLooper looper2(LIBINT_OPERATOR::coulomb, 0, basisController, auxBasisController, screeningAct);
    auto const calcCVector2 = [&CvecActPriv, &densityMatrix](const unsigned int& i, const unsigned int& j,
                                                             const unsigned int& K, Eigen::VectorXd& intValues,
                                                             const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecActPriv[threadId][K] += perm * intValues(0) * densityMatrix(i, j);
    };
    looper2.loop(calcCVector2, maxDens.maxCoeff());

    TwoElecThreeCenterIntLooper looper3(LIBINT_OPERATOR::coulomb, 0, basisController, auxBasisControllerEnv, screeningAct);

    auto const calcCVector3 = [&CvecEnvPriv, &densityMatrix](const unsigned int& i, const unsigned int& j,
                                                             const unsigned int& K, Eigen::VectorXd& intValues,
                                                             const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecEnvPriv[threadId][K] += perm * intValues(0) * densityMatrix(i, j);
    };
    looper3.loop(calcCVector3, maxDens.maxCoeff());

    TwoElecThreeCenterIntLooper looper4(LIBINT_OPERATOR::coulomb, 0, basisControllerEnv, auxBasisControllerEnv, screeningEnv);
    auto const calcCVector4 = [&CvecEnvPriv, &densityMatrixEnv](const unsigned int& i, const unsigned int& j,
                                                                const unsigned int& K, Eigen::VectorXd& intValues,
                                                                const unsigned int threadId) {
      const double perm = (i == j ? 1.0 : 2.0);
      CvecEnvPriv[threadId][K] += perm * intValues(0) * densityMatrixEnv(i, j);
    };
    looper4.loop(calcCVector4, maxDens.maxCoeff());
    /*
     * Has to be merged to a supersystem C^J because calculating the inverse of the
     * 2c-integral matrix in parts gives nonsensical results.
     */
    Eigen::VectorXd CvecAct(nAuxBasFunc);
    CvecAct.setZero();
    Eigen::VectorXd CvecEnv(nAuxBasFuncEnv);
    CvecEnv.setZero();
#ifdef _OPENMP
    for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
      CvecAct += CvecActPriv[i];
      CvecEnv += CvecEnvPriv[i];
    }
#else
    CvecAct += CvecActPriv[0];
    CvecEnv += CvecEnvPriv[0];
#endif
    Eigen::VectorXd Cvec(CvecAct.rows() + CvecEnv.rows());
    Cvec << CvecAct, CvecEnv;
    Eigen::VectorXd Dvec = ri_j_IntController->getLLTMetric().solve(Cvec);

    /*
     * The resulting Vector Dvec is the equivalent of the C^J vector in the paper.
     * The following set of loops calculates the complete 3rd part of eq ([1].[6])
     * (easily recognizable by the sign).
     *
     * Check RICoulombPotential::getGeomGradients()!
     *
     * Has to be done twice in order to get the combinations
     * auxBasisEnv/Basis and auxBasis/BasisEnv.
     *
     */
    Libint& libint = Libint::getInstance();
    libint.initialize(LIBINT_OPERATOR::coulomb, 1, 2);
#ifdef _OPENMP
    // create a vector of matrices for each thread
    std::vector<Eigen::MatrixXd> eriContrPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
    // or just one
    std::vector<Eigen::MatrixXd> eriContrPriv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

    for (unsigned int auxJ = 0; auxJ < auxBasisFunctions.size(); auxJ++) {
      const unsigned int nAuxJ = auxBasisFunctions[auxJ]->getNContracted();
      unsigned int offauxJ = auxBasisController->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
      for (unsigned int auxK = 0; auxK < auxBasisFunctionsEnv.size(); auxK++) {
        Eigen::MatrixXd intDerivs;
#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int nAuxK = auxBasisFunctionsEnv[auxK]->getNContracted();
        unsigned int offauxK = auxBasisControllerEnv->extendedIndex(auxK);
        if (libint.compute(LIBINT_OPERATOR::coulomb, 1, *auxBasisFunctions[auxJ], *auxBasisFunctionsEnv[auxK], intDerivs)) {
          double perm = 1.0;
          Eigen::MatrixXd prefac =
              0.5 * perm * Dvec.segment(offauxJ, nAuxJ) * Dvec.segment(offauxK + nAuxBasFunc, nAuxK).transpose();
          for (unsigned int iAtom = 0; iAtom < 2; ++iAtom) {
            unsigned int nAtom;
            switch (iAtom) {
              case (0):
                nAtom = auxMapping[auxJ];
                break;
              case (1):
                nAtom = auxMappingEnv[auxK] + atomsAct.size();
                break;
            }
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom * 3 + iDirection).data(), nAuxK, nAuxJ);
              eriContrPriv[threadId](nAtom, iDirection) -= tmp.transpose().cwiseProduct(prefac).sum();
            }
          }
        }
      }
    }

    for (unsigned int auxJ = 0; auxJ < auxBasisFunctionsEnv.size(); auxJ++) {
      const unsigned int nAuxJ = auxBasisFunctionsEnv[auxJ]->getNContracted();
      unsigned int offauxJ = auxBasisControllerEnv->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
      for (unsigned int auxK = 0; auxK < auxBasisFunctions.size(); auxK++) {
        Eigen::MatrixXd intDerivs;
#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int nAuxK = auxBasisFunctions[auxK]->getNContracted();
        unsigned int offauxK = auxBasisController->extendedIndex(auxK);
        if (libint.compute(LIBINT_OPERATOR::coulomb, 1, *auxBasisFunctionsEnv[auxJ], *auxBasisFunctions[auxK], intDerivs)) {
          double perm = 1.0;

          Eigen::MatrixXd prefac =
              0.5 * perm * Dvec.segment(offauxJ + nAuxBasFunc, nAuxJ) * Dvec.segment(offauxK, nAuxK).transpose();
          for (unsigned int iAtom = 0; iAtom < 2; ++iAtom) {
            unsigned int nAtom;
            switch (iAtom) {
              case (0):
                nAtom = auxMappingEnv[auxJ] + atomsAct.size();
                break;
              case (1):
                nAtom = auxMapping[auxK];
                break;
            }
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom * 3 + iDirection).data(), nAuxK, nAuxJ);
              eriContrPriv[threadId](nAtom, iDirection) -= tmp.transpose().cwiseProduct(prefac).sum();
            }
          }
        }
      }
    }

    for (unsigned int auxJ = 0; auxJ < auxBasisFunctions.size(); auxJ++) {
      const unsigned int nAuxJ = auxBasisFunctions[auxJ]->getNContracted();
      unsigned int offauxJ = auxBasisController->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
      for (unsigned int auxK = 0; auxK < auxBasisFunctions.size(); auxK++) {
        Eigen::MatrixXd intDerivs;
#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        const unsigned int nAuxK = auxBasisFunctions[auxK]->getNContracted();
        unsigned int offauxK = auxBasisController->extendedIndex(auxK);
        if (libint.compute(LIBINT_OPERATOR::coulomb, 1, *auxBasisFunctions[auxJ], *auxBasisFunctions[auxK], intDerivs)) {
          double perm = 1.0;

          Eigen::MatrixXd prefac = 0.5 * perm * Dvec.segment(offauxJ, nAuxJ) * Dvec.segment(offauxK, nAuxK).transpose();
          for (unsigned int iAtom = 0; iAtom < 2; ++iAtom) {
            unsigned int nAtom;
            switch (iAtom) {
              case (0):
                nAtom = auxMapping[auxJ];
                break;
              case (1):
                nAtom = auxMapping[auxK];
                break;
            }
            for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
              Eigen::Map<Eigen::MatrixXd> tmp(intDerivs.col(iAtom * 3 + iDirection).data(), nAuxK, nAuxJ);
              eriContrPriv[threadId](nAtom, iDirection) -= tmp.transpose().cwiseProduct(prefac).sum();
            }
          }
        }
      }
    }

    libint.finalize(LIBINT_OPERATOR::coulomb, 1, 2);

    /*
     * The following set of loops calculates the complete 1st and 2nd part of eq ([1].[6])
     * (easily recognizable by the sign).
     *
     * Check RICoulombPotential::getGeomGradients()!
     *
     * Has to be done twice in order to get the combinations
     * auxBasisEnv/Basis and auxBasis/BasisEnv.
     */
    TwoElecThreeCenterIntLooper derivLooper2(LIBINT_OPERATOR::coulomb, 1, basisControllerEnv, auxBasisController, screeningEnv);
    auto const add3cDerivs2 = [&eriContrPriv, &Dvec, &densityMatrixEnv, &mappingEnv, &auxMapping, &basisControllerEnv,
                               &auxBasisController, &nAtomsAct](const unsigned int& i, const unsigned int& j,
                                                                const unsigned int& K, Eigen::VectorXd& intValues,
                                                                const unsigned int threadId) {
      double perm = (i == j ? 1.0 : 2.0);
      for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
        for (unsigned int iAtom = 0; iAtom < 3; ++iAtom) {
          unsigned int nAtom;
          switch (iAtom) {
            case (0):
              nAtom = auxMapping[auxBasisController->reducedIndex(K)];
              break;
            case (1):
              nAtom = mappingEnv[basisControllerEnv->reducedIndex(i)] + nAtomsAct;
              break;
            case (2):
              nAtom = mappingEnv[basisControllerEnv->reducedIndex(j)] + nAtomsAct;
              break;
          }
          eriContrPriv[threadId](nAtom, iDirection) +=
              perm * Dvec[K] * densityMatrixEnv(i, j) * intValues(iAtom * 3 + iDirection);
        }
      }
    };
    derivLooper2.loop(add3cDerivs2, maxDens.maxCoeff());

    TwoElecThreeCenterIntLooper derivLooper1(LIBINT_OPERATOR::coulomb, 1, basisController, auxBasisControllerEnv, screeningAct);
    auto const add3cDerivs1 = [&eriContrPriv, &Dvec, &densityMatrix, &mapping, &auxMappingEnv, &basisController,
                               &auxBasisControllerEnv, &nAtomsAct,
                               &nAuxBasFunc](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                             Eigen::VectorXd& intValues, const unsigned int threadId) {
      double perm = (i == j ? 1.0 : 2.0);
      for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
        for (unsigned int iAtom = 0; iAtom < 3; ++iAtom) {
          unsigned int nAtom;
          switch (iAtom) {
            case (0):
              nAtom = auxMappingEnv[auxBasisControllerEnv->reducedIndex(K)] + nAtomsAct;
              break;
            case (1):
              nAtom = mapping[basisController->reducedIndex(i)];
              break;
            case (2):
              nAtom = mapping[basisController->reducedIndex(j)];
              break;
          }
          eriContrPriv[threadId](nAtom, iDirection) +=
              perm * Dvec[K + nAuxBasFunc] * densityMatrix(i, j) * intValues(iAtom * 3 + iDirection);
        }
      }
    };
    derivLooper1.loop(add3cDerivs1, maxDens.maxCoeff());

    TwoElecThreeCenterIntLooper derivLooper3(LIBINT_OPERATOR::coulomb, 1, basisController, auxBasisController, screeningAct);
    auto const add3cDerivs3 = [&eriContrPriv, &Dvec, &densityMatrix, &mapping, &auxMapping, &basisController,
                               &auxBasisController](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                                    Eigen::VectorXd& intValues, const unsigned int threadId) {
      double perm = (i == j ? 1.0 : 2.0);
      for (unsigned int iDirection = 0; iDirection < 3; ++iDirection) {
        for (unsigned int iAtom = 0; iAtom < 3; ++iAtom) {
          unsigned int nAtom;
          switch (iAtom) {
            case (0):
              nAtom = auxMapping[auxBasisController->reducedIndex(K)];
              break;
            case (1):
              nAtom = mapping[basisController->reducedIndex(i)];
              break;
            case (2):
              nAtom = mapping[basisController->reducedIndex(j)];
              break;
          }
          eriContrPriv[threadId](nAtom, iDirection) +=
              perm * Dvec[K] * densityMatrix(i, j) * intValues(iAtom * 3 + iDirection);
        }
      }
    };
    derivLooper3.loop(add3cDerivs3, maxDens.maxCoeff());

#ifdef _OPENMP
    for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
      eriContr += eriContrPriv[i];
    }
#else
    eriContr = eriContrPriv[0];
#endif

    for (unsigned int i = 0; i < atomsAct.size(); ++i) {
      for (unsigned int j = 0; j < 3; ++j) {
        activeSystemGradientContr(i, j) += eriContr(i, j);
      }
    }
  }
  // TODO: The above scheme also includes inner system Coulomb gradient contributions once per environment system,
  // hence this has to be explicitly calculated and subtracted in the end. This is not incredibly expensive,
  // but unnecessary and, alas, very hard to extract from the above algorithm.
  RICoulombPotential<SCFMode> coulPot(
      actSystem, actSystem->template getElectronicStructure<SCFMode>()->getDensityMatrixController(),
      RI_J_IntegralControllerFactory::getInstance().produce(
          actSystem->getAtomCenteredBasisController(),
          actSystem->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB)),
      actSystem->getSettings().basis.integralThreshold, actSystem->getSettings().basis.integralIncrementThresholdStart,
      actSystem->getSettings().basis.integralIncrementThresholdEnd, actSystem->getSettings().basis.incrementalSteps);
  auto ElEl = coulPot.getGeomGradients();
  return activeSystemGradientContr - (_envSystems.size() * ElEl);
}

template class CoulombInteractionPotential<Options::SCF_MODES::RESTRICTED>;
template class CoulombInteractionPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
