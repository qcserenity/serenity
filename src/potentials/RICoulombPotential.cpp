/**
 * @file   RICoulombPotential.cpp
 *
 * @date   Nov 24, 2016
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
#include "potentials/RICoulombPotential.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"        //Gradients.
#include "data/ElectronicStructure.h"                 //Density matrix.
#include "data/OrbitalController.h"                   //Gradients.
#include "integrals/RI_J_IntegralControllerFactory.h" //Construct RIJ integral controller.
#include "integrals/wrappers/Libint.h"                //getSharedPtr()
#include "io/FormattedOutputStream.h"                 //Filtered output streams.
#include "misc/Timing.h"                              //Timings.
#include "system/SystemController.h"                  //gradients.

namespace Serenity {

template<Options::SCF_MODES SCFMode>
RICoulombPotential<SCFMode>::RICoulombPotential(std::shared_ptr<SystemController> actSystem,
                                                std::shared_ptr<DensityMatrixController<SCFMode>> dMat,
                                                std::shared_ptr<RI_J_IntegralController> ri_j_IntController,
                                                const double prescreeningThreshold, double prescreeningIncrementStart,
                                                double prescreeningIncrementEnd, unsigned int incrementSteps)
  : Potential<SCFMode>(dMat->getDensityMatrix().getBasisController()),
    _actSystem(actSystem),
    _dMatController(dMat),
    _ri_j_IntController(ri_j_IntController),
    _fullpotential(nullptr),
    _outOfDate(true),
    _incrementHelper(std::make_shared<IncrementalFockMatrix<SCFMode>>(
        dMat, (prescreeningThreshold != 0) ? prescreeningThreshold : this->_basis->getPrescreeningThreshold(),
        prescreeningIncrementStart, prescreeningIncrementEnd, incrementSteps, "RI-J Coulomb")) {
  this->_basis->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  this->_dMatController->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  _fullpotential = std::make_shared<FockMatrix<SCFMode>>(FockMatrix<SCFMode>(this->_basis));
  auto& temp = *_fullpotential;
  for_spin(temp) {
    temp_spin.setZero();
  };
  _screening = prescreeningIncrementStart;
};

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& RICoulombPotential<SCFMode>::getMatrix() {
  Timings::takeTime("Active System -    Coulomb Pot.");
  if (_outOfDate) {
    DensityMatrix<SCFMode> densityMatrix(this->_basis);
    std::vector<std::shared_ptr<FockMatrix<SCFMode>>> matrices = {_fullpotential};
    _incrementHelper->updateDensityAndThreshold(densityMatrix, _screening, matrices);
    this->addToMatrix(*_fullpotential, densityMatrix);
    _outOfDate = false;
  }
  Timings::timeTaken("Active System -    Coulomb Pot.");
  return *_fullpotential;
}

template<Options::SCF_MODES SCFMode>
double RICoulombPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  if (_outOfDate)
    this->getMatrix();
  Timings::takeTime("Active System -    Coulomb Pot.");
  auto& pot = *_fullpotential;
  double energy = 0.0;
  for_spin(pot, P) {
    energy += 0.5 * pot_spin.cwiseProduct(P_spin).sum();
  };
  Timings::timeTaken("Active System -    Coulomb Pot.");
  return energy;
};

template<Options::SCF_MODES SCFMode>
void RICoulombPotential<SCFMode>::addToMatrix(FockMatrix<SCFMode>& F, const DensityMatrix<SCFMode>& deltaP) {
  unsigned nThreads = omp_get_max_threads();
  std::vector<MatrixInBasis<RESTRICTED>> fock_threads(nThreads, MatrixInBasis<RESTRICTED>(this->_basis));
  MatrixInBasis<RESTRICTED> totalDensity = deltaP.total();
  auto densptr = totalDensity.data();

  unsigned nx = _ri_j_IntController->getAuxBasisController()->getNBasisFunctions();
  unsigned nb = totalDensity.rows();

  auto maxDens = totalDensity.shellWiseAbsMax();
  auto maxDensPtr = maxDens.data();
  unsigned ns = maxDens.rows();

  Eigen::MatrixXd sumMat = Eigen::MatrixXd::Zero(nx, nThreads);
  auto sumptr = sumMat.data();

  // Distribute function for first contraction.
  auto distribute1 = [&](unsigned i, unsigned j, unsigned K, double integral, unsigned threadId) {
    double perm = (i == j ? 1.0 : 2.0);
    sumptr[nx * threadId + K] += perm * integral * densptr[i * nb + j];
  };

  // Prescreening function for first contraction.
  auto prescreen1 = [&](const unsigned i, const unsigned j, unsigned int, const double schwarz) {
    unsigned long ij = i * ns + j;
    return (maxDensPtr[ij] * schwarz < _screening);
  };

  // Perform first contraction.
  _ri_j_IntController->loopOver3CInts(distribute1, prescreen1);

  // Solve linear system.
  Eigen::VectorXd sumPMuNu_DMuNu = sumMat.rowwise().sum();
  sumMat.resize(1, 1);
  Eigen::VectorXd coefficients = _ri_j_IntController->getLLTMetric().solve(sumPMuNu_DMuNu).eval();

  // Distribute function for second contraction.
  auto coeff = coefficients.data();
  auto coeffMax = _ri_j_IntController->getAuxBasisController()->shellWiseAbsMax(coefficients);
  auto distribute2 = [&](unsigned i, unsigned j, unsigned K, double integral, unsigned threadId) {
    // Only perform half the needed contractions and symmetrize afterwards.
    fock_threads[threadId].data()[i * nb + j] += integral * coeff[K];
  };

  // Prescreening function for second contraction.
  auto prescreen2 = [&](const unsigned, const unsigned, const unsigned K, const double schwarz) {
    return (coeffMax(K) * schwarz < _screening);
  };

  // Perform second contraction.
  _ri_j_IntController->loopOver3CInts(distribute2, prescreen2);

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

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd RICoulombPotential<SCFMode>::getGeomGradients() {
  auto actSystem = _actSystem.lock();
  auto atoms = actSystem->getAtoms();
  unsigned int nAtoms = atoms.size();
  Matrix<double> eriContr(nAtoms, 3);
  Matrix<double> exchangeContr(nAtoms, 3);
  eriContr.setZero();
  exchangeContr.setZero();

  // Collecting all the necessities, Mapping Basis and AuxBasis to atoms

  auto basisController = actSystem->getAtomCenteredBasisController();
  auto auxBasisController = actSystem->getAtomCenteredBasisController(Options::BASIS_PURPOSES::AUX_COULOMB);
  bool normAux = !(auxBasisController->isAtomicCholesky());

  auto basisIndicesRed = actSystem->getAtomCenteredBasisController()->getBasisIndicesRed();
  auto mapping = actSystem->getAtomCenteredBasisController()->getAtomIndicesOfBasisShells();

  auto& auxBasis = auxBasisController->getBasis();
  auto nAuxBasFunc = auxBasisController->getNBasisFunctions();
  auto auxBasisIndicesRed = auxBasisController->getBasisIndicesRed();
  auto auxMapping = auxBasisController->getAtomIndicesOfBasisShells();

  DensityMatrix<RESTRICTED> densityMatrix(actSystem->template getElectronicStructure<SCFMode>()->getDensityMatrix().total());

  // Get the existing RI_J Controller from Config (there should be one from the first SCF)

  auto ri_j_IntController = RI_J_IntegralControllerFactory::getInstance().produce(basisController, auxBasisController);

  auto inverseM = ri_j_IntController->getInverseM();

  /*
   * Constructing the different parts of eq. ([2].[6]).
   *
   * Starting with the C vectors. The first step combines the
   * Density Matrix elements ij with the respecting 3c Integral <J|ij>.
   * This is technically a Matrix-Vector multiplication.
   *
   */
#ifdef _OPENMP
  // create a vector of vectors for each thread
  std::vector<Eigen::VectorXd> CvecPriv(omp_get_max_threads(), Eigen::VectorXd::Zero(nAuxBasFunc));
#else
  // or just one
  std::vector<Eigen::VectorXd> CvecPriv(1, Eigen::VectorXd::Zero(nAuxBasFunc));
#endif

  TwoElecThreeCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, basisController, auxBasisController,
                                     _incrementHelper->getPrescreeningThreshold());

  Eigen::VectorXd Cvec(nAuxBasFunc);
  Cvec.setZero();
  auto const calcCVector = [&CvecPriv, &densityMatrix](const unsigned int& i, const unsigned int& j, const unsigned int& K,
                                                       Eigen::VectorXd& intValues, const unsigned int threadId) {
    const double perm = (i == j ? 1.0 : 2.0);
    CvecPriv[threadId][K] += perm * intValues(0) * densityMatrix(i, j);
  };
  looper.loop(calcCVector);

#ifdef _OPENMP
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    Cvec += CvecPriv[i];
  }
#else
  Cvec += CvecPriv[0];
#endif

  /*
   * The vector resulting from this Matrix * Vector multiplication is then
   * part of another Matrix-Vector Multiplication with the inverted 2c Integral Matrix (inverseM)
   */

  Eigen::VectorXd Dvec = inverseM * Cvec;

  /*
   * The resulting Vector Dvec is the equivalent of the C^J vector in the paper.
   * The following set of loops calculates the complete 3rd part of eq ([2].[6])
   * (easily recognizable by the sign).
   *
   * This is done analogously to the calculation of the inverseM in the RI_J_IntegralController,
   * but with the added cartesian directions. Note that to get a certain part of an Eigen Matrix,
   * the method .block() is invoked, the equivalent for a vector is .segment().
   *
   */
  Libint& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::coulomb, 1, 2);
#ifdef _OPENMP
  std::vector<Eigen::MatrixXd> eriContrPriv(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
  std::vector<Eigen::MatrixXd> eriContrPrivAlt(omp_get_max_threads(), Eigen::MatrixXd::Zero(nAtoms, 3));
#else
  // or just one
  std::vector<Eigen::MatrixXd> eriContrPriv(1, Eigen::MatrixXd::Zero(nAtoms, 3));
#endif

  for (unsigned int auxJ = 0; auxJ < auxBasis.size(); auxJ++) {
    const unsigned int nAuxJ = auxBasis[auxJ]->getNContracted();
    unsigned int offauxJ = auxBasisController->extendedIndex(auxJ);
#pragma omp parallel for schedule(static, 1)
    for (unsigned int auxK = 0; auxK <= auxJ; auxK++) {
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      Eigen::MatrixXd intDerivs;
      const unsigned int nAuxK = auxBasis[auxK]->getNContracted();
      unsigned int offauxK = auxBasisController->extendedIndex(auxK);

      if (libint.compute(LIBINT_OPERATOR::coulomb, 1, *auxBasis[auxJ], *auxBasis[auxK], intDerivs, normAux)) {
        double perm = (auxJ == auxK ? 1.0 : 2.0);

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
   *
   * The following looper calculates the complete 1st and 2nd part of eq ([2].[6])
   * (easily recognizable by the sign).
   *
   */

  TwoElecThreeCenterIntLooper derivLooper(LIBINT_OPERATOR::coulomb, 1, basisController, auxBasisController,
                                          _incrementHelper->getPrescreeningThreshold());

  auto const add3cDerivs = [&eriContrPriv, &Dvec, &densityMatrix, &mapping, &auxMapping, &basisController, &auxBasisController,
                            &eriContrPrivAlt](const unsigned int& i, const unsigned int& j, const unsigned int& K,
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
        eriContrPriv[threadId](nAtom, iDirection) += perm * Dvec[K] * densityMatrix(i, j) * intValues(iAtom * 3 + iDirection);
        eriContrPrivAlt[threadId](nAtom, iDirection) +=
            perm * Dvec[K] * densityMatrix(i, j) * intValues(iAtom * 3 + iDirection);
      }
    }
  };
  derivLooper.loop(add3cDerivs);

#ifdef _OPENMP
  for (unsigned int i = 0; i < (unsigned int)omp_get_max_threads(); ++i) {
    eriContr += eriContrPriv[i];
  }
#else
  eriContr = eriContrPriv[0];
#endif
  Eigen::MatrixXd eriContrAlt = Eigen::MatrixXd::Zero(nAtoms, 3);
  for (const auto& deriv : eriContrPrivAlt) {
    eriContrAlt += deriv;
  }
  return std::move(eriContr);
}

template class RICoulombPotential<Options::SCF_MODES::RESTRICTED>;
template class RICoulombPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
