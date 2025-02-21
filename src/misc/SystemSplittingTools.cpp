/**
 * @file SystemSplittingTools.cpp
 *
 * @date Sep 25, 2018
 * @author Moritz Bensberg
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
#include "misc/SystemSplittingTools.h"
/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/MullikenPopulationCalculator.h"
#include "basis/AtomCenteredBasisController.h"
#include "basis/Basis.h"
#include "basis/BasisFunctionMapper.h" //Coefficient matrix resorting.
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Geometry.h"
#include "grid/AtomCenteredGridController.h"
#include "integrals/OneElectronIntegralController.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutputStream.h"
#include "settings/Settings.h"
#include "system/SystemController.h"

namespace Serenity {
template<Options::SCF_MODES SCFMode>
MatrixInBasis<SCFMode> SystemSplittingTools<SCFMode>::projectMatrixIntoNewBasis(
    const MatrixInBasis<SCFMode>& oldMatrix, std::shared_ptr<BasisController> newBasis,
    std::shared_ptr<MatrixInBasis<Options::SCF_MODES::RESTRICTED>> newOverlap) {
  auto& libint = Libint::getInstance();
  Eigen::MatrixXd overlapB;
  if (newOverlap) {
    overlapB = *newOverlap;
  }
  else {
    overlapB = libint.compute1eInts(LIBINT_OPERATOR::overlap, newBasis, newBasis);
  }
  Eigen::MatrixXd overlapAB = libint.compute1eInts(LIBINT_OPERATOR::overlap, oldMatrix.getBasisController(), newBasis);
  // Calculate inverse of AO overlap integrals in basis B. Use SVD,
  // i.e. S_B = U * D * V^T, since overlapB could be ill-conditioned
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(overlapB, Eigen::ComputeThinU | Eigen::ComputeThinV);
  svd.setThreshold(1e-6);
  // Calculate projection operator P_{BA}
  Eigen::MatrixXd projectionOperator = svd.solve(overlapAB);
  MatrixInBasis<SCFMode> newMatrix(newBasis);
  for_spin(newMatrix, oldMatrix) {
    newMatrix_spin = (projectionOperator * oldMatrix_spin * projectionOperator.transpose()).eval();
  };
  return newMatrix;
}
template<Options::SCF_MODES SCFMode>
Eigen::SparseMatrix<int>
SystemSplittingTools<SCFMode>::reduceMatrixToMullikenNetPopulationMap(const Eigen::MatrixXd& matrix, double mnpThreshold) {
  std::vector<Eigen::Triplet<int>> tripletList;
  for (unsigned int x = 0; x < matrix.cols(); ++x) {
    for (unsigned int coef = 0; coef < matrix.rows(); ++coef) {
      if (matrix(coef, x) * matrix(coef, x) > mnpThreshold)
        tripletList.push_back(Eigen::Triplet<int>(coef, x, 1));
    } // for coef
  }   // for x
  Eigen::SparseMatrix<int> mnpMap(matrix.rows(), matrix.cols());
  mnpMap.setFromTriplets(tripletList.begin(), tripletList.end());
  return mnpMap;
}

template<Options::SCF_MODES SCFMode>
unsigned int SystemSplittingTools<SCFMode>::matchAtom(std::shared_ptr<Geometry> geometry, std::shared_ptr<Atom> atom) {
  unsigned int actIndex = 0;
  for (const auto& actAtom : geometry->getAtoms()) {
    if (*actAtom == *makeAtomFromDummy(atom)) {
      return actIndex;
    }
    ++actIndex;
  }
  ++actIndex;
  // if not found actIndex > nActAtoms
  return actIndex;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<DensityMatrix<SCFMode>>
SystemSplittingTools<SCFMode>::buildNonOrthogonalDensityMatrix(std::shared_ptr<SystemController> environmentSystem,
                                                               SpinPolarizedData<SCFMode, std::vector<bool>> distantOrbitals) {
  const auto& envCoeff = environmentSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  double nOcc = (SCFMode == Options::SCF_MODES::UNRESTRICTED) ? 1.0 : 2.0;
  auto nonOrthoDensMatPtr = std::make_shared<DensityMatrix<SCFMode>>(environmentSystem->getBasisController());
  auto& nonOrthoDensMat = *nonOrthoDensMatPtr;
  for_spin(envCoeff, distantOrbitals, nonOrthoDensMat) {
    nonOrthoDensMat_spin.setZero();
    unsigned int nNonOrtho = 0;
    for (const auto& nonOrtho : distantOrbitals_spin)
      if (nonOrtho)
        ++nNonOrtho;
    if (nNonOrtho > 0) {
      Eigen::MatrixXd nonOrthoCoeff(environmentSystem->getBasisController()->getNBasisFunctions(), nNonOrtho);
      unsigned int nAssigned = 0;
      for (unsigned int coeffIndex = 0; coeffIndex < distantOrbitals_spin.size(); ++coeffIndex) {
        if (distantOrbitals_spin[coeffIndex]) {
          nonOrthoCoeff.block(0, nAssigned, nonOrthoCoeff.rows(), 1) = envCoeff_spin.col(coeffIndex);
          ++nAssigned;
        } // if distantOrbitals_spin[coeffIndex]
      }   // for coeffIndex
      // Build density matrix as C*C^T * 2.0 or 1.0 (Restricted/Unrestricted)
      nonOrthoDensMat_spin = nOcc * (nonOrthoCoeff * nonOrthoCoeff.transpose());
    } // if nNonOrtho > 0
  };
  return nonOrthoDensMatPtr;
}

template<Options::SCF_MODES SCFMode>
SpinPolarizedData<SCFMode, std::vector<bool>>
SystemSplittingTools<SCFMode>::selectDistantOrbitals(SPMatrix<SCFMode>& orbitalPopulations,
                                                     std::shared_ptr<SystemController> activeSystem,
                                                     std::shared_ptr<SystemController> environmentSystem,
                                                     double basisFunctionRatio, double borderAtomThreshold) {
  SpinPolarizedData<SCFMode, std::vector<bool>> distantOrbitals;
  unsigned int nEnvAtoms = environmentSystem->getGeometry()->getNAtoms();
  /*
   * 1. Match atoms and select distant atoms which have their number of shells reduced significantly.
   */
  auto actAtoms = activeSystem->getGeometry()->getAtoms();
  const auto& actSystemBasisLoadingData = activeSystem->getAtomCenteredBasisController()->getBasisLoadingData();

  std::vector<bool> distantAtoms(nEnvAtoms, true);
  unsigned int envAtomIndex = 0;
  for (const auto& envAtom : environmentSystem->getGeometry()->getAtoms()) {
    unsigned int matchedIndex = matchAtom(activeSystem->getGeometry(), envAtom);
    if (matchedIndex < actAtoms.size()) {
      // Check remaining basis function shells in the active system.
      double remainingBasisFunctionRatio = (double)actSystemBasisLoadingData[matchedIndex].shells.sum() /
                                           (double)actSystemBasisLoadingData[matchedIndex].shells.size();
      if (remainingBasisFunctionRatio >= basisFunctionRatio)
        distantAtoms[envAtomIndex] = false;
    } // if
    ++envAtomIndex;
  } // for envAtom
  /*
   * 2. Select the orbitals which are localized on not distant atom
   */
  // search for orbitals located on "distant" atoms
  // by getting the sum of Mulliken populations on not "distant" atoms.
  // If this sum is lower than the given threshold, the complete orbital is considered to be distant.
  const auto& nOccEnv = environmentSystem->getNOccupiedOrbitals<SCFMode>();
  unsigned int numberOfDistantOrbitals = 0;
  for_spin(distantOrbitals, nOccEnv, orbitalPopulations) {
    distantOrbitals_spin.resize(nOccEnv_spin, false);
    // Loop over all environment orbitals and atoms and evaluate the population on the non-distant atoms.
    // Sanity check for dimension of populations and environment atoms.
    assert((int)orbitalPopulations_spin.rows() == (int)distantAtoms.size() &&
           "Dimensions of atoms and their populations are not fitting!");
    for (unsigned int iEnvOrb = 0; iEnvOrb < nOccEnv_spin; ++iEnvOrb) {
      double popOnNotDistantAtoms = 0.0;
      for (unsigned int atomIndex = 0; atomIndex < distantAtoms.size(); ++atomIndex) {
        if (!distantAtoms[atomIndex]) {
          popOnNotDistantAtoms += std::fabs(orbitalPopulations_spin(atomIndex, iEnvOrb));
        } // if !distantAtoms[atomIndex]
      }   // for atomIndex
      if (popOnNotDistantAtoms < borderAtomThreshold) {
        distantOrbitals_spin[iEnvOrb] = true;
        ++numberOfDistantOrbitals;
      }
    } // for iEnvOrb
  };
  OutputControl::nOut << std::endl;
  OutputControl::nOut << "Number of orbitals considered to be distant (alpha+beta): " << numberOfDistantOrbitals << std::endl;
  OutputControl::nOut << std::endl;
  return distantOrbitals;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<ElectronicStructure<SCFMode>> SystemSplittingTools<SCFMode>::resortBasisSetOfElectronicStructure(
    std::shared_ptr<ElectronicStructure<SCFMode>> electronicStructure, std::shared_ptr<BasisController> newBasisController,
    std::shared_ptr<OneElectronIntegralController> oneElectronIntegralController) {
  auto oldBasisController = electronicStructure->getDensityMatrix().getBasisController();
  BasisFunctionMapper basisFuncMapperOldToNew(oldBasisController);
  std::shared_ptr<Eigen::SparseMatrix<double>> projection = basisFuncMapperOldToNew.getSparseProjection(newBasisController);
  auto orbitalController = electronicStructure->getMolecularOrbitals();
  auto oldCoefficientMatrix = orbitalController->getCoefficients();
  auto oldEigenvalues = orbitalController->getEigenvalues();
  auto oldCoreOrbitals = orbitalController->getOrbitalFlags();
  auto newCoefficientMatrixPtr = std::make_unique<CoefficientMatrix<SCFMode>>(newBasisController);
  auto newEigenvaluesPtr =
      std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(newBasisController->getNBasisFunctions());
  auto newCoreOrbitalsPtr =
      std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>(newBasisController->getNBasisFunctions());
  auto nOcc = electronicStructure->getNOccupiedOrbitals();
  CoefficientMatrix<SCFMode>& newCoefficientMatrix = *newCoefficientMatrixPtr;
  SpinPolarizedData<SCFMode, Eigen::VectorXd>& newEigenvalues = *newEigenvaluesPtr;
  SpinPolarizedData<SCFMode, Eigen::VectorXi>& newCoreOrbitals = *newCoreOrbitalsPtr;
  unsigned int nOrbEigenvalues =
      std::min(newBasisController->getNBasisFunctions(), oldBasisController->getNBasisFunctions());
  for_spin(oldCoefficientMatrix, oldEigenvalues, newCoefficientMatrix, newEigenvalues, newCoreOrbitals, oldCoreOrbitals) {
    newCoefficientMatrix_spin = *projection * oldCoefficientMatrix_spin;
    newEigenvalues_spin.head(nOrbEigenvalues) = oldEigenvalues_spin.head(nOrbEigenvalues);
    newCoreOrbitals_spin.head(nOrbEigenvalues) = oldCoreOrbitals_spin.head(nOrbEigenvalues);
  };
  auto newOrbitalController = std::make_shared<OrbitalController<SCFMode>>(
      std::move(newCoefficientMatrixPtr), newBasisController, std::move(newEigenvaluesPtr), std::move(newCoreOrbitalsPtr));
  auto newElectronicStructure =
      std::make_shared<ElectronicStructure<SCFMode>>(newOrbitalController, oneElectronIntegralController, nOcc);
  newElectronicStructure->getDensityMatrixController()->updateDensityMatrix();
  return newElectronicStructure;
}

template<Options::SCF_MODES SCFMode>
std::pair<std::shared_ptr<ElectronicStructure<SCFMode>>, std::shared_ptr<ElectronicStructure<SCFMode>>>
SystemSplittingTools<SCFMode>::splitElectronicStructure(std::shared_ptr<SystemController> system,
                                                        SpinPolarizedData<SCFMode, std::vector<bool>>& activeOrbitals) {
  // get data
  const auto basisController = system->getBasisController();
  const unsigned int nBasisFunctions = basisController->getNBasisFunctions();
  const auto& superCoefficients = system->getActiveOrbitalController<SCFMode>()->getCoefficients();
  const auto& superEigenvalues = system->getActiveOrbitalController<SCFMode>()->getEigenvalues();
  const auto& superCoreOrbitals = system->getActiveOrbitalController<SCFMode>()->getOrbitalFlags();
  // Build coefficient matrix and eigenvalue vectors
  auto actCoeffPtr = std::make_unique<CoefficientMatrix<SCFMode>>(basisController);
  auto actEigenPtr = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(nBasisFunctions);
  auto actCoreOPtr =
      std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>(Eigen::VectorXi::Constant(nBasisFunctions, 3));
  auto envCoeffPtr = std::make_unique<CoefficientMatrix<SCFMode>>(basisController);
  auto envEigenPtr = std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(nBasisFunctions);
  auto envCoreOPtr =
      std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>(Eigen::VectorXi::Constant(nBasisFunctions, 3));

  // Split coefficients
  auto& actCoefficients = *actCoeffPtr;
  auto& envCoefficients = *envCoeffPtr;
  for_spin(activeOrbitals, actCoefficients, superCoefficients, envCoefficients) {
    actCoefficients_spin.setZero();
    envCoefficients_spin.setZero();
    unsigned int aIndex = 0;
    unsigned int eIndex = 0;
    // occupied orbitals
    for (unsigned int i = 0; i < activeOrbitals_spin.size(); ++i) {
      if (activeOrbitals_spin[i]) {
        actCoefficients_spin.col(aIndex) = superCoefficients_spin.col(i);
        ++aIndex;
      }
      else {
        envCoefficients_spin.col(eIndex) = superCoefficients_spin.col(i);
        ++eIndex;
      }
    }
    // virtual orbital
    for (unsigned int i = activeOrbitals_spin.size(); i < nBasisFunctions; ++i) {
      actCoefficients_spin.col(aIndex) = superCoefficients_spin.col(i);
      envCoefficients_spin.col(eIndex) = superCoefficients_spin.col(i);
      ++aIndex;
      ++eIndex;
    }
  };
  // Split eigenvalues
  auto& actEigenvalues = *actEigenPtr;
  auto& envEigenvalues = *envEigenPtr;
  auto& actCoreOrbitals = *actCoreOPtr;
  auto& envCoreOrbitals = *envCoreOPtr;
  for_spin(activeOrbitals, superEigenvalues, actEigenvalues, envEigenvalues, actCoreOrbitals, envCoreOrbitals, superCoreOrbitals) {
    actEigenvalues_spin = Eigen::VectorXd::Constant(nBasisFunctions, std::numeric_limits<double>::infinity());
    envEigenvalues_spin = Eigen::VectorXd::Constant(nBasisFunctions, std::numeric_limits<double>::infinity());
    unsigned int aIndex = 0;
    unsigned int eIndex = 0;
    // occupied orbitals
    for (unsigned int i = 0; i < activeOrbitals_spin.size(); ++i) {
      if (activeOrbitals_spin[i]) {
        actEigenvalues_spin(aIndex) = superEigenvalues_spin(i);
        actCoreOrbitals_spin(aIndex) = superCoreOrbitals_spin(i);
        ++aIndex;
      }
      else {
        envEigenvalues_spin(eIndex) = superEigenvalues_spin(i);
        envCoreOrbitals_spin(eIndex) = superCoreOrbitals_spin(i);
        ++eIndex;
      }
    }
    // virtual orbital
    for (unsigned int i = activeOrbitals_spin.size(); i < nBasisFunctions; ++i) {
      actEigenvalues_spin(aIndex) = superEigenvalues_spin(i);
      envEigenvalues_spin(eIndex) = superEigenvalues_spin(i);
      actCoreOrbitals_spin(aIndex) = superCoreOrbitals_spin(i);
      envCoreOrbitals_spin(eIndex) = superCoreOrbitals_spin(i);
      ++aIndex;
      ++eIndex;
    }
  };
  // Number of occupied orbitals
  auto occupations = getNOccupiedOrbitals(activeOrbitals);
  // Build orbital controller
  auto activeOrbitalSet = std::make_shared<OrbitalController<SCFMode>>(std::move(actCoeffPtr), basisController,
                                                                       std::move(actEigenPtr), std::move(actCoreOPtr));
  auto actES = std::make_shared<ElectronicStructure<SCFMode>>(activeOrbitalSet, system->getOneElectronIntegralController(),
                                                              occupations.first);
  auto environmentOrbitalSet = std::make_shared<OrbitalController<SCFMode>>(
      std::move(envCoeffPtr), basisController, std::move(envEigenPtr), std::move(envCoreOPtr));
  auto envES = std::make_shared<ElectronicStructure<SCFMode>>(
      environmentOrbitalSet, system->getOneElectronIntegralController(), occupations.second);
  // update density matrices
  actES->getDensityMatrixController()->updateDensityMatrix();
  envES->getDensityMatrixController()->updateDensityMatrix();
  assert(actES->getDensityMatrix().getBasisController() == basisController);
  assert(envES->getDensityMatrix().getBasisController() == basisController);
  return std::pair<std::shared_ptr<ElectronicStructure<SCFMode>>, std::shared_ptr<ElectronicStructure<SCFMode>>>(actES, envES);
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<Geometry>>
SystemSplittingTools<SCFMode>::splitGeometry(std::shared_ptr<SystemController> system,
                                             const SpinPolarizedData<SCFMode, Eigen::VectorXi>& assignment,
                                             bool prioFirst, double locThreshold, unsigned int nFrag) {
  std::vector<std::vector<std::shared_ptr<Atom>>> partitioniedAtoms(nFrag, std::vector<std::shared_ptr<Atom>>());
  auto supersystemAtoms = system->getGeometry()->getAtoms();
  auto nOcc = system->template getNOccupiedOrbitals<SCFMode>();
  const unsigned int nAtoms = supersystemAtoms.size();
  Eigen::MatrixXd fragWiseAtomPopulations = Eigen::MatrixXd::Zero(nAtoms, nFrag);
  const auto& coefficients = system->template getActiveOrbitalController<SCFMode>()->getCoefficients();
  const auto& overlapMatrix = system->getOneElectronIntegralController()->getOverlapIntegrals();
  const auto& atomToBasis = system->getAtomCenteredBasisController()->getBasisIndices();
  for_spin(nOcc, assignment, coefficients) {
    auto orbitalWisePopulations = MullikenPopulationCalculator<RESTRICTED>::calculateAtomwiseOrbitalPopulations(
        coefficients_spin.leftCols(nOcc_spin), overlapMatrix, atomToBasis);
    for (unsigned int iOcc = 0; iOcc < nOcc_spin; ++iOcc) {
      fragWiseAtomPopulations.col(assignment_spin(iOcc)) += orbitalWisePopulations.col(iOcc);
    }
  };
  double occupation = (SCFMode == RESTRICTED) ? 2.0 : 1.0;
  for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
    unsigned int iFrag;
    fragWiseAtomPopulations.row(iAtom).maxCoeff(&iFrag);
    // If the first system is supposed to be prioritized in the selection,
    // check if the given threshold is exceeded.
    if (prioFirst && fragWiseAtomPopulations(iAtom, 0) > locThreshold / occupation)
      iFrag = 0;
    partitioniedAtoms[iFrag].push_back(supersystemAtoms[iAtom]);
  } // for iAtom
  std::vector<std::shared_ptr<Geometry>> geoms;
  for (auto& atomList : partitioniedAtoms) {
    geoms.push_back(std::make_shared<Geometry>(atomList));
  }
  return geoms;
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Atom> SystemSplittingTools<SCFMode>::makeDummyAtom(std::shared_ptr<Atom> atom) {
  auto oldAtomType = atom->getAtomType();
  auto atomType = std::make_shared<AtomType>(
      oldAtomType->getName() + ":", oldAtomType->getNuclearCharge(), oldAtomType->getMass(),
      oldAtomType->getBraggSlaterRadius(), oldAtomType->getVanDerWaalsRadius(), oldAtomType->getUFFRadius(),
      oldAtomType->getNCoreElectrons(), oldAtomType->getOccupations(), oldAtomType->getChemicalHardness(), true);
  std::string basisLabel = atom->getPrimaryBasisLabel();
  auto basisFunctions = atom->getBasisFunctions();

  std::pair<std::string, std::vector<std::shared_ptr<Shell>>> atomBasisFunctions(basisLabel, basisFunctions);

  auto dummyAtom = std::make_shared<Atom>(atomType, atom->getX(), atom->getY(), atom->getZ(), atomBasisFunctions);

  return (dummyAtom);
}

template<Options::SCF_MODES SCFMode>
std::shared_ptr<Atom> SystemSplittingTools<SCFMode>::makeAtomFromDummy(std::shared_ptr<Atom> atom) {
  auto atomType = atom->getAtomType();
  std::string nonDummyName = atomType->getName();
  if (atom->isDummy())
    nonDummyName = nonDummyName.substr(0, nonDummyName.size() - 1);
  auto newAtomType = std::make_shared<AtomType>(nonDummyName, atomType->getPSEPosition(), atomType->getMass(),
                                                atomType->getBraggSlaterRadius(), atomType->getVanDerWaalsRadius(),
                                                atomType->getUFFRadius(), atomType->getNCoreElectrons(),
                                                atomType->getOccupations(), atomType->getChemicalHardness());
  std::string basisLabel = atom->getPrimaryBasisLabel();
  auto basisFunctions = atom->getBasisFunctions();
  std::pair<std::string, std::vector<std::shared_ptr<Shell>>> atomBasisFunctions(basisLabel, basisFunctions);
  auto newAtom = std::make_shared<Atom>(newAtomType, atom->getX(), atom->getY(), atom->getZ(), atomBasisFunctions);
  return newAtom;
}

template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>
SystemSplittingTools<SCFMode>::getNOccupiedOrbitals(SpinPolarizedData<SCFMode, std::vector<bool>>& activeOrbitals) {
  std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>> occupations(
      SpinPolarizedData<SCFMode, unsigned int>(0), SpinPolarizedData<SCFMode, unsigned int>(0));
  auto& actOcc = occupations.first;
  auto& envOcc = occupations.second;
  for_spin(activeOrbitals, actOcc, envOcc) {
    for (const auto& orbBool : activeOrbitals_spin) {
      if (orbBool) {
        ++actOcc_spin;
      }
      else {
        ++envOcc_spin;
      }
    }
  };
  return occupations;
}

template<Options::SCF_MODES SCFMode>
std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>
SystemSplittingTools<SCFMode>::getNElectrons(
    std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>> occ) {
  auto nOrbsAct = occ.first;
  auto nOrbsEnv = occ.second;
  SpinPolarizedData<SCFMode, unsigned int> actElectrons = 0;
  SpinPolarizedData<SCFMode, unsigned int> envElectrons = 0;
  unsigned int nElectronsPerOrbital = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2 : 1;
  for_spin(nOrbsAct, nOrbsEnv, actElectrons, envElectrons) {
    actElectrons_spin += nElectronsPerOrbital * nOrbsAct_spin;
    envElectrons_spin += nElectronsPerOrbital * nOrbsEnv_spin;
  };
  return std::pair<SpinPolarizedData<SCFMode, unsigned int>, SpinPolarizedData<SCFMode, unsigned int>>(actElectrons,
                                                                                                       envElectrons);
}

template<>
int SystemSplittingTools<Options::SCF_MODES::RESTRICTED>::getSpin(SpinPolarizedData<Options::SCF_MODES::RESTRICTED, unsigned int> nOcc) {
  (void)nOcc;
  return 0;
}
template<>
int SystemSplittingTools<Options::SCF_MODES::UNRESTRICTED>::getSpin(
    SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, unsigned int> nOcc) {
  return (int)nOcc.alpha - (int)nOcc.beta;
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode>
SystemSplittingTools<SCFMode>::getMatrixBlock(const SPMatrix<SCFMode>& matrix, const std::vector<unsigned int> atomsA,
                                              const std::vector<unsigned int> atomsB,
                                              const std::vector<std::pair<unsigned int, unsigned int>>& atomBasisIndices) {
  // get the block sizes
  unsigned int nBasisFunctionsA = 0;
  unsigned int nBasisFunctionsB = 0;
  for (const auto atomIndexA : atomsA) {
    unsigned int blockRowStart = atomBasisIndices[atomIndexA].first;
    unsigned int blockRowEnd = atomBasisIndices[atomIndexA].second;
    unsigned int nRows = blockRowEnd - blockRowStart;
    nBasisFunctionsA += nRows;
  }
  for (const auto atomIndexB : atomsB) {
    unsigned int blockColsStart = atomBasisIndices[atomIndexB].first;
    unsigned int blockColsEnd = atomBasisIndices[atomIndexB].second;
    unsigned int nCols = blockColsEnd - blockColsStart;
    nBasisFunctionsB += nCols;
  }
  SPMatrix<SCFMode> toReturn(nBasisFunctionsA, nBasisFunctionsB);
  // Extract blocks from the matrix
  unsigned int rowIndex = 0;
  for (const auto& atomIndexA : atomsA) {
    const unsigned int blockRowStart = atomBasisIndices[atomIndexA].first;
    const unsigned int blockRowEnd = atomBasisIndices[atomIndexA].second;
    assert(blockRowEnd >= blockRowStart);
    const unsigned int nRows = blockRowEnd - blockRowStart;
    unsigned int columnIndex = 0;
    for (const auto& atomIndexB : atomsB) {
      const unsigned int blockColStart = atomBasisIndices[atomIndexB].first;
      const unsigned int blockColEnd = atomBasisIndices[atomIndexB].second;
      assert(blockColEnd >= blockColStart);
      const unsigned int nCols = blockColEnd - blockColStart;
      for_spin(toReturn, matrix) {
        toReturn_spin.block(rowIndex, columnIndex, nRows, nCols) =
            matrix_spin.block(blockRowStart, blockColStart, nRows, nCols);
      };
      columnIndex += nCols;
    }
    rowIndex += nRows;
  }
  return toReturn;
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode> SystemSplittingTools<SCFMode>::getMatrixBlockShellWise(const MatrixInBasis<SCFMode>& matrix,
                                                                         const Eigen::SparseVector<int> shellsA,
                                                                         const Eigen::SparseVector<int> shellsB) {
  const auto basisController = matrix.getBasisController();
  // get the block sizes
  unsigned int nBasisFunctionsA = 0;
  unsigned int nBasisFunctionsB = 0;
  const auto& basis = basisController->getBasis();
  for (Eigen::SparseVector<int>::InnerIterator it(shellsA); it; ++it) {
    nBasisFunctionsA += basis[it.row()]->getNContracted();
  }
  for (Eigen::SparseVector<int>::InnerIterator it(shellsB); it; ++it) {
    nBasisFunctionsB += basis[it.row()]->getNContracted();
  }
  SPMatrix<SCFMode> toReturn(nBasisFunctionsA, nBasisFunctionsB);
  // Extract blocks from the matrix
  unsigned int rowIndex = 0;
  for (Eigen::SparseVector<int>::InnerIterator itA(shellsA); itA; ++itA) {
    const unsigned int blockRowStart = basisController->extendedIndex(itA.row());
    const unsigned int nRows = basis[itA.row()]->getNContracted();
    unsigned int columnIndex = 0;
    for (Eigen::SparseVector<int>::InnerIterator itB(shellsB); itB; ++itB) {
      const unsigned int blockColStart = basisController->extendedIndex(itB.row());
      const unsigned int nCols = basis[itB.row()]->getNContracted();
      for_spin(toReturn, matrix) {
        toReturn_spin.block(rowIndex, columnIndex, nRows, nCols) =
            matrix_spin.block(blockRowStart, blockColStart, nRows, nCols);
      };
      columnIndex += nCols;
    } // itB
    rowIndex += nRows;
  } // itA
  return toReturn;
}

template<Options::SCF_MODES SCFMode>
void SystemSplittingTools<SCFMode>::diagonalizationInNonRedundantPAOBasis(
    const Eigen::MatrixXd& R_ij, const Eigen::MatrixXd& s_AO, const Eigen::MatrixXd& f_AO,
    double paoOrthogonalizationThreshold, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& transformation) {
  // Calculate the metric of the PAOs.
  Eigen::MatrixXd paoOverlapMatrix = R_ij.transpose() * s_AO * R_ij;
  // Remove linear dependencies.
  // Diagonalize overlap matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(paoOverlapMatrix);
  auto U = es.eigenvectors();
  auto s = es.eigenvalues();
  int n = s.rows();
  Eigen::VectorXd sigma(n);
  sigma.setZero();
  Eigen::VectorXd sigmaInvers(n);
  sigmaInvers.setZero();
  unsigned int notZero = 0;
  for (int i = n - 1; i >= 0; --i) {
    if (s(i) > paoOrthogonalizationThreshold) {
      sigma[i] = s[i];
      sigmaInvers[i] = 1 / s[i];
      ++notZero;
    }
    else {
      i = 0;
    }
  }
  if (notZero == 0)
    throw SerenityError("Redundant diagonalization failed!");
  // Transform the Fock matrix into this basis and diagonalize it.
  Eigen::MatrixXd orthogonalizedPAOs = (U * sigmaInvers.array().sqrt().matrix().asDiagonal()).rightCols(notZero).eval();
  Eigen::MatrixXd linearIndependentPAOs = R_ij * orthogonalizedPAOs;
  Eigen::MatrixXd f_pao = (linearIndependentPAOs.transpose() * f_AO * linearIndependentPAOs).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ef(f_pao);
  // Save eigenvalues
  eigenvalues = ef.eigenvalues();
  // Save transformation to linear independent PAO-pair basis that diagonalize the fock matrix (in this block).
  transformation = (orthogonalizedPAOs * ef.eigenvectors()).eval();
}

template<Options::SCF_MODES SCFMode>
void SystemSplittingTools<SCFMode>::diagonalizationInNonRedundantPAOBasis(const Eigen::MatrixXd& s_PAO,
                                                                          const Eigen::MatrixXd& f_PAO,
                                                                          double paoOrthogonalizationThreshold,
                                                                          Eigen::VectorXd& eigenvalues,
                                                                          Eigen::MatrixXd& transformation) {
  // Remove linear dependencies.
  // Diagonalize overlap matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s_PAO);
  auto U = es.eigenvectors();
  auto s = es.eigenvalues();
  int n = s.rows();
  Eigen::VectorXd sigma(n);
  sigma.setZero();
  Eigen::VectorXd sigmaInvers(n);
  sigmaInvers.setZero();
  unsigned int notZero = 0;
  for (int i = n - 1; i >= 0; --i) {
    if (s(i) > paoOrthogonalizationThreshold) {
      sigma[i] = s[i];
      sigmaInvers[i] = 1 / s[i];
      ++notZero;
    }
    else {
      i = 0;
    }
  }
  if (notZero == 0)
    throw SerenityError("Redundant diagonalization failed!");
  // Transform the Fock matrix into this basis and diagonalize it.
  Eigen::MatrixXd orthogonalizedPAOs = (U * sigmaInvers.array().sqrt().matrix().asDiagonal()).rightCols(notZero).eval();
  Eigen::MatrixXd f_pao = (orthogonalizedPAOs.transpose() * f_PAO * orthogonalizedPAOs).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ef(f_pao);
  // Save eigenvalues
  eigenvalues = ef.eigenvalues();
  // Save transformation to linear independent PAO-pair basis that diagonalize the fock matrix (in this block).
  transformation = (orthogonalizedPAOs * ef.eigenvectors()).eval();
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<Eigen::MatrixXd>>
SystemSplittingTools<SCFMode>::getProjectedSubsystems(std::shared_ptr<SystemController> activeSystem,
                                                      std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                      double truncationThreshold) {
  std::vector<std::shared_ptr<Eigen::MatrixXd>> overlapMatrices;
  auto basisContA = activeSystem->getBasisController();
  for (auto env : environmentSystems) {
    auto& libint = Libint::getInstance();
    auto s_AB = std::make_shared<Eigen::MatrixXd>(
        libint.compute1eInts(LIBINT_OPERATOR::overlap, env->getBasisController(), basisContA));
    double totalInterSystemBasisOverlap = 0.0;
    for (unsigned int k = 0; k < s_AB->rows(); ++k) {
      for (unsigned int l = 0; l < s_AB->cols(); ++l) {
        totalInterSystemBasisOverlap += std::abs((*s_AB)(k, l));
      } // for l
    }   // for k
    if (totalInterSystemBasisOverlap <= truncationThreshold) {
      // If the system does not get projected, forget the basis overlap.
      s_AB = nullptr;
    }
    overlapMatrices.push_back(s_AB);
  } // for env
  return overlapMatrices;
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>>
SystemSplittingTools<SCFMode>::getEnvironmentDensityControllers(std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                                                bool topDown) {
  // Build the density matrix controllers of the environment systems in the same
  // spin-polarization as the active system.
  // These are needed for the Coulomb and exchange contribution to the supersystem fock operator of systems,
  // which are not in the system pair.
  std::vector<std::shared_ptr<DensityMatrixController<SCFMode>>> environmentDensityControllers;
  for (auto sys : environmentSystems) {
    if (sys->getSCFMode() == SCFMode || topDown) {
      assert(sys->hasElectronicStructure<SCFMode>());
      environmentDensityControllers.push_back(sys->template getElectronicStructure<SCFMode>()->getDensityMatrixController());
    }
    else {
      if (sys->getSCFMode() == Options::SCF_MODES::RESTRICTED) {
        // Build unrestricted DensityMatrixController
        DensityMatrix<SCFMode> uDensMat(sys->getBasisController());
        for_spin(uDensMat) {
          uDensMat_spin = 0.5 * sys->template getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getDensityMatrix();
        };
        environmentDensityControllers.push_back(std::make_shared<DensityMatrixController<SCFMode>>(uDensMat));
      }
      else if (sys->getSCFMode() == Options::SCF_MODES::UNRESTRICTED) {
        // Build restricted DensityMatrixController
        DensityMatrix<SCFMode> rDensMat(sys->getBasisController());
        for_spin(rDensMat) {
          rDensMat_spin =
              sys->template getElectronicStructure<Options::SCF_MODES::UNRESTRICTED>()->getDensityMatrix().total();
        };
        environmentDensityControllers.push_back(std::make_shared<DensityMatrixController<SCFMode>>(rDensMat));
      }
      else {
        assert(false);
      }
    }
  } // for sys
  return environmentDensityControllers;
}

template<Options::SCF_MODES SCFMode>
void SystemSplittingTools<SCFMode>::splitSupersystemBasedOnAssignment(std::shared_ptr<SystemController> supersystem,
                                                                      std::vector<std::shared_ptr<SystemController>> fragments,
                                                                      const SpinPolarizedData<SCFMode, Eigen::VectorXi>& assignment) {
  unsigned int nFragments = fragments.size();
  const auto nOccSuper = supersystem->getNOccupiedOrbitals<SCFMode>();
  for (unsigned int iFrag = 0; iFrag < nFragments; ++iFrag) {
    auto subsystem = fragments[iFrag];
    // Generate the final orbital selection based on the criteria constructed above.
    SpinPolarizedData<SCFMode, std::vector<bool>> orbitalSelection;
    for_spin(nOccSuper, orbitalSelection, assignment) {
      orbitalSelection_spin = std::vector<bool>(nOccSuper_spin, false);
      for (unsigned int iOrb = 0; iOrb < nOccSuper_spin; ++iOrb) {
        if ((int)iFrag == assignment_spin(iOrb))
          orbitalSelection_spin[iOrb] = true;
      } // for iOrb
    };
    // Update basis, electronic structure and geometry and save them on disk.
    subsystem->getGeometry()->addAsDummy(*supersystem->getGeometry());
    subsystem->getGeometry()->deleteIdenticalAtoms();
    subsystem->setBasisController(nullptr);
    subsystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_COULOMB);
    subsystem->setBasisController(nullptr, Options::BASIS_PURPOSES::AUX_CORREL);

    // Select the part of the electronic structure.
    auto subsystemES = SystemSplittingTools<SCFMode>::splitElectronicStructure(supersystem, orbitalSelection).first;

    if (subsystem->getAtomCenteredBasisController()->getBasisLabel() !=
        supersystem->getAtomCenteredBasisController()->getBasisLabel()) {
      OutputControl::dOut << "NOTE: The basis label of the subsystem " << subsystem->getSystemName() << std::endl;
      OutputControl::dOut << "      is different from the basis label of the supersystem!" << std::endl;
      OutputControl::dOut << "      The density matrix will be projected into the subsystem basis." << std::endl;
      OutputControl::dOut << "      Orbitals will not be available until the next Fock-matrix diagonalization." << std::endl;
      OutputControl::dOut << "      Note that the projection may be inaccurate." << std::endl;
      const DensityMatrix<SCFMode>& selectedDensity = subsystemES->getDensityMatrix();
      DensityMatrix<SCFMode> projectedDensity =
          SystemSplittingTools<SCFMode>::projectMatrixIntoNewBasis(selectedDensity, subsystem->getBasisController());
      subsystemES = std::make_shared<ElectronicStructure<SCFMode>>(subsystem->getOneElectronIntegralController(),
                                                                   subsystemES->getNOccupiedOrbitals(),
                                                                   subsystemES->getMolecularOrbitals()->getNCoreOrbitals());
      subsystemES->getDensityMatrixController()->setDensityMatrix(projectedDensity);
    }
    else {
      // The same basis set is used for the subsystem and the supersystem.
      // The coefficients may need to be sorted, but no projection is necessary!
      subsystemES = SystemSplittingTools<SCFMode>::resortBasisSetOfElectronicStructure(
          subsystemES, subsystem->getBasisController(), subsystem->getOneElectronIntegralController());
    }

    subsystem->setSCFMode(SCFMode);
    subsystem->setElectronicStructure<SCFMode>(subsystemES);
    subsystem->getGeometry()->printToFile(subsystem->getHDF5BaseName(), subsystem->getSettings().identifier);
    subsystemES->toHDF5(subsystem->getHDF5BaseName(), subsystem->getSettings().identifier);
    // Calculate charges and spin of the system, since it may have changed during the partitioning.
    auto nElSubSpin = SystemSplittingTools<SCFMode>::getNElectrons(orbitalSelection).first;
    int totEffCharge = subsystem->getGeometry()->getTotalEffectiveCharge();
    unsigned int nElSub = 0;
    for_spin(nElSubSpin) {
      nElSub += nElSubSpin_spin;
    };
    int subSpin = SystemSplittingTools<SCFMode>::getSpin(nElSubSpin);
    subsystem->setSpin(subSpin);
    subsystem->setCharge(totEffCharge - nElSub);
    // Print information about the system selection to the output.
    OutputControl::nOut << "--------------------------------------------------------" << std::endl;
    OutputControl::nOut << "  Partitioning for subsystem " << subsystem->getSystemName() << std::endl;
    OutputControl::nOut << "  Charge: " << subsystem->getCharge() << std::endl;
    OutputControl::nOut << "  Spin:   " << subsystem->getSpin() << std::endl;
    OutputControl::vOut << "Orbital selection:" << std::endl;
    std::vector<std::string> preafixes = {"Alpha", "Beta"};
    unsigned int counter = 0;
    for_spin(orbitalSelection) {
      if (SCFMode == UNRESTRICTED)
        OutputControl::vOut << preafixes[counter] << " orbitals:" << std::endl;
      unsigned int printCounter = 0;
      for (unsigned int iOrb = 0; iOrb < orbitalSelection_spin.size(); ++iOrb) {
        if (orbitalSelection_spin[iOrb]) {
          OutputControl::vOut << iOrb << " ";
          ++printCounter;
          if (printCounter % 15 == 0)
            OutputControl::vOut << std::endl;
        } // if orbitalSelection_spin[iOrb]
      }   // for iOrb
      if (printCounter % 15 != 0)
        OutputControl::vOut << std::endl;
      OutputControl::nOut << "Number of occ. orbitals: " << printCounter << std::endl;
      ++counter;
    };
  } // for iSub
  OutputControl::nOut << "--------------------------------------------------------" << std::endl;
}

template class SystemSplittingTools<Options::SCF_MODES::RESTRICTED>;
template class SystemSplittingTools<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
