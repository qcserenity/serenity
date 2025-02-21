/**
 * @file   FCIDumpFileWriter.cpp
 *
 * @date   Jan 30, 2024
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
#include "io/FCIDumpFileWriter.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"                        //Access to density and coefficient matrix.
#include "data/ElectronicStructure.h"                        // Get occupations from density matrix controller.
#include "data/OrbitalController.h"                          //Access to the coefficient matrix.
#include "data/OrbitalPair.h"                                //Pseudo orbital pairs/quartets for integral prescreening.
#include "data/SparseMapsController.h"                       //Integral prescreening.
#include "data/grid/BasisFunctionOnGridControllerFactory.h"  //DOIs.
#include "data/grid/DifferentialOverlapIntegralCalculator.h" //Prescreening.
#include "data/matrices/DensityMatrixController.h"           //Orbital occupations.
#include "data/matrices/MatrixInBasis.h"                     //Data type definition.
#include "integrals/MO3CenterIntegralController.h"           //Integral transformation and evaluation.
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h" //Integral transformation.
#include "io/FormattedOutputStream.h"                               //Filtered output.
#include "misc/SerenityError.h"                                     //Errors.
#include "misc/SystemSplittingTools.h"                              //Extract local metric.
#include "misc/Timing.h"                                            //Timings.
#include "misc/WarningTracker.h"
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h" //Pseudo orbital pairs/quartets for integral prescreening.
#include "potentials/ERIPotential.h" // Core one particle matrix excluding the active space occupied orbitals.
#include "settings/Settings.h"       //Access to the system settings.
#include "system/SystemController.h" //Access to system properties.
#include "tasks/FDETask.h"           //Core one particle matrix and core energy evaluation.
/* Include Std and External Headers */
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
FCIDumpFileWriter<SCFMode>::FCIDumpFileWriter(std::shared_ptr<SystemController> activeSystem,
                                              std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                              const FCIDumpFileWriterTaskSettings& settings,
                                              std::shared_ptr<MatrixInBasis<SCFMode>> oneParticleIntegrals)
  : _settings(settings), _activeSystem(activeSystem), _environmentSystems(environmentSystems) {
  if (oneParticleIntegrals) {
    _oneParticleIntegrals = std::make_unique<MatrixInBasis<SCFMode>>(*oneParticleIntegrals);
  }
  if (activeSystem->getSettings().method != Options::ELECTRONIC_STRUCTURE_THEORIES::HF) {
    throw SerenityError("The active system must be described with HF!");
  }
}

template<Options::SCF_MODES SCFMode>
void FCIDumpFileWriter<SCFMode>::writeFCIDumpFile(std::ostream& outputStream,
                                                  const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  writeHeader(outputStream, orbitalRange);
  writeERIPart(outputStream, orbitalRange);
  writeCorePart(outputStream, orbitalRange);
  writeCoreEnergy(outputStream);
  OutputControl::nOut << "  FCI dump file successfully written to output stream!" << std::endl;
}
template<Options::SCF_MODES SCFMode>
void FCIDumpFileWriter<SCFMode>::writeCorePart(std::ostream& outputStream,
                                               const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  /*
   * Example output
   *  -1.9733392198594029                 1           1           0           0
   *  -0.18208843465986255                2           1           0           0
   */
  if (SCFMode != RESTRICTED) {
    throw SerenityError("Writing FCI dump files for unrestricted orbitals is not supported at the moment!");
  }
  const auto hcore = getCore(orbitalRange);
  OutputControl::nOut << "  Writing core hamiltonian                               ...";
  OutputControl::nOut.flush();
  unsigned int counter = 1;
  for_spin(orbitalRange, hcore) {
    const unsigned int nOrbs = orbitalRange_spin.size();
    for (unsigned int i = 0; i < nOrbs; ++i) {
      const unsigned int iRange = orbitalRange_spin[i];
      for (unsigned int j = 0; j <= i; ++j) {
        const unsigned int jRange = orbitalRange_spin[j];
        outputStream << std::setw(18) << hcore_spin(iRange, jRange);
        outputStream << std::setw(6) << i + counter;
        outputStream << std::setw(6) << j + counter;
        outputStream << std::setw(6) << 0 << std::setw(6) << 0 << "\n";
      }
    }
    counter += nOrbs;
  };
  outputStream.flush();
  OutputControl::nOut << " done" << std::endl;
}
template<Options::SCF_MODES SCFMode>
void FCIDumpFileWriter<SCFMode>::writeHeader(std::ostream& outputStream,
                                             const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  /*
   * Example head:
   * &FCI NORB=             9 ,NELEC= 8,MS2=0,
   * ORBSYM=1,1,1,1,1,1,1,1,1,
   * ISYM=1,
   * &END
   */
  OutputControl::nOut.flush();
  unsigned int nOrbitals = 0;
  unsigned int nElectrons = 0;
  const auto occupations = _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getOccupations();
  const unsigned int alphaElectronExcess = _activeSystem->getSpin();
  OutputControl::nOut << "  Writing FCI Dump header                                ...";
  for_spin(orbitalRange, occupations) {
    nOrbitals += orbitalRange_spin.size();
    for (const auto& iOrb : orbitalRange_spin) {
      nElectrons += occupations_spin(iOrb);
    }
  };
  outputStream << "&FCI NORB= " << nOrbitals << ", NELEC= " << nElectrons << ", MS2= " << alphaElectronExcess << ",\n";
  outputStream << "ORBSYM=";
  for (unsigned int i = 0; i < nOrbitals; ++i) {
    outputStream << " 1,";
  }
  outputStream << "\n";
  outputStream << "ISYM= 1\n";
  outputStream << "&END" << std::endl;
  OutputControl::nOut << " done" << std::endl;
}

template<Options::SCF_MODES SCFMode>
SPMatrix<SCFMode> FCIDumpFileWriter<SCFMode>::getCore(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  const auto& aoIntegrals = getOneParticleIntegrals(orbitalRange);
  const CoefficientMatrix<SCFMode>& coefficients = _activeSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  SPMatrix<SCFMode> moIntegrals;
  for_spin(aoIntegrals, coefficients, moIntegrals) {
    moIntegrals_spin = coefficients_spin.transpose() * aoIntegrals_spin * coefficients_spin;
  };
  return moIntegrals;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode> FCIDumpFileWriter<SCFMode>::getFockMatrix() {
  auto electronicStructure = _activeSystem->template getElectronicStructure<SCFMode>();
  ;
  if (electronicStructure->checkFock()) {
    return electronicStructure->getFockMatrix();
  }
  else {
    WarningTracker::printWarning("WARNING: No Fock matrix available for the active system! Read from disk instead!", true);
    try {
      electronicStructure->fockFromHDF5(_activeSystem->getHDF5BaseName(), _activeSystem->getSystemIdentifier());
      return electronicStructure->getFockMatrix();
    }
    catch (...) {
      WarningTracker::printWarning("WARNING: No Fock matrix on disk! Reconstruct Fock matrix from scratch!", true);
    }
  } // else (electronicStructure->checkFock())
  FDETask<SCFMode> fdeTask(_activeSystem, _environmentSystems);
  fdeTask.settings.embedding = _settings.embedding;
  fdeTask.settings.skipSCF = true;
  fdeTask.run();
  return electronicStructure->getFockMatrix();
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>
FCIDumpFileWriter<SCFMode>::getOneParticleIntegrals(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  if (!_lastActiveSpace || orbitalRange != *_lastActiveSpace) {
    this->updateOneParticleIntegrals(orbitalRange);
  }
  return *_oneParticleIntegrals;
}

template<Options::SCF_MODES SCFMode>
void FCIDumpFileWriter<SCFMode>::updateOneParticleIntegrals(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  _oneParticleIntegrals = std::make_unique<FockMatrix<SCFMode>>(_activeSystem->getBasisController());
  auto& h = *_oneParticleIntegrals;
  auto& coreEnergy = *_coreEnergy;
  auto& totalEnergy = *_totalUncorrelatedEnergy;
  auto& eActSpace = *_totalInitialActiveSpaceEnergy;
  _lastActiveSpace = std::make_unique<SpinPolarizedData<SCFMode, std::vector<unsigned int>>>(orbitalRange);
  /**
   * 1. Calculate full Fock operator f.
   * 2. Calculate 2e- contribution of active occupied orbitals gOccAct.
   * 3. Result = f - gOccAct
   */
  auto f = this->getFockMatrix();
  auto molecularOrbitals = _activeSystem->getActiveOrbitalController<SCFMode>();
  auto originalOccupations = _activeSystem->getElectronicStructure<SCFMode>()->getDensityMatrixController()->getOccupations();
  auto nOrbitals = _activeSystem->getBasisController()->getNBasisFunctions();
  SpinPolarizedData<SCFMode, Eigen::VectorXd> occupations(Eigen::VectorXd::Zero(nOrbitals));
  for_spin(occupations, orbitalRange, originalOccupations) {
    for (const auto& iOrb : orbitalRange_spin) {
      occupations_spin(iOrb) = originalOccupations_spin(iOrb);
    }
  };
  OutputControl::nOut << "\n  Calculating G(active-space)                            ...";
  OutputControl::nOut.flush();
  OutputControl::dOut << std::endl;
  auto systemSettings = _activeSystem->getSettings();
  auto activeSpaceDensityMatrixController = std::make_shared<DensityMatrixController<SCFMode>>(molecularOrbitals, occupations);
  ERIPotential<SCFMode> twoElectronInteraction(
      _activeSystem, activeSpaceDensityMatrixController, 1.0, systemSettings.basis.integralThreshold,
      systemSettings.basis.integralIncrementThresholdStart, systemSettings.basis.integralIncrementThresholdEnd,
      systemSettings.basis.incrementalSteps);
  auto gOccAct = twoElectronInteraction.getMatrix();
  OutputControl::nOut << " done" << std::endl;

  h = f - gOccAct;
  totalEnergy = _activeSystem->template getElectronicStructure<SCFMode>()->getEnergy();
  eActSpace = twoElectronInteraction.getEnergy(activeSpaceDensityMatrixController->getDensityMatrix());
  auto activeSpaceDensityMatrix = activeSpaceDensityMatrixController->getDensityMatrix();
  for_spin(activeSpaceDensityMatrix, h) {
    eActSpace += activeSpaceDensityMatrix_spin.cwiseProduct(h_spin).sum();
  };
  coreEnergy = totalEnergy - eActSpace;
}

template<Options::SCF_MODES SCFMode>
void FCIDumpFileWriter<SCFMode>::writeCoreEnergy(std::ostream& outputStream) {
  outputStream << std::setw(18) << *_coreEnergy << std::setw(6) << 0 << std::setw(6) << 0 << std::setw(6) << 0
               << std::setw(6) << 0 << std::endl;
}

template<Options::SCF_MODES SCFMode>
void FCIDumpFileWriter<SCFMode>::writeERIPart(std::ostream& outputStream,
                                              const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  if (SCFMode != RESTRICTED) {
    /*
     * The code below needs some adjustments for unrestricted reference orbitals.
     */
    throw SerenityError("Writing FCI dump files for unrestricted orbitals is not supported at the moment!");
  }
  /*
   * The following code may look a bit cryptic. The basic idea is this:
   * 1. Create sparse mappings (SparseMapsController):
   *   * orbitals -> auxiliary functions (K); based on Mulliken population analysis mapping orbitals to atoms.
   *   * orbitals -> basis function shells; based on prescreening of the absolute coefficient value.
   *   * extend the maps above to make sure that if we calculate an integral (ik|jl), all i(extended) = i u k u j u l.
   *   To achieve this, we provide the sparse maps controller with the orbital coefficients directly, essentially
   * telling it that all these orbitals are "occupied" orbitals, at least for the purpose of the integral calculation.
   *
   * 2. Calculate the integrals (ik|K) over auxiliary functions and transform them to the MO basis. This is handled by
   * the MO3CenterIntegralController.
   * 3. Extract the integrals (ik|K) and (jl|Q) from the list and calculate the full integral (ik|jl).
   * 4. Write the integrals to the output file.
   */
  Timings::takeTime("Local Cor. -   Int. Transform.");
  outputStream << std::scientific << std::setprecision(10);
  const auto allCoefficients = _activeSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  auto auxiliaryBasisController = _activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
  for_spin(allCoefficients, orbitalRange) {
    const unsigned int nOrbs = orbitalRange_spin.size();
    auto selectedCoefficients = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(allCoefficients_spin.rows(), nOrbs));
    unsigned int counter = 0;
    for (const auto& iOrb : orbitalRange_spin) {
      selectedCoefficients->col(counter) = allCoefficients_spin.col(iOrb).eval();
      counter++;
    }
    auto orb2OrbMap = std::make_shared<SparseMap>(getOrbital2OrbitalMap(selectedCoefficients));
    auto orbitalPairs = getPseudoOrbitalPairs(orb2OrbMap);
    const Eigen::VectorXd kThresholds = Eigen::VectorXd::Constant(nOrbs, _settings.mullikenThreshold);
    const Eigen::VectorXd shellThresholds = Eigen::VectorXd::Constant(nOrbs, _settings.orbitalToShellThreshold);
    auto sparseMapsController = std::make_shared<SparseMapsController>(
        _activeSystem, selectedCoefficients, selectedCoefficients, orb2OrbMap, orbitalPairs,
        std::vector<std::shared_ptr<OrbitalPair>>{}, kThresholds, shellThresholds);
    MO3CenterIntegralController mo3CenterIntegralController(
        auxiliaryBasisController, _activeSystem->getBasisController(), sparseMapsController, selectedCoefficients,
        selectedCoefficients, _activeSystem->getSystemName(), _activeSystem->getSystemIdentifier());
    const Eigen::SparseVector<int> allAuxFunctions =
        Eigen::VectorXi::Constant(
            _activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL)->getReducedNBasisFunctions(), 1)
            .sparseView();
    const MO3CenterIntegrals& klK = mo3CenterIntegralController.getMO3CenterInts(MO3CENTER_INTS::kl_K, allAuxFunctions);
    const auto occToK = mo3CenterIntegralController.getSparseMapsController()->getExtendedOccToAuxShellMap();
    MatrixInBasis<RESTRICTED> M(auxiliaryBasisController);
    Ao2MoExchangeIntegralTransformer::calculateTwoCenterIntegrals(M);
    const auto& reducedOccIndices = mo3CenterIntegralController.getIndices(ORBITAL_TYPE::OCCUPIED);

    OutputControl::nOut << "  Writing electron repulsion integrals                   ...";
    OutputControl::nOut.flush();
    // The pair construction already asserts that i >= j
    for (const auto& pair : orbitalPairs) {
      const unsigned int i = pair->i;
      const unsigned int j = pair->j;
      std::map<unsigned int, unsigned int> kTokIndexMap; // Key: k, value: index in kIndices list.
      unsigned int kCounter = 0;
      for (const auto& kSet : pair->coupledPairs) {
        unsigned int k = kSet->getK();
        if (kTokIndexMap.find(k) == kTokIndexMap.end()) {
          kTokIndexMap.insert(std::make_pair(k, kCounter));
          ++kCounter;
        }
      }

      const Eigen::SparseVector<int> pairDomainToK = (occToK.col(i) + occToK.col(j)).pruned();
      const Eigen::LLT<Eigen::MatrixXd> lltMetric =
          SystemSplittingTools<RESTRICTED>::getMatrixBlockShellWise(M, pairDomainToK, pairDomainToK).llt();
      const unsigned int nLocalAux = lltMetric.cols();
      const Eigen::MatrixXd ik_K = Ao2MoExchangeIntegralTransformer::get_ikK(
          auxiliaryBasisController, nLocalAux, pairDomainToK, reducedOccIndices, kTokIndexMap, i, klK); // nOrb x nAux
      const Eigen::MatrixXd jl_K = Ao2MoExchangeIntegralTransformer::get_ikK(
          auxiliaryBasisController, nLocalAux, pairDomainToK, reducedOccIndices, kTokIndexMap, j, klK); // nOrb x nAux
      const Eigen::MatrixXd jl_V_T = lltMetric.solve(jl_K.transpose()).eval();
      const Eigen::MatrixXd ik_jl = ik_K * jl_V_T; // nOrb x nOrb
      for (auto& kPair : kTokIndexMap) {
        const unsigned int k = kPair.first;
        const unsigned int kIndex = kPair.second;
        if (k > i) {
          break;
        }
        for (auto& lPair : kTokIndexMap) {
          const unsigned int l = lPair.first;
          const unsigned int lIndex = lPair.second;
          if (l > j) {
            break;
          }
          /*
           * Example body:
           *  0.32468447981663390               1           1           1           1
           * -3.3373939961268860E-003           2           1           1           1
           *  4.9814559309284842E-002           2           1           2           1
           */
          if (std::fabs(ik_jl(kIndex, lIndex)) < _settings.integralSizeCutOff) {
            continue;
          }
          outputStream << std::setw(18) << ik_jl(kIndex, lIndex);
          outputStream << std::setw(6) << i + 1;
          outputStream << std::setw(6) << k + 1;
          outputStream << std::setw(6) << j + 1;
          outputStream << std::setw(6) << l + 1 << "\n";
        }
      }
      outputStream.flush();
    }
    OutputControl::nOut << " done" << std::endl;
  };
  Timings::timeTaken("Local Cor. -   Int. Transform.");
}

template<Options::SCF_MODES SCFMode>
std::vector<std::shared_ptr<OrbitalPair>>
FCIDumpFileWriter<SCFMode>::getPseudoOrbitalPairs(std::shared_ptr<SparseMap> orb2OrbMap) {
  std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs;
  Eigen::MatrixXi orbitalPairIndices = Eigen::MatrixXi::Constant(orb2OrbMap->cols(), orb2OrbMap->cols(), -1);
  unsigned int counter = 0;
  // Probably, we have to take all ij combinations since they will not naturally truncate in the integrals (ik|jl).
  for (unsigned int i = 0; i < orb2OrbMap->cols(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      orbitalPairs.push_back(std::make_shared<OrbitalPair>(i, j, 1e-7, 1e-4, 1.0));
      orbitalPairIndices(i, j) = counter;
      orbitalPairIndices(j, i) = counter;
      counter++;
    }
  }
  /*
   * We will calculate the integral (ik|jl). The indices ij are independent, i.e., integrals will not vanish if ij do
   * show no differential overlap. However, we can select k and l such that they must have differential overlap with at
   * least i or j. The indices for k and l are collected in the kSets below.
   */
  for (auto& pair : orbitalPairs) {
    std::vector<std::shared_ptr<CouplingOrbitalSet>> kSets;
    const Eigen::SparseVector<int> allK = orb2OrbMap->col(pair->i) + orb2OrbMap->col(pair->j);
    for (Eigen::SparseVector<int>::InnerIterator itK(allK); itK; ++itK) {
      const unsigned int k = itK.row();
      auto ikPairIndex = orbitalPairIndices(pair->i, k);
      auto jkPairIndex = orbitalPairIndices(pair->j, k);
      if (ikPairIndex >= 0 && jkPairIndex >= 0) {
        kSets.push_back(std::make_shared<CouplingOrbitalSet>(pair, orbitalPairs[ikPairIndex], orbitalPairs[jkPairIndex], k));
      }
    }
    pair->coupledPairs = kSets;
  }
  return orbitalPairs;
}

template<Options::SCF_MODES SCFMode>
SparseMap FCIDumpFileWriter<SCFMode>::getOrbital2OrbitalMap(std::shared_ptr<Eigen::MatrixXd> coefficients) {
  const unsigned int nOrbitals = coefficients->cols();
  auto basFuncOnGridController = BasisFunctionOnGridControllerFactory::produce(
      _activeSystem->getSettings(), _activeSystem->getBasisController(), _activeSystem->getGridController());
  Eigen::SparseMatrix<int> basisFunctionsToOrbitalMap =
      SystemSplittingTools<Options::SCF_MODES::RESTRICTED>::reduceMatrixToMullikenNetPopulationMap(
          *coefficients, _settings.doiNetThreshold);
  Eigen::MatrixXd dois;
  DifferentialOverlapIntegralCalculator::calculateDOI(*coefficients, *coefficients, basisFunctionsToOrbitalMap,
                                                      basisFunctionsToOrbitalMap, basFuncOnGridController, dois);
  std::vector<Eigen::Triplet<int>> triplets;
  for (unsigned int iOrb = 0; iOrb < nOrbitals; ++iOrb) {
    for (unsigned int jOrb = 0; jOrb <= iOrb; ++jOrb) {
      if (dois(iOrb, jOrb) > _settings.doiIntegralPrescreening) {
        triplets.emplace_back(iOrb, jOrb, 1);
        triplets.emplace_back(jOrb, iOrb, 1);
      }
    }
  }
  SparseMap orb2Orb(nOrbitals, nOrbitals);
  orb2Orb.setFromTriplets(triplets.begin(), triplets.end());
  return orb2Orb;
}

template<Options::SCF_MODES SCFMode>
double FCIDumpFileWriter<SCFMode>::getTotalCoreEnergy(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  this->getOneParticleIntegrals(orbitalRange);
  return *_coreEnergy;
}
template<Options::SCF_MODES SCFMode>
double
FCIDumpFileWriter<SCFMode>::getUncorrelatedActiveSpaceEnergy(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange) {
  this->getOneParticleIntegrals(orbitalRange);
  return *_totalInitialActiveSpaceEnergy;
}
template<Options::SCF_MODES SCFMode>
double FCIDumpFileWriter<SCFMode>::getTotalUncorrelatedEnergy() {
  if (_lastActiveSpace) {
    this->getOneParticleIntegrals(*_lastActiveSpace);
  }
  else {
    this->getOneParticleIntegrals({});
  }
  return *_totalUncorrelatedEnergy;
}

template class FCIDumpFileWriter<Options::SCF_MODES::RESTRICTED>;
template class FCIDumpFileWriter<Options::SCF_MODES::UNRESTRICTED>;

} // namespace Serenity