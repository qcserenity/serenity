/**
 * @file DLPNO_CCSD_T0.cpp
 *
 * @author Moritz Bensberg
 * @date Nov 5, 2019
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
#include "postHF/CC/DLPNO_CCSD_T0.h"
/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/PNOConstructor.h"                   //PNO construction for distant orbital pairs.
#include "analysis/PAOSelection/TNOConstructor.h"                   //TNO construction.
#include "basis/Basis.h"                                            //Aux basis size.
#include "data/OrbitalTriple.h"                                     //Orbital triple definition.
#include "data/OrbitalTripleSet.h"                                  //Set-wise representation of the orbital triples
#include "data/matrices/MatrixInBasis.h"                            //Container for the aux. metric.
#include "integrals/MO3CenterIntegralController.h"                  //Three center MO integrals.
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h" //Calculation of aux. metric.
#include "io/FormattedOutputStream.h"                               //Filtered output.
#include "misc/Timing.h"                                            //Timings.
#include "postHF/LocalCorrelation/LocalCorrelationController.h"     //Definition of the local correlation controller.
#include "system/SystemController.h"                                //Aux. basis controller.
/* Include Std and External Headers */
#include <Eigen/SparseCore> // Input for integral getter.

namespace Serenity {

double DLPNO_CCSD_T0::calculateEnergyCorrection(std::shared_ptr<LocalCorrelationController> localCorrelationController) {
  takeTime("Semi-Canonical Triples Correction");
  // Generate PNO basis for distant orbital pairs in case this was not done already. If it was, this will change nothing.
  auto pnoConstructor = localCorrelationController->producePNOConstructor();
  auto distantOrbitalPairs = localCorrelationController->getOrbitalPairs(OrbitalPairTypes::DISTANT_TRIPLES);
  pnoConstructor->transformToPNOBasis(distantOrbitalPairs);
  // Update the orbital pair lists.
  localCorrelationController->selectDistantOrbitalPairs();
  // Get the orbital triples and print some basic information.
  auto orbitalTripleSets = localCorrelationController->getOrbitalTripleSets();
  auto activeSystem = localCorrelationController->getActiveSystemController();
  auto auxBasisController = activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
  auto tnoConstructor = localCorrelationController->produceTNOConstructor();
  auto mo3CenterIntegralController = localCorrelationController->getTriplesMO3CenterIntegralController();
  auto triples = localCorrelationController->getOrbitalTriples();
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  OutputControl::mOut << " Semi-Canonical Triples Correction" << std::endl;
  OutputControl::nOut << "  Number of close pairs      "
                      << localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE).size() << std::endl;
  OutputControl::nOut << "  Number of distant pairs    "
                      << localCorrelationController->getOrbitalPairs(OrbitalPairTypes::DISTANT_TRIPLES).size() << std::endl;
  OutputControl::nOut << "  Number of triples          " << triples.size() << std::endl;
  OutputControl::nOut << "  TNO threshold              " << localCorrelationController->getSettings().tnoThreshold
                      << std::endl;
  OutputControl::nOut << "  MKN strong scaling         "
                      << localCorrelationController->getSettings().crudeStrongTripFactor << std::endl;
  OutputControl::nOut << "  MKN weak scaling           "
                      << localCorrelationController->getSettings().crudeWeakTripFactor << std::endl;
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  // Pre-calculate the coulomb metric.
  MatrixInBasis<RESTRICTED> metric(auxBasisController);
  Ao2MoExchangeIntegralTransformer::calculateTwoCenterIntegrals(metric);
  double finalEnergy = 0.0;
  double averageNTNOs = 0;
  double averageNAux = 0;
  for (auto triplesSet : orbitalTripleSets) {
    // Make sure that the integrals are available in order to avoid any parallel construction of
    // the integral sets, since these functions are not thread save.
    // TODO use locks in integral calculation.
    const Eigen::SparseVector<int> auxSuperDomain = triplesSet->getTotalFittingDomain();
    const MO3CenterIntegrals& iaK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain);
    const MO3CenterIntegrals& abK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ab_K, auxSuperDomain);
    const MO3CenterIntegrals& klK = mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::kl_K, auxSuperDomain);
    unsigned int nThreads = 1;
#ifdef _OPENMP
    nThreads = omp_get_max_threads();
    Eigen::setNbThreads(1);
#endif

    std::vector<double> energies(nThreads, 0.0);
    std::vector<unsigned int> nTNOsTotal(nThreads, 0);
    std::vector<unsigned int> nAuxTotal(nThreads, 0);
    bool takeTiming = (nThreads == 1) ? true : false;
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iT = 0; iT < triplesSet->size(); ++iT) {
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      std::shared_ptr<OrbitalTriple>& triple = (*triplesSet)[iT];
      /*
       * The actual calculations happen in these three lines:
       * 1. The TNO basis is constructed.
       * 2. The integrals are calculated in this basis.
       * 3. The triple energy increment is calculated.
       */
      if (takeTiming)
        Timings::takeTime("TNO construction");
      tnoConstructor->transformToTNOBasis(triple);
      if (takeTiming)
        Timings::timeTaken("TNO construction");
      if (takeTiming)
        Timings::takeTime("Triples Integrals");
      triple->calculateIntegrals(auxBasisController, mo3CenterIntegralController, iaK, abK, klK, metric);
      nTNOsTotal[threadId] += triple->getNTNOs();
      nAuxTotal[threadId] += triple->getNLocalAuxiliaryFunctions();
      if (takeTiming)
        Timings::timeTaken("Triples Integrals");
      if (takeTiming)
        Timings::takeTime("Triples Energy Increment");
      energies[threadId] += triple->getTripleEnergy();
      if (takeTiming)
        Timings::timeTaken("Triples Energy Increment");
    } // for iT
    finalEnergy += energies[0];
    averageNTNOs += nTNOsTotal[0];
    averageNAux += nAuxTotal[0];
#ifdef _OPENMP
    Eigen::setNbThreads(0);
    for (unsigned int thread = 1; thread < nThreads; ++thread) {
      finalEnergy += energies[thread];
      averageNTNOs += nTNOsTotal[thread];
      averageNAux += nAuxTotal[thread];
    }
#endif
  } // for orbitalTripleSet
  averageNTNOs = (triples.size()) ? averageNTNOs / triples.size() : 0;
  averageNAux = (triples.size()) ? averageNAux / triples.size() : 0;

  OutputControl::nOut << "  Average number of TNOs:                   " << averageNTNOs << std::endl;
  OutputControl::nOut << "  Average number of aux. functions:         " << averageNAux << std::endl;
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  timeTaken(1, "Semi-Canonical Triples Correction");

  // TMP
  for (const auto& triple : triples) {
    if (std::abs(triple->getTripleEnergy()) > 1e-2) {
      std::cout << triple->getI() << " " << triple->getJ() << " " << triple->getK() << "  " << triple->getTripleEnergy()
                << std::endl;
    }
  }

  return finalEnergy;
}

} /* namespace Serenity */
