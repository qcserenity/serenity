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
#include "analysis/PAOSelection/TNOConstructor.h"                   //TNO construction.
#include "data/OrbitalTriple.h"                                     //Orbital triple definition.
#include "data/matrices/MatrixInBasis.h"                            //Container for the aux. metric.
#include "integrals/MO3CenterIntegralController.h"                  //Three center MO integrals.
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h" //Calculation of aux. metric.
#include "io/FormattedOutputStream.h"                               //Filtered output.
#include "misc/Timing.h"                                            //Timinigs.
#include "postHF/LocalCorrelation/LocalCorrelationController.h"     //Definition of the local correlation controller.
#include "system/SystemController.h"                                //Aux. basis controller.
/* Include Std and External Headers */
#include <Eigen/Dense>      //Dummy input for integral getter.
#include <Eigen/SparseCore> //Dummy input for integral getter.

namespace Serenity {

double DLPNO_CCSD_T0::calculateEnergyCorrection(std::shared_ptr<LocalCorrelationController> localCorrelationController) {
  std::vector<std::shared_ptr<OrbitalTriple>> triples = localCorrelationController->getOrbitalTriples();

  auto activeSystem = localCorrelationController->getActiveSystemController();
  auto auxBasisController = activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL);
  auto tnoConstructor = localCorrelationController->produceTNOConstructor();
  auto mo3CenterIntegralController = localCorrelationController->getMO3CenterIntegralController();
  // Make sure that the integrals are available in order to avoid any parallel construction of
  // the integral sets, since these functions are not thread save.
  // TODO use locks in integral calculation.
  Eigen::SparseVector<int> auxSuperDomain = Eigen::VectorXi::Constant(auxBasisController->getBasis().size(), 1).sparseView();
  mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ia_K, auxSuperDomain);
  mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::ab_K, auxSuperDomain);
  mo3CenterIntegralController->getMO3CenterInts(MO3CENTER_INTS::kl_K, auxSuperDomain);
  // Make sure that the prescreening maps exist in order to prevent multiple threads to build
  // the same map simultaneously.
  // TODO use locks in map construction.
  mo3CenterIntegralController->getSparseMapsController()->getTripletOccToAuxShellMap(true);
  mo3CenterIntegralController->getSparseMapsController()->getTripletOccToAuxShellMap(false);

  // Pre-calculate the coulomb metric.
  MatrixInBasis<RESTRICTED> metric(auxBasisController);
  Ao2MoExchangeIntegralTransformer::calculateTwoCenterIntegrals(metric);
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

  takeTime("Semi-Canonical Triples Correction");
  std::vector<double> energies = {0.0};
#ifdef _OPENMP
  unsigned int nThreads = omp_get_max_threads();
  for (unsigned int thread = 1; thread < nThreads; ++thread)
    energies.push_back(0.0);
  Eigen::setNbThreads(1);
#endif
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iT = 0; iT < triples.size(); ++iT) {
#ifdef _OPENMP
    const unsigned int threadId = omp_get_thread_num();
#else
    const unsigned int threadId = 0;
#endif
    std::shared_ptr<OrbitalTriple>& triple = triples[iT];
    /*
     * The actual calculations happen in these three lines:
     * 1. The TNO basis is constructed.
     * 2. The integrals are calculated in this basis.
     * 3. The triple energy increment is calculated.
     */
    tnoConstructor->transformToTNOBasis(triple);
    triple->calculateIntegrals(auxBasisController, mo3CenterIntegralController, metric);
    energies[threadId] += triple->getTripleEnergy();
  } // for iT
  double finalEnergy = energies[0];
#ifdef _OPENMP
  Eigen::setNbThreads(0);
  for (unsigned int thread = 1; thread < nThreads; ++thread)
    finalEnergy += energies[thread];
#endif
  timeTaken(1, "Semi-Canonical Triples Correction");
  return finalEnergy;
}

} /* namespace Serenity */
