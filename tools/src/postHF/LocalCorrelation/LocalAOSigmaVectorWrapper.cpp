/**
 * @file LocalAOSigmaVectorWrapper.cpp
 *
 * @author Moritz Bensberg
 * @date Oct 21, 2019
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
#include "postHF/LocalCorrelation/LocalAOSigmaVectorWrapper.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"                            //Tests only.
#include "data/OrbitalPair.h"                                  //PAO domains of diagonal orbital pairs.
#include "data/PAOController.h"                                //Density matrix construction in AO basis.
#include "data/SingleSubstitution.h"                           //Singles definition.
#include "postHF/LRSCF/LRSCFController.h"                      //Sigma vector construction.
#include "postHF/LRSCF/Sigmavectors/ExchangeSigmavector.h"     //Sigma vector construction.
#include "postHF/LRSCF/Sigmavectors/RI/RICoulombSigmavector.h" //Sigma vector construction.
#include "settings/Settings.h"                                 //Settings defintion.
#include "system/SystemController.h"                           //Number of occupied orbitals and tests.
#include "tasks/LRSCFTask.h"                                   //LRSCFTaskSettings for the LRSCFController.

namespace Serenity {

MatrixInBasis<RESTRICTED>
LocalAOSigmaVectorWrapper::calculatePerturbedDensityMatrix(std::shared_ptr<SystemController> system,
                                                           std::vector<std::shared_ptr<SingleSubstitution>> singles,
                                                           std::shared_ptr<PAOController> paoController, bool testRun) {
  assert(singles.size() > 0);
  unsigned int nOcc = system->getNOccupiedOrbitals<RESTRICTED>();
  const Eigen::MatrixXd occCoefficients = system->getActiveOrbitalController<RESTRICTED>()->getCoefficients().leftCols(nOcc);
  std::vector<MatrixInBasis<RESTRICTED>> ps = {};
  unsigned int nThreads = 1;
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
  Eigen::setNbThreads(1);
#endif
  for (unsigned int i = 0; i < nThreads; ++i) {
    MatrixInBasis<RESTRICTED> pIncrement(system->getBasisController());
    pIncrement.setZero();
    ps.push_back(pIncrement);
  } // for i
#pragma omp parallel for schedule(static)
  for (unsigned int iSingle = 0; iSingle < singles.size(); ++iSingle) {
#ifdef _OPENMP
    const unsigned int threadId = omp_get_thread_num();
#else
    const unsigned int threadId = 0;
#endif
    const auto& single = singles[iSingle];
    auto& pIncrement = ps[threadId];
    Eigen::MatrixXd virtCoeff = (Eigen::MatrixXd)paoController->getAllPAOs() *
                                single->getDiagonalPair()->domainProjection.transpose() * single->toPAODomain;
    if (testRun)
      virtCoeff = system->getActiveOrbitalController<RESTRICTED>()->getCoefficients().rightCols(single->t_i.rows());
    const Eigen::MatrixXd pno_x_t = virtCoeff * single->t_i;
    pIncrement += (occCoefficients.col(single->i) * pno_x_t.transpose()).transpose();
  } // for iSingle
  MatrixInBasis<RESTRICTED> p = ps[0];
  for (unsigned int i = 1; i < nThreads; ++i) {
    p += ps[i];
  }
  return p;
}

MatrixInBasis<RESTRICTED>
LocalAOSigmaVectorWrapper::getSigmaVector_AO_AO(std::shared_ptr<SystemController> system,
                                                std::vector<std::shared_ptr<SingleSubstitution>> singles,
                                                std::shared_ptr<PAOController> paoController, bool testRun) {
  LRSCFTaskSettings lrscfSettings;
  auto lrscf = std::make_shared<LRSCFController<RESTRICTED>>(system, lrscfSettings);
  std::vector<std::shared_ptr<LRSCFController<RESTRICTED>>> lrscfs = {lrscf};
  RICoulombSigmavector<RESTRICTED> coulombSigmaVector(lrscfs);
  ExchangeSigmavector<RESTRICTED> exchangeSigmaVector(lrscfs);
  const MatrixInBasis<RESTRICTED> p = calculatePerturbedDensityMatrix(system, singles, paoController, testRun);
  std::vector<std::vector<MatrixInBasis<RESTRICTED>>> p_set = {{p}};
  // Note: The way the sigma vectors are defined makes this here quite awkward ...
  auto p_ptr1 = std::make_unique<std::vector<std::vector<MatrixInBasis<RESTRICTED>>>>(p_set);
  auto p_ptr2 = std::make_unique<std::vector<std::vector<MatrixInBasis<RESTRICTED>>>>(p_set);

  // I split the construction and calculation of the sigma vector up, because the compliler becomes confused otherwise.
  // It would try to cast the MatrixInBasis to an Eigen::MatrixXd and then realize that it can not reconstruct a
  // MatrixInBasis from an Eigen::MatrixXd ...
  MatrixInBasis<RESTRICTED> sigma_AO_AO(system->getBasisController());
  sigma_AO_AO = 2.0 * (*coulombSigmaVector.calcF(0, 0, std::move(p_ptr1)))[0][0].transpose();
  sigma_AO_AO -= (*exchangeSigmaVector.calcF(0, 0, std::move(p_ptr2)))[0][0].transpose();
  return sigma_AO_AO;
}

} /* namespace Serenity */
