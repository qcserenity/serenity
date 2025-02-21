/**
 * @file CC2HelperFunctions.cpp
 * @author Niklas Niemeyer
 *
 * @date Jul 25, 2023
 *
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
#include "postHF/LRSCF/RICC2/CC2HelperFunctions.h"
/* Include Serenity Internal Headers */
#include "io/FormattedOutputStream.h"
#include "io/HDF5.h"
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "postHF/LRSCF/RICC2/CC2Controller.h"
#include "postHF/LRSCF/Tools/EigenvalueSolver.h"
#include "postHF/LRSCF/Tools/NonlinearEigenvalueSolver.h"
#include "postHF/LRSCF/Tools/NonlinearResponseSolver.h"
#include "system/SystemController.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculateCISEigenvectors(LRSCFTaskSettings& settings,
                                                           std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                           Eigen::VectorXd& eigenvalues, unsigned nDimension,
                                                           const Eigen::VectorXd& diagonal, unsigned initialSubspace,
                                                           SigmaCalculator sigmaCalculator) {
  // Integral-direct sigmavectors skip diagonal blocks of the response matrix
  // for coupled calculations. In order to have the RI-CIS calculation performed
  // correctly, the global method of this task is temporarily set to TDA.
  auto old_method = settings.method;
  settings.method = Options::LR_METHOD::TDA;

  printSubSectionTitle("RI-CIS");
  double cisthresh = (old_method == Options::LR_METHOD::CISD) ? settings.conv : settings.preopt;
  EigenvalueSolver eigensolver(nDimension, settings.nEigen, diagonal, cisthresh, settings.maxCycles,
                               settings.maxSubspaceDimension, initialSubspace, Options::LR_METHOD::TDA, sigmaCalculator);

  eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(eigensolver.getEigenvectors());
  eigenvalues = eigensolver.getEigenvalues();

  settings.method = old_method;
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculateRightEigenvectors(
    LRSCFTaskSettings& settings, std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
    Eigen::VectorXd& eigenvalues, const Eigen::VectorXd& diagonal, NonlinearSigmaCalculator sigmaCalculator,
    Options::LRSCF_TYPE type, std::function<void(std::vector<Eigen::MatrixXd>&, Eigen::VectorXd&)> writeToDisk) {
  printSubSectionTitle("Right Eigenvector");
  NonlinearEigenvalueSolver nlEigenSolver(settings.nEigen, diagonal, settings.conv, settings.preopt, settings.diis,
                                          settings.diisStore,
                                          (type == Options::LRSCF_TYPE::COUPLED && !settings.fullFDEc) ? 1 : settings.maxCycles,
                                          (type == Options::LRSCF_TYPE::COUPLED) ? true : false, settings.method,
                                          (*eigenvectors)[0], eigenvalues, sigmaCalculator, writeToDisk);
  nlEigenSolver.solve();
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculateLeftEigenvectors(LRSCFTaskSettings& settings,
                                                            std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                            Eigen::VectorXd& eigenvalues, const Eigen::VectorXd& diagonal,
                                                            NonlinearSigmaCalculator sigmaCalculator,
                                                            Options::LRSCF_TYPE type) {
  Eigen::VectorXd rightEigenvalues = eigenvalues;
  printSubSectionTitle("Left Eigenvector");
  NonlinearEigenvalueSolver nlEigenSolver(settings.nEigen, diagonal, settings.conv, settings.preopt, settings.diis,
                                          settings.diisStore,
                                          (type == Options::LRSCF_TYPE::COUPLED && !settings.fullFDEc) ? 1 : settings.maxCycles,
                                          (type == Options::LRSCF_TYPE::COUPLED) ? true : false, settings.method,
                                          (*eigenvectors)[1], eigenvalues, sigmaCalculator);
  nlEigenSolver.solve();

  if ((eigenvalues - rightEigenvalues).cwiseAbs().maxCoeff() > settings.conv && type != Options::LRSCF_TYPE::COUPLED) {
    printf("  Careful: left and right eigenvalues are not converged to %5.2e Hartree!\n", settings.conv);
  }
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::normalizeBothEigenvectors(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                            const LRSCFTaskSettings& settings,
                                                            std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                            Eigen::VectorXd& eigenvalues) {
  Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(settings.nEigen, 2);
  long voStart = 0;
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
    };
    std::vector<Eigen::MatrixXd> eigenvectorsI(2);
    eigenvectorsI[0] = (*eigenvectors)[0].middleRows(voStart, nvno);
    eigenvectorsI[1] = (*eigenvectors)[1].middleRows(voStart, nvno);

    norms += ilrscf->getCC2Controller()->normalizeEigenvectors(eigenvectorsI, eigenvalues);
    voStart += nvno;
  }

  for (unsigned iEigen = 0; iEigen < settings.nEigen; ++iEigen) {
    (*eigenvectors)[0].col(iEigen) *= 1 / std::sqrt(norms(iEigen, 0));
    (*eigenvectors)[1].col(iEigen) *= std::sqrt(norms(iEigen, 0)) / norms(iEigen, 1);
  }
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::normalizeRightEigenvectors(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                             const LRSCFTaskSettings& settings,
                                                             std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                             Eigen::VectorXd& eigenvalues) {
  // Normalize only right eigenvectors.
  Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(settings.nEigen, 1);
  long voStart = 0;
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
    };

    std::vector<Eigen::MatrixXd> eigenvectorsI(1);
    eigenvectorsI[0] = (*eigenvectors)[0].middleRows(voStart, nvno);

    norms += ilrscf->getCC2Controller()->normalizeEigenvectors(eigenvectorsI, eigenvalues);
    voStart += nvno;
  }
  for (unsigned iEigen = 0; iEigen < settings.nEigen; ++iEigen) {
    (*eigenvectors)[0].col(iEigen) *= 1 / std::sqrt(norms(iEigen, 0));
    (*eigenvectors)[1].col(iEigen) *= 1 / std::sqrt(norms(iEigen, 0));
  }
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::prepareVectors(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                 const LRSCFTaskSettings& settings,
                                                 std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                 std::shared_ptr<std::vector<Eigen::MatrixXd>>& transitiondensities) {
  // Get dimensions.
  long nvnospace = 0, momospace = 0;
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    for_spin(nv, no) {
      nvnospace += nv_spin * no_spin;
      momospace += (nv_spin + no_spin) * (nv_spin + no_spin);
    };
  }

  // From this point on it is assumed that eigenvectors[0] are the right, and eigenvectors[1] are
  // the left eigenvectors[1]. If no eigenvectors have been determined yet, just fill up with null
  // matrices because they will not be used further (this can only occur for response property
  // calculations where no eigenvalues have been determined beforehand).
  if (!eigenvectors) {
    eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd>>();
  }
  for (unsigned iSize = eigenvectors->size(); iSize < 2; ++iSize) {
    eigenvectors->push_back(Eigen::MatrixXd::Zero(nvnospace, settings.nEigen));
  }

  // eigenvectors[2] : State multiplier (col 0 = ground state).
  eigenvectors->push_back(Eigen::MatrixXd::Zero(nvnospace, settings.ccexdens ? (1 + settings.nEigen) : 1));
  // eigenvectors[3] : Transition multiplier.
  eigenvectors->push_back(Eigen::MatrixXd::Zero(nvnospace, settings.nEigen));

  transitiondensities = std::make_shared<std::vector<Eigen::MatrixXd>>(0);

  // transitiondensities[0] : right transition density [Eta(R) + Xi(M)].
  transitiondensities->push_back(Eigen::MatrixXd::Zero(momospace, settings.nEigen));
  // transitiondensities[1] : left transition density [Xi(L)].
  transitiondensities->push_back(Eigen::MatrixXd::Zero(momospace, settings.nEigen));
  // transitiondensities[2] : state (col 0 = ground state) density matrices [Xi(lambda,N)].
  transitiondensities->push_back(Eigen::MatrixXd::Zero(momospace, settings.ccexdens ? (1 + settings.nEigen) : 1));
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculateStateMultipliers(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                            const LRSCFTaskSettings& settings,
                                                            std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                            Eigen::VectorXd& eigenvalues,
                                                            NonlinearSigmaCalculator sigmaCalculator,
                                                            const Eigen::VectorXd& diagonal) {
  printSubSectionTitle("State Multiplier");
  long nvnospace = (*eigenvectors)[2].rows();

  long voStart = 0;
  Eigen::MatrixXd rhs(nvnospace, settings.ccexdens ? (1 + settings.nEigen) : 1);
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
    };

    std::vector<Eigen::MatrixXd> eigenvectorsI(2);
    eigenvectorsI[0] = (*eigenvectors)[0].middleRows(voStart, nvno);
    eigenvectorsI[1] = (*eigenvectors)[1].middleRows(voStart, nvno);

    auto rhsI = ilrscf->getCC2Controller()->calculateExcitedStateLagrangeMultiplier(eigenvectorsI, eigenvalues);
    rhs.middleRows(voStart, rhsI.rows()) = rhsI;

    voStart += nvno;
  }
  NonlinearResponseSolver responsesolver(diagonal, settings.conv, settings.maxCycles, settings.maxSubspaceDimension,
                                         {0}, {rhs}, sigmaCalculator);
  (*eigenvectors)[2] = responsesolver.getEigenvectors()[0];

  // Distribute ground-state Lagrange multiplier.
  voStart = 0;
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
    };
    ilrscf->getCC2Controller()->setGroundStateLagrangeMultiplier((*eigenvectors)[2].col(0).middleRows(voStart, nvno));
    voStart += nvno;
  }
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculateStateDensities(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                          std::shared_ptr<std::vector<Eigen::MatrixXd>>& transitiondensities,
                                                          std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                          Eigen::VectorXd& eigenvalues) {
  printSubSectionTitle("State Densities");
  long voStart = 0, moStart = 0;
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0, momo = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
      momo += (nv_spin + no_spin) * (nv_spin + no_spin);
    };

    std::vector<Eigen::MatrixXd> eigenvectorsI(3);
    eigenvectorsI[0] = (*eigenvectors)[0].middleRows(voStart, nvno);
    eigenvectorsI[1] = (*eigenvectors)[1].middleRows(voStart, nvno);
    eigenvectorsI[2] = (*eigenvectors)[2].middleRows(voStart, nvno);

    std::vector<Eigen::MatrixXd> transitiondensitiesI(3);
    ilrscf->getCC2Controller()->calculateExcitedStateDensities(eigenvectorsI, eigenvalues, transitiondensitiesI);
    (*transitiondensities)[2].middleRows(moStart, momo) = transitiondensitiesI[2];

    std::string fileName = ilrscf->getSys()->getSystemPath() + ilrscf->getSys()->getSystemName() + "_cc2_dens.";
    if (lrscf.size() > 1)
      fileName += "fdec.";
    fileName += (SCFMode == RESTRICTED) ? "res." : "unres.";
    fileName += "h5";

    // H5F_ACC_TRUNC means create a new file or overwrite an existing file
    HDF5::H5File file(fileName, H5F_ACC_TRUNC);
    HDF5::save_scalar_attribute(file, "ID", ilrscf->getSys()->getSystemIdentifier());
    HDF5::save(file, "State Densities", transitiondensitiesI[2]);
    file.close();

    moStart += momo;
    voStart += nvno;
  }
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculateTransitionMultipliers(
    const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf, const LRSCFTaskSettings& settings,
    std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors, Eigen::VectorXd& eigenvalues,
    NonlinearSigmaCalculator sigmaCalculator, const Eigen::VectorXd& diagonal) {
  printSubSectionTitle("Transition Multiplier");
  long nvnospace = (*eigenvectors)[2].rows();

  long voStart = 0;
  std::vector<Eigen::MatrixXd> rhs(settings.nEigen, Eigen::MatrixXd(nvnospace, 1));
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
    };

    std::vector<Eigen::MatrixXd> eigenvectorsI(1);
    eigenvectorsI[0] = (*eigenvectors)[0].middleRows(voStart, nvno);

    auto rhsI = ilrscf->getCC2Controller()->calculateTransitionMomentLagrangeMultiplier(eigenvectorsI, eigenvalues);
    for (unsigned iEigen = 0; iEigen < settings.nEigen; ++iEigen) {
      rhs[iEigen].middleRows(voStart, nvno) = rhsI[iEigen];
    }

    voStart += nvno;
  }

  std::vector<double> frequencies(settings.nEigen);
  for (unsigned iEigen = 0; iEigen < settings.nEigen; ++iEigen) {
    frequencies[iEigen] = -eigenvalues(iEigen);
  }

  NonlinearResponseSolver responsesolver(diagonal, settings.conv, settings.maxCycles, settings.maxSubspaceDimension,
                                         frequencies, rhs, sigmaCalculator);
  auto multiplier = responsesolver.getEigenvectors();
  for (unsigned iEigen = 0; iEigen < settings.nEigen; ++iEigen) {
    (*eigenvectors)[3].col(iEigen) = multiplier[iEigen].col(0);
  }
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculateTransitionDensities(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                               std::shared_ptr<std::vector<Eigen::MatrixXd>>& transitiondensities,
                                                               std::shared_ptr<std::vector<Eigen::MatrixXd>>& eigenvectors,
                                                               Eigen::VectorXd& eigenvalues) {
  printSubSectionTitle("Transition Densities");
  long voStart = 0;
  long moStart = 0;
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0, momo = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
      momo += (nv_spin + no_spin) * (nv_spin + no_spin);
    };

    std::vector<Eigen::MatrixXd> eigenvectorsI(4);
    eigenvectorsI[0] = (*eigenvectors)[0].middleRows(voStart, nvno);
    eigenvectorsI[1] = (*eigenvectors)[1].middleRows(voStart, nvno);
    eigenvectorsI[2] = (*eigenvectors)[2].middleRows(voStart, nvno);
    eigenvectorsI[3] = (*eigenvectors)[3].middleRows(voStart, nvno);

    std::vector<Eigen::MatrixXd> transitiondensitiesI(2);
    ilrscf->getCC2Controller()->calculateTransitionDensities(eigenvectorsI, eigenvalues, transitiondensitiesI);

    (*transitiondensities)[0].middleRows(moStart, momo) = transitiondensitiesI[0];
    (*transitiondensities)[1].middleRows(moStart, momo) = transitiondensitiesI[1];

    std::string fileName = ilrscf->getSys()->getSystemPath() + ilrscf->getSys()->getSystemName() + "_cc2_dens.";
    if (lrscf.size() > 1)
      fileName += "fdec.";
    fileName += (SCFMode == RESTRICTED) ? "res." : "unres.";
    fileName += "h5";

    // H5F_ACC_RDWR means opening an existing file with read-write access
    HDF5::H5File file(fileName, H5F_ACC_RDWR);
    HDF5::save(file, "Right Transition Densities", transitiondensitiesI[0]);
    HDF5::save(file, "Left Transition Densities", transitiondensitiesI[1]);
    file.close();

    moStart += momo;
    voStart += nvno;
  }
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculatePerturbedAmplitudes(
    const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf, const LRSCFTaskSettings& settings,
    std::shared_ptr<std::vector<Eigen::MatrixXd>>& solutionvectors, const std::vector<double>& frequencies,
    NonlinearSigmaCalculator sigmaCalculator, const Eigen::VectorXd& diagonal, Eigen::MatrixXd& dipoles) {
  printSubSectionTitle("Perturbation Vectors");

  long voStart = 0, moStart = 0;
  std::vector<Eigen::MatrixXd> pert(frequencies.size(), Eigen::MatrixXd(0, 6));
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0, momo = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
      momo += (nv_spin + no_spin) * (nv_spin + no_spin);
    };

    auto pertI = ilrscf->getCC2Controller()->getPerturbation(frequencies, dipoles.middleRows(moStart, momo));

    for (unsigned iFreq = 0; iFreq < frequencies.size(); ++iFreq) {
      pert[iFreq].conservativeResize(pert[iFreq].rows() + pertI[iFreq].rows(), Eigen::NoChange);
      pert[iFreq].middleRows(voStart, pertI[iFreq].rows()) = pertI[iFreq];
    }
    moStart += momo;
    voStart += nvno;
  }

  printSubSectionTitle("Perturbed Amplitudes");
  NonlinearResponseSolver responsesolver(diagonal, settings.conv, settings.maxCycles, settings.maxSubspaceDimension,
                                         frequencies, pert, sigmaCalculator);
  solutionvectors = std::make_shared<std::vector<Eigen::MatrixXd>>(responsesolver.getEigenvectors());
}

template<Options::SCF_MODES SCFMode>
void CC2HelperFunctions<SCFMode>::calculatePerturbedDensities(
    const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf, const LRSCFTaskSettings& settings,
    std::shared_ptr<std::vector<Eigen::MatrixXd>>& solutionvectors,
    std::shared_ptr<std::vector<Eigen::MatrixXd>>& perturbeddensities, const std::vector<double>& frequencies,
    Eigen::MatrixXd& dipoles, std::vector<Eigen::Matrix3d>& Fdipdip, std::vector<Eigen::Matrix3d>& Fdipmag) {
  printSubSectionTitle("Perturbed Densities and F-Contractions");
  const unsigned nFreqs = settings.frequencies.size();
  long voStart = 0;
  long moStart = 0;
  std::vector<Eigen::MatrixXd> densities(frequencies.size(), Eigen::MatrixXd(0, 6));
  Fdipdip = std::vector<Eigen::Matrix3d>(nFreqs, Eigen::Matrix3d::Zero());
  Fdipmag = std::vector<Eigen::Matrix3d>(nFreqs, Eigen::Matrix3d::Zero());
  for (auto& ilrscf : lrscf) {
    auto no = ilrscf->getNOccupied();
    auto nv = ilrscf->getNVirtual();
    long nvno = 0, momo = 0;
    for_spin(nv, no) {
      nvno += nv_spin * no_spin;
      momo += (nv_spin + no_spin) * (nv_spin + no_spin);
    };

    std::vector<Eigen::MatrixXd> solutionI;
    std::vector<Eigen::MatrixXd> dipI;
    std::vector<Eigen::MatrixXd> magI;
    for (unsigned iFreq = 0; iFreq < frequencies.size(); ++iFreq) {
      solutionI.push_back((*solutionvectors)[iFreq].middleRows(voStart, nvno));
      dipI.push_back(solutionI[iFreq].leftCols(3));
      magI.push_back(solutionI[iFreq].rightCols(3));
    }

    // Get perturbed densities.
    std::vector<Eigen::MatrixXd> densI;
    ilrscf->getCC2Controller()->calcPerturbedDensities(solutionI, frequencies, densI, dipoles.middleRows(moStart, momo));

    for (unsigned iFreq = 0; iFreq < frequencies.size(); ++iFreq) {
      densities[iFreq].conservativeResize(densities[iFreq].rows() + densI[iFreq].rows(), Eigen::NoChange);
      densities[iFreq].middleRows(moStart, densI[iFreq].rows()) = densI[iFreq];
    }

    Eigen::Ref<Eigen::MatrixXd> electricsI = dipoles.middleRows(moStart, momo).leftCols(3);
    Eigen::Ref<Eigen::MatrixXd> magneticsI = dipoles.middleRows(moStart, momo).rightCols(3);

    printf("  F-matrix contractions.\n");
    auto F = ilrscf->getCC2Controller()->getFContractions(dipI, magI, settings.frequencies, electricsI, magneticsI,
                                                          settings.gauge);
    auto FdipI = F.first;
    auto FmagI = F.second;

    for (unsigned iFreq = 0; iFreq < nFreqs; ++iFreq) {
      Fdipdip[iFreq] += FdipI[iFreq];
      Fdipmag[iFreq] += FmagI[iFreq];
    }

    moStart += momo;
    voStart += nvno;
  }

  perturbeddensities = std::make_shared<std::vector<Eigen::MatrixXd>>(4, Eigen::MatrixXd::Zero(moStart, nFreqs * 3));
  for (unsigned iFreq = 0; iFreq < nFreqs; ++iFreq) {
    unsigned plusIndex = iFreq;
    unsigned minusIndex = (frequencies[iFreq] != 0) ? iFreq + nFreqs : iFreq;
    // [0] : electric dipole at plus frequency.
    (*perturbeddensities)[0].middleCols(iFreq * 3, 3) = densities[plusIndex].leftCols(3);
    // [1] : electric dipole at minus frequency.
    (*perturbeddensities)[1].middleCols(iFreq * 3, 3) = densities[minusIndex].leftCols(3);
    // [2] : magnetic dipole at plus frequency.
    (*perturbeddensities)[2].middleCols(iFreq * 3, 3) = densities[plusIndex].rightCols(3);
    // [3] : magnetic dipole at minus frequency.
    (*perturbeddensities)[3].middleCols(iFreq * 3, 3) = densities[minusIndex].rightCols(3);
  }
}

template class CC2HelperFunctions<Options::SCF_MODES::RESTRICTED>;
template class CC2HelperFunctions<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
