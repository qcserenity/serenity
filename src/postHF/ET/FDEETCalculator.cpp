/**
 * @file   FDEETCalculator.cpp
 *
 * @date   Apr. 20, 2020
 * @author Patrick Eschenbach
 * @copyright \n
 *  This file is part of the program Serenity.\n\n
 *  Serenity is free software: you can redistribute it and/or modify
 *  it under the terms of the LGNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.\n\n
 *  Serenity is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.\n\n
 *  You should have received a copy of the LGNU Lesser General
 *  Public License along with Serenity.
 *  If not, see <http://www.gnu.org/licenses/>.\n

/* Include Class Header*/
#include "postHF/ET/FDEETCalculator.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "data/grid/ElectrostaticPotentialOnGridController.h"
#include "dft/dispersionCorrection/DispersionCorrectionCalculator.h"
#include "energies/EnergyContributions.h"
#include "geometry/MolecularSurfaceController.h"
#include "io/CubeFileWriter.h"
#include "io/Filesystem.h"
#include "io/FormattedOutputStream.h"
#include "io/HDF5.h"
#include "parameters/Constants.h"
#include "potentials/ERIPotential.h"
#include "potentials/ExchangePotential.h"
#include "potentials/FuncPotential.h"
#include "potentials/HCorePotential.h"
#include "potentials/ZeroPotential.h"
#include "potentials/bundles/DFTPotentials.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/PlotTask.h"
/* Include Std and External Headers */
#include <algorithm>

namespace Serenity {

FDEETCalculator::FDEETCalculator() : _dMatFiles(nullptr) {
}

std::vector<Eigen::MatrixXd> FDEETCalculator::getMoOverlap() {
  return _overlapMO;
}

void FDEETCalculator::setMoOverlap(std::vector<Eigen::MatrixXd> overlap) {
  _overlapMO = overlap;
}

void FDEETCalculator::calculateMoOverlap(const unsigned nStates, const unsigned nSysPerState,
                                         const std::vector<unsigned> jointIndex,
                                         const std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>>& stateCoeffs,
                                         const std::vector<std::vector<std::shared_ptr<SystemController>>>& stateVector,
                                         const std::shared_ptr<std::vector<unsigned>>& indexVector,
                                         const Eigen::MatrixXi nOccState) {
  auto indexVec = *indexVector;
  // initialize and fill overlap super matrix
  std::vector<Eigen::MatrixXd> overlapMO(2);
  overlapMO[0] = Eigen::MatrixXd::Zero(nOccState.col(0).sum(), nOccState.col(0).sum());
  overlapMO[1] = Eigen::MatrixXd::Zero(nOccState.col(1).sum(), nOccState.col(1).sum());
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, int> dummy;
  auto& libint = Libint::getInstance();
  std::vector<unsigned> iCount = {0, 0};
  std::vector<unsigned> jCount = {0, 0};
  std::vector<unsigned> iStart = {0, 0};
  std::vector<unsigned> jStart = {0, 0};
  std::vector<unsigned> iBlockI = {0, 0};
  std::vector<unsigned> jBlockI = {0, 0};
  std::vector<unsigned> iBlockJ = {0, 0};
  std::vector<unsigned> jBlockJ = {0, 0};
  unsigned iSpin = 0;
  for (unsigned iState = 0; iState < nStates; ++iState) {
    for (unsigned jState = 0; jState < nStates; ++jState) {
      for (unsigned iSys = 0; iSys < nSysPerState; ++iSys) {
        auto sysI = stateVector[indexVec[iState]][iSys];
        std::vector<Eigen::MatrixXd> coeffI = {*(stateCoeffs[iState][0]), *(stateCoeffs[iState][1])};
        auto nOccI = sysI->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
        unsigned nBFI = sysI->getBasisController()->getNBasisFunctions();
        for (unsigned jSys = 0; jSys < nSysPerState; ++jSys) {
          auto sysJ = stateVector[indexVec[jState]][jSys];
          std::vector<Eigen::MatrixXd> coeffJ = {*(stateCoeffs[jState][0]), *(stateCoeffs[jState][1])};
          auto nOccJ = sysJ->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
          unsigned nBFJ = sysJ->getBasisController()->getNBasisFunctions();
          auto overlapAO =
              libint.compute1eInts(LIBINT_OPERATOR::overlap, sysJ->getBasisController(), sysI->getBasisController());
          iSpin = 0;
          for_spin(nOccI, nOccJ) {
            /* joint calculation */
            if (jointIndex.size() == 0) {
              overlapMO[iSpin].block(iCount[iSpin], jCount[iSpin], nOccI_spin, nOccJ_spin) =
                  coeffI[iSpin].block(iBlockI[iSpin], jBlockI[iSpin], nBFI, nOccI_spin).transpose() * overlapAO *
                  coeffJ[iSpin].block(iBlockJ[iSpin], jBlockJ[iSpin], nBFJ, nOccJ_spin);
            }
            /* disjoint calculation */
            else {
              bool containsISys = std::find(jointIndex.begin(), jointIndex.end(), iSys + 1) != jointIndex.end();
              bool containsJSys = std::find(jointIndex.begin(), jointIndex.end(), jSys + 1) != jointIndex.end();
              bool diag = (iSys == jSys) ? 1 : 0;
              bool offDiag = (iSys != jSys and containsISys and containsJSys) ? 1 : 0;
              if (diag or offDiag) {
                overlapMO[iSpin].block(iCount[iSpin], jCount[iSpin], nOccI_spin, nOccJ_spin) =
                    coeffI[iSpin].block(iBlockI[iSpin], jBlockI[iSpin], nBFI, nOccI_spin).transpose() * overlapAO *
                    coeffJ[iSpin].block(iBlockJ[iSpin], jBlockJ[iSpin], nBFJ, nOccJ_spin);
              }
              else {
                overlapMO[iSpin].block(iCount[iSpin], jCount[iSpin], nOccI_spin, nOccJ_spin) =
                    Eigen::MatrixXd::Zero(nOccI_spin, nOccJ_spin);
              }
            }
            jCount[iSpin] += nOccJ_spin;
            iBlockJ[iSpin] += nBFJ;
            jBlockJ[iSpin] += nOccJ_spin;
            ++iSpin;
          };
        } /* jSys */
        iSpin = 0;
        for_spin(nOccI) {
          iCount[iSpin] += nOccI_spin;
          jCount[iSpin] = jStart[iSpin];
          iBlockI[iSpin] += nBFI;
          jBlockI[iSpin] += nOccI_spin;
          iBlockJ[iSpin] = 0;
          jBlockJ[iSpin] = 0;
          ++iSpin;
        };
      } /* iSys */
      iSpin = 0;
      for_spin(dummy) {
        (void)dummy_spin;
        jStart[iSpin] += nOccState(jState, iSpin);
        iCount[iSpin] = iStart[iSpin];
        jCount[iSpin] = jStart[iSpin];
        iBlockI[iSpin] = 0;
        jBlockI[iSpin] = 0;
        iBlockJ[iSpin] = 0;
        jBlockJ[iSpin] = 0;
        ++iSpin;
      };
    } /* jState */
    iSpin = 0;
    for_spin(dummy) {
      (void)dummy_spin;
      jStart[iSpin] = 0;
      iStart[iSpin] += nOccState(iState, iSpin);
      iCount[iSpin] = iStart[iSpin];
      jCount[iSpin] = jStart[iSpin];
      ++iSpin;
    };
  } /* iState */
  this->setMoOverlap(overlapMO);
  OutputControl::vOut << "Dim. overlapMO alpha: " << overlapMO[0].rows() << "x" << overlapMO[0].cols() << std::endl;
  OutputControl::vOut << "Dim. overlapMO beta:  " << overlapMO[1].rows() << "x" << overlapMO[1].cols() << std::endl;
  OutputControl::vOut << std::endl;
}

std::vector<Eigen::MatrixXd> FDEETCalculator::getPseudoInverse() {
  return _pseudoInverse;
}

void FDEETCalculator::setPseudoInverse(std::vector<Eigen::MatrixXd> pseudo) {
  _pseudoInverse = pseudo;
}

void FDEETCalculator::calcBlockedMoorePenrose(const std::vector<Eigen::MatrixXd>& matrix, Eigen::MatrixXi nOccState,
                                              unsigned nStatesCouple, double integralThreshold) {
  // Calculate pseudo inverse (has same dimensions as overlapMO) and determinants of transition-overlap matrix
  std::vector<Eigen::MatrixXd> pseudoInverse(2);
  pseudoInverse[0] = Eigen::MatrixXd::Zero(nOccState.col(0).sum(), nOccState.col(0).sum());
  pseudoInverse[1] = Eigen::MatrixXd::Zero(nOccState.col(1).sum(), nOccState.col(1).sum());
  Eigen::MatrixXd determinants = Eigen::MatrixXd::Ones(nStatesCouple, nStatesCouple);
  std::vector<unsigned> iCount = {0, 0};
  std::vector<unsigned> jCount = {0, 0};
  unsigned valsBelowThresh = 0;
  unsigned vals = 0;
  // Calculate pseudo-inverse and derterminant of overlapMO blocks, sorted by spin
  for (unsigned iSpin = 0; iSpin < 2; ++iSpin) {
    for (unsigned iState = 0; iState < nStatesCouple; ++iState) {
      for (unsigned jState = 0; jState < nStatesCouple; ++jState) {
        // Block-wise Jacobi SVD for each state combination in overlapMO;
        Eigen::MatrixXd block =
            matrix[iSpin].block(iCount[iSpin], jCount[iSpin], nOccState(iState, iSpin), nOccState(jState, iSpin));
        Eigen::MatrixXd singMatrix = Eigen::MatrixXd::Zero(block.rows(), block.cols());
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(block, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::VectorXd sing = svd.singularValues();
        determinants(iState, jState) *= sing.prod() * fabs(svd.matrixU().determinant()) * fabs(svd.matrixV().determinant());
        // Invert singular values that are above a certain threshold and add to diagonal matrix
        double invThreshold = integralThreshold * sing.size() * sing.maxCoeff();
        OutputControl::vOut << "  Inversion threshold: " << invThreshold << std::endl;
        for (unsigned i = 0; i < sing.size(); ++i) {
          vals += 1;
          if (fabs(sing(i)) < invThreshold)
            valsBelowThresh += 1;
          singMatrix(i, i) = fabs(sing(i)) >= invThreshold ? 1 / sing(i) : 0.0;
        }
        // calculate blocks of pseudo inverse
        pseudoInverse[iSpin].block(iCount[iSpin], jCount[iSpin], nOccState(iState, iSpin), nOccState(jState, iSpin)) =
            svd.matrixU() * singMatrix * svd.matrixV().transpose();
        jCount[iSpin] += nOccState(jState, iSpin);
      } /* jState */
      iCount[iSpin] += nOccState(iState, iSpin);
      jCount[iSpin] = 0;
    } /* iState */
  }   /* iSpin */
  this->setPseudoInverse(pseudoInverse);
  this->setDeterminants(determinants);
  OutputControl::vOut << "Dim. pseudoInverse alpha: " << pseudoInverse[0].rows() << "x" << pseudoInverse[0].cols()
                      << std::endl;
  OutputControl::vOut << "Dim. pseudoInverse beta:  " << pseudoInverse[1].rows() << "x" << pseudoInverse[1].cols()
                      << std::endl;
  OutputControl::n.printf("  %5i of %5i Singular values were below the inversion threshold.\n\n", valsBelowThresh, vals);
}

Eigen::MatrixXd FDEETCalculator::getDeterminants() {
  return _determinants;
}

void FDEETCalculator::setDeterminants(Eigen::MatrixXd dets) {
  _determinants = dets;
}

std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> FDEETCalculator::getTransDensMats() {
  return _densMats;
}

void FDEETCalculator::setTransDensMats(std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> mat) {
  _densMats = mat;
}

void FDEETCalculator::calcTransDensMats(const std::vector<Eigen::MatrixXd>& pseudoInverse,
                                        std::shared_ptr<SystemController> superSystem,
                                        const std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>>& coeffs,
                                        Eigen::MatrixXi nOccState, unsigned nStatesCouple) {
  unsigned iSpin = 0;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, int> dummy;
  // Build density matrices of only zeros
  std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>> densMats(
      nStatesCouple, std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>(
                         nStatesCouple, DensityMatrix<Options::SCF_MODES::UNRESTRICTED>(superSystem->getBasisController())));
  std::vector<unsigned> iCount = {0, 0};
  std::vector<unsigned> jCount = {0, 0};
  // get transition density matrix P for the state combinations (ij)
  // P(ij) = C(i) * S(ij) * C(j)^T                            Eq.(1)
  // for each state combination and spin
  for (unsigned iState = 0; iState < nStatesCouple; ++iState) {
    for (unsigned jState = 0; jState < nStatesCouple; ++jState) {
      auto& densMat = densMats[iState][jState];
      iSpin = 0;
      // loop over alpha and beta
      for_spin(densMat) {
        // Eq. (1)
        densMat_spin =
            *coeffs[iState][iSpin] *
            pseudoInverse[iSpin].block(iCount[iSpin], jCount[iSpin], nOccState(iState, iSpin), nOccState(jState, iSpin)) *
            coeffs[jState][iSpin]->transpose();
        ++iSpin;
      };
      // reset loop and block vars
      iSpin = 0;
      for_spin(dummy) {
        (void)dummy_spin;
        jCount[iSpin] += nOccState(jState, iSpin);
        ++iSpin;
      };
    }
    iSpin = 0;
    for_spin(dummy) {
      (void)dummy_spin;
      iCount[iSpin] += nOccState(iState, iSpin);
      jCount[iSpin] = 0;
      ++iSpin;
    };
  }
  this->setTransDensMats(densMats);
}

void FDEETCalculator::calcTransDensMatsDisk(const std::vector<Eigen::MatrixXd>& pseudoInverse,
                                            std::shared_ptr<SystemController> superSystem,
                                            const std::vector<std::vector<std::shared_ptr<Eigen::MatrixXd>>>& coeffs,
                                            Eigen::MatrixXi nOccState, unsigned nStatesCouple) {
  std::vector<std::string> dmatFiles;
  unsigned iSpin = 0;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, int> dummy;
  // Build density matrices of only zeros and print to file
  auto transDmatController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(
      DensityMatrix<Options::SCF_MODES::UNRESTRICTED>(superSystem->getBasisController()));
  for (unsigned iState = 0; iState < nStatesCouple; ++iState) {
    for (unsigned jState = 0; jState < nStatesCouple; ++jState) {
      transDmatController->toHDF5(superSystem->getSystemPath() + "t_dmat" + iState + jState,
                                  superSystem->getSystemIdentifier());
      std::string path = superSystem->getSystemPath() + "t_dmat" + iState + jState + ".dmat.unres.h5";
      dmatFiles.push_back(path.c_str());
    }
  }
  _dMatFiles = std::make_shared<std::vector<std::string>>(dmatFiles);
  std::vector<unsigned> iCount = {0, 0};
  std::vector<unsigned> jCount = {0, 0};
  unsigned ij = 0;
  for (unsigned iState = 0; iState < nStatesCouple; ++iState) {
    for (unsigned jState = 0; jState < nStatesCouple; ++jState) {
      DensityMatrix<Options::SCF_MODES::UNRESTRICTED> densMat(superSystem->getBasisController());
      std::string pathUnres = (*_dMatFiles)[ij];
      ++ij;
      std::vector<Eigen::MatrixXd> temp(2, Eigen::MatrixXd::Zero(superSystem->getBasisController()->getNBasisFunctions(),
                                                                 superSystem->getBasisController()->getNBasisFunctions()));
      HDF5::Filepath name(pathUnres);
      HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
      HDF5::dataset_exists(file, "densityMatrix_alpha");
      HDF5::load(file, "densityMatrix_alpha", temp[0]);
      HDF5::dataset_exists(file, "densityMatrix_beta");
      HDF5::load(file, "densityMatrix_beta", temp[1]);
      file.close();
      unsigned iSpin = 0;
      for_spin(densMat) {
        densMat_spin = temp[iSpin];
        ++iSpin;
      };
      iSpin = 0;
      // loop over alpha and beta
      for_spin(densMat) {
        // P(ij) = C(i) * S(ij) * C(j)^T                            Eq.(1)
        auto& temp =
            *coeffs[iState][iSpin] *
            pseudoInverse[iSpin].block(iCount[iSpin], jCount[iSpin], nOccState(iState, iSpin), nOccState(jState, iSpin)) *
            coeffs[jState][iSpin]->transpose();
        // remove numerical noise from matrix
        densMat_spin = (temp + temp.transpose()) * 0.5;
        ++iSpin;
      };
      auto transDmatController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(densMat);
      transDmatController->toHDF5(superSystem->getSystemPath() + "t_dmat" + iState + jState,
                                  superSystem->getSystemIdentifier());
      // reset loop and block vars
      iSpin = 0;
      for_spin(dummy) {
        (void)dummy_spin;
        jCount[iSpin] += nOccState(jState, iSpin);
        ++iSpin;
      };
    }
    iSpin = 0;
    for_spin(dummy) {
      (void)dummy_spin;
      iCount[iSpin] += nOccState(iState, iSpin);
      jCount[iSpin] = 0;
      ++iSpin;
    };
  }
}

std::shared_ptr<std::vector<std::string>> FDEETCalculator::getDensityMatrixFiles() {
  return _dMatFiles;
}

std::shared_ptr<EnergyComponentController>
FDEETCalculator::calcDFTEnergy(std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatC,
                               std::shared_ptr<SystemController> superSystem, bool isDiagonal, bool useHFCoupling) {
  auto D = dMatC->getDensityMatrix();
  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> Dtranspose(superSystem->getBasisController());
  DensityMatrix<Options::SCF_MODES::UNRESTRICTED> Dsymm(superSystem->getBasisController());
  for_spin(Dsymm, Dtranspose, D) {
    Dtranspose_spin = D_spin.transpose();
    Eigen::MatrixXd tmp = 1.0 * D_spin;
    Dsymm_spin = 0.5 * (tmp + tmp.transpose());
  };
  dMatC->setDensityMatrix(Dsymm);
  // Hcore potential
  auto hcore = std::make_shared<HCorePotential<Options::SCF_MODES::UNRESTRICTED>>(superSystem);
  // XC potential
  auto functional = superSystem->getSettings().customFunc.basicFunctionals.size()
                        ? Functional(superSystem->getSettings().customFunc)
                        : resolveFunctional(superSystem->getSettings().dft.functional);
  auto Vxc = std::make_shared<FuncPotential<UNRESTRICTED>>(superSystem->getSharedPtr(), dMatC,
                                                           superSystem->getGridController(), functional);

  // J potential
  auto J = std::make_shared<ERIPotential<Options::SCF_MODES::UNRESTRICTED>>(
      superSystem, dMatC, (isDiagonal) ? functional.getHfExchangeRatio() : 0.0,
      superSystem->getSettings().basis.integralThreshold, superSystem->getSettings().basis.integralIncrementThresholdStart,
      superSystem->getSettings().basis.integralIncrementThresholdEnd, superSystem->getSettings().basis.incrementalSteps,
      true, functional.getLRExchangeRatio(), functional.getRangeSeparationParameter());
  // PCM potential
  std::shared_ptr<Potential<Options::SCF_MODES::UNRESTRICTED>> pcm(
      new ZeroPotential<Options::SCF_MODES::UNRESTRICTED>(superSystem->getBasisController()));
  // DFT potential
  auto dftpot = std::make_shared<DFTPotentials<Options::SCF_MODES::UNRESTRICTED>>(
      hcore, J, Vxc, pcm, superSystem->getGeometry(), dMatC, superSystem->getSettings().basis.integralThreshold);
  auto ec = std::make_shared<EnergyComponentController>();
  dftpot->getFockMatrix(dMatC->getDensityMatrix(), ec);
  // Exchange for off-diagonal elements
  if ((functional.isHybrid() || useHFCoupling) && !isDiagonal) {
    dMatC->setDensityMatrix(D);
    double hfRatio = (useHFCoupling) ? 1.0 : functional.getHfExchangeRatio();
    auto K = std::make_shared<ExchangePotential<Options::SCF_MODES::UNRESTRICTED>>(
        superSystem, dMatC, hfRatio, superSystem->getSettings().basis.integralThreshold,
        superSystem->getSettings().basis.integralIncrementThresholdStart,
        superSystem->getSettings().basis.integralIncrementThresholdEnd,
        superSystem->getSettings().basis.incrementalSteps, true, true);
    ec->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXACT_EXCHANGE, K->getEnergy(Dtranspose));
    if (useHFCoupling)
      ec->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXCHANGE_CORRELATION_NO_HF, 0.0);
  }
  // Dispersion
  auto dispCorr = DispersionCorrectionCalculator::calcDispersionEnergyCorrection(
      superSystem->getSettings().dft.dispersion, superSystem->getGeometry(), superSystem->getSettings().dft.functional);
  ec->addOrReplaceComponent(ENERGY_CONTRIBUTIONS::KS_DFT_DISPERSION_CORRECTION, dispCorr);
  // Energy contributions
  if (GLOBAL_PRINT_LEVEL == Options::GLOBAL_PRINT_LEVELS::VERBOSE)
    ec->printAllComponents();
  return ec;
}

Eigen::MatrixXd FDEETCalculator::getHamiltonian() {
  return _hamiltonian;
}

void FDEETCalculator::setHamiltonian(Eigen::MatrixXd newH) {
  _hamiltonian = newH;
}

void FDEETCalculator::calcHamiltonian(const std::vector<std::vector<DensityMatrix<Options::SCF_MODES::UNRESTRICTED>>>& transDensMats,
                                      std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple,
                                      bool useHFCoupling) {
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nStatesCouple, nStatesCouple);
  for (unsigned iState = 0; iState < nStatesCouple; ++iState) {
    for (unsigned jState = 0; jState < nStatesCouple; ++jState) {
      auto& densMat = transDensMats[iState][jState];
      unsigned iSpin = 0;
      for_spin(densMat) {
        OutputControl::vOut << "Dim. P(" << iState + 1 << jState + 1 << ") and spin " << iSpin << ": "
                            << densMat_spin.rows() << "x" << densMat_spin.cols() << std::endl;
        iSpin += 1;
      };
      OutputControl::vOut << std::endl;
      auto dMatController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(densMat);
      OutputControl::vOut << "Energy Components of Diabatic State Combination: " << iState + 1 << " and " << jState + 1
                          << std::endl;
      OutputControl::vOut << std::string(100, '-') << std::endl;
      auto ec = this->calcDFTEnergy(dMatController, superSystem, iState == jState, useHFCoupling);
      if (useHFCoupling && iState == jState) {
        double correction = this->calcCorrection(superSystem, ec, dMatController);
        H.row(iState) += correction * Eigen::VectorXd::Ones(nStatesCouple);
        H.col(jState) += correction * Eigen::VectorXd::Ones(nStatesCouple);
        H(iState, jState) = 0.0;
      }
      H(iState, jState) += ec->getTotalEnergy();
      OutputControl::vOut << std::endl;
    }
  }
  OutputControl::nOut << "...Success" << std::endl;
  this->setHamiltonian(H);
}

void FDEETCalculator::calcHamiltonianDisk(std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple,
                                          bool useHFCoupling) {
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nStatesCouple, nStatesCouple);
  for (unsigned iState = 0; iState < nStatesCouple; ++iState) {
    for (unsigned jState = 0; jState < nStatesCouple; ++jState) {
      DensityMatrix<Options::SCF_MODES::UNRESTRICTED> densMat(superSystem->getBasisController());
      std::string pathUnres = superSystem->getSystemPath() + "t_dmat" + iState + jState + ".dmat.unres.h5";
      std::vector<Eigen::MatrixXd> temp(2, Eigen::MatrixXd::Zero(superSystem->getBasisController()->getNBasisFunctions(),
                                                                 superSystem->getBasisController()->getNBasisFunctions()));
      HDF5::Filepath name(pathUnres);
      HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY);
      HDF5::dataset_exists(file, "densityMatrix_alpha");
      HDF5::load(file, "densityMatrix_alpha", temp[0]);
      HDF5::dataset_exists(file, "densityMatrix_beta");
      HDF5::load(file, "densityMatrix_beta", temp[1]);
      file.close();
      unsigned iSpin = 0;
      for_spin(densMat) {
        densMat_spin = temp[iSpin];
        ++iSpin;
      };
      auto dMatController = std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(densMat);
      OutputControl::vOut << "Energy Components of Diabatic State Combination: " << iState + 1 << " and " << jState + 1
                          << std::endl;
      OutputControl::vOut << std::string(100, '-') << std::endl;
      auto ec = this->calcDFTEnergy(dMatController, superSystem, iState == jState, useHFCoupling);
      if (useHFCoupling && iState == jState) {
        double correction = this->calcCorrection(superSystem, ec, dMatController);
        H.row(iState) += correction * Eigen::VectorXd::Ones(nStatesCouple);
        H.col(jState) += correction * Eigen::VectorXd::Ones(nStatesCouple);
        H(iState, jState) = 0.0;
      }
      H(iState, jState) += ec->getTotalEnergy();
      OutputControl::vOut << std::endl;
    }
  }
  OutputControl::nOut << "...Success" << std::endl;
  this->setHamiltonian(H);
}

double FDEETCalculator::calcCorrection(std::shared_ptr<SystemController> superSystem,
                                       std::shared_ptr<EnergyComponentController> ec,
                                       std::shared_ptr<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>> dMatController) {
  double correction = 0.0;
  auto functional = superSystem->getSettings().customFunc.basicFunctionals.size()
                        ? Functional(superSystem->getSettings().customFunc)
                        : resolveFunctional(superSystem->getSettings().dft.functional);
  if (functional.isHybrid()) {
    // XC contribution
    correction += 0.5 * ec->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXCHANGE_CORRELATION_NO_HF);
    correction += 0.5 * ec->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXACT_EXCHANGE);
    // Exact exchange contribution
    correction -=
        0.5 * (ec->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXACT_EXCHANGE) / functional.getHfExchangeRatio());
  }
  else {
    // XC contribution
    correction += 0.5 * ec->getEnergyComponent(ENERGY_CONTRIBUTIONS::KS_DFT_EXCHANGE_CORRELATION_NO_HF);
    // Exact exchange contribution
    auto K = std::make_shared<ExchangePotential<Options::SCF_MODES::UNRESTRICTED>>(
        superSystem, dMatController, 1.0, superSystem->getSettings().basis.integralThreshold,
        superSystem->getSettings().basis.integralIncrementThresholdStart,
        superSystem->getSettings().basis.integralIncrementThresholdEnd, superSystem->getSettings().basis.incrementalSteps);
    correction -= 0.5 * K->getEnergy(dMatController->getDensityMatrix());
  }
  return correction;
}

void FDEETCalculator::solveEigenValueProblem(Eigen::MatrixXd H, Eigen::MatrixXd determinants, unsigned nStatesCouple,
                                             bool analyticCoupling, bool isPrime) {
  Eigen::MatrixXd C;
  Eigen::VectorXd E;
  Eigen::setNbThreads(1);
  // Analytical electronic coupling
  if (nStatesCouple == 1) {
    C = Eigen::MatrixXd::Ones(1, 1);
    E = H(0, 0) * Eigen::VectorXd::Ones(1);
  }
  else {
    if (nStatesCouple == 2 && analyticCoupling) {
      auto SPrime = determinants;
      Eigen::MatrixXd HPrime = H;
      if (!isPrime) {
        SPrime(0, 0) = 1.0;
        SPrime(1, 1) = 1.0;
        SPrime(0, 1) = determinants(0, 1) / sqrt(determinants(0, 0) * determinants(1, 1));
        SPrime(1, 0) = determinants(1, 0) / sqrt(determinants(0, 0) * determinants(1, 1));
        HPrime = H.cwiseProduct(SPrime);
      }
      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver(HPrime, SPrime);
      C = solver.eigenvectors().real();
      E = solver.eigenvalues().real();
      double coupling =
          (1 / (1 - pow(SPrime(0, 1), 2))) * fabs(HPrime(0, 1) - SPrime(0, 1) * ((HPrime(0, 0) + HPrime(1, 1)) / 2));
      double lree = sqrt(pow(HPrime(0, 0) - HPrime(1, 1), 2) / (1 - pow(SPrime(0, 1), 2)) + 4 * pow(coupling, 2));
      this->setAnalyticalCoupling(coupling);
      this->setLongRangeExEnergy(lree);
    }
    else {
      Eigen::MatrixXd HPrime = H;
      if (!isPrime)
        HPrime = H.cwiseProduct(determinants);
      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver(HPrime, determinants);
      C = solver.eigenvectors().real();
      E = solver.eigenvalues().real();
    }
    // sort eigenvectors and eigenvalues in ascending order
    for (unsigned i = 0; i < E.size(); ++i) {
      unsigned iMin;
      E.tail(E.size() - i).minCoeff(&iMin);
      E.row(i).swap(E.row(iMin + i));
      C.col(i).swap(C.col(iMin + i));
    }
    E = (E.head(E.size())).eval();
    C = (C.leftCols(E.size())).eval();
    // normalize eigenvectors
    for (unsigned iCol = 0; iCol < C.cols(); ++iCol) {
      double norm = sqrt(C.col(iCol).transpose() * determinants * C.col(iCol));
      C.col(iCol) = C.col(iCol) / norm;
    }
  }
  this->setEigenValues(E);
  this->setLinCoeffs(C);
  Eigen::setNbThreads(0);
}

Eigen::VectorXd FDEETCalculator::getEigenValues() {
  return _eigenValues;
}

void FDEETCalculator::setEigenValues(Eigen::VectorXd newValues) {
  _eigenValues = newValues;
}

Eigen::MatrixXd FDEETCalculator::getLinCoeffs() {
  return _linCoeffs;
}

void FDEETCalculator::setLinCoeffs(Eigen::MatrixXd newCoeffs) {
  _linCoeffs = newCoeffs;
}

double FDEETCalculator::getAnalyticalCoupling() {
  return _analyticalCoupling;
}

void FDEETCalculator::setAnalyticalCoupling(double newC) {
  _analyticalCoupling = newC;
}

double FDEETCalculator::getLongRangeExEnergy() {
  return _longRangeExEnergy;
}

void FDEETCalculator::setLongRangeExEnergy(double newE) {
  _longRangeExEnergy = newE;
}

void FDEETCalculator::printResults(unsigned nStatesCouple, bool analyticCoupling) {
  auto overlapMO = this->getMoOverlap();
  auto pseudoInverse = this->getPseudoInverse();
  auto determinants = this->getDeterminants();
  auto H = this->getHamiltonian();
  auto C = this->getLinCoeffs();
  auto E = this->getEigenValues();
  auto coupling = this->getAnalyticalCoupling();
  auto lree = this->getLongRangeExEnergy();
  printSubSectionTitle("Results for Coupling Problem");
  printTableHead("Normalized Linear Combination Coefficients for Construction of Adiabatic States");
  OutputControl::nOut << std::string(10, ' ');
  for (unsigned col = 0; col < C.cols(); ++col) {
    OutputControl::nOut << "Psi_" << col << std::string(5, ' ');
  }
  OutputControl::nOut.flush();
  OutputControl::nOut << std::endl;
  for (unsigned row = 0; row < C.rows(); ++row) {
    OutputControl::nOut << "Phi_" << row << std::string(2, ' ');
    OutputControl::nOut.flush();
    for (unsigned col = 0; col < C.cols(); ++col) {
      OutputControl::n.printf("%10.5f", C(row, col));
    }
    OutputControl::nOut << std::endl;
  }
  OutputControl::nOut
      << std::endl
      << "Note: Psi_i are the adiabatic states, where Phi_i are the diabatic states used for the linear "
         "combination.\n"
      << std::endl;
  printTableHead("Hamiltonian/a.u.");
  OutputControl::nOut << std::fixed << std::setprecision(5) << H << std::endl;
  OutputControl::nOut << std::endl;
  printTableHead("Overlap of Diabatic States/a.u.");
  OutputControl::nOut << std::fixed << std::setprecision(5) << determinants << std::endl;
  OutputControl::nOut << std::endl;
  printTableHead("Eigenvalues/a.u.");
  OutputControl::nOut << std::fixed << std::setprecision(5) << E << std::endl;
  if (nStatesCouple == 2 && analyticCoupling) {
    printSubSectionTitle("Analytical Results for Coupling of 2 States");
    OutputControl::nOut << std::string(100, '-') << std::endl;
    OutputControl::n.printf("%-80s %18.10f \n", "Analytical electronic coupling/a.u.:", coupling);
    OutputControl::n.printf("%-80s %18.10f \n", "Analytical electronic coupling/eV:", coupling * HARTREE_TO_EV);
    OutputControl::n.printf("%-80s %18.10f \n", "Long-range excitation-energy/a.u.:", lree);
    OutputControl::n.printf("%-80s %18.10f \n", "Long-range excitation-energy/eV:", lree * HARTREE_TO_EV);
    OutputControl::nOut << std::string(100, '-') << std::endl;
  }
}

void FDEETCalculator::printTransformedMatrices(Eigen::MatrixXd HPrime, Eigen::MatrixXd S, Eigen::MatrixXd C) {
  printSubSectionTitle("Results for Coupling Problem in Transformed Basis");
  printTableHead("Normalized Linear Combination Coefficients for Construction of Adiabatic States");
  OutputControl::nOut << std::string(10, ' ');
  for (unsigned col = 0; col < C.cols(); ++col) {
    OutputControl::nOut << "Psi_" << col << std::string(5, ' ');
  }
  OutputControl::nOut.flush();
  OutputControl::nOut << std::endl;
  for (unsigned row = 0; row < C.rows(); ++row) {
    OutputControl::nOut << "Phi_" << row << std::string(2, ' ');
    OutputControl::nOut.flush();
    for (unsigned col = 0; col < C.cols(); ++col) {
      OutputControl::n.printf("%10.5f", C(row, col));
    }
    OutputControl::nOut << std::endl;
  }
  OutputControl::nOut
      << std::endl
      << "Note: Psi_i are the new adiabatic states, where Phi_i are the coupled adiabatic states used for the linear "
         "combination.\n"
      << std::endl;
  printTableHead("Transformed Hamiltonian/a.u.");
  OutputControl::nOut << std::fixed << std::setprecision(5) << HPrime << std::endl;
  OutputControl::nOut << std::endl
                      << "Note: The element-wise product of energy contributions and overlap is shown.\n"
                      << std::endl;
  printTableHead("Overlap of Coupled Adiabatic States/a.u.");
  OutputControl::nOut << std::fixed << std::setprecision(5) << S << std::endl;
}

void FDEETCalculator::printTransDensContributions(std::shared_ptr<SystemController> superSystem, unsigned nStatesCouple) {
  printSubSectionTitle("Writing Transition Densities to Cube");
  struct PlotTaskSettings plotSettings;
  plotSettings.gridSpacing = {0.12, 0.12, 0.12};
  auto cubeWriter = CubeFileWriter(superSystem->getSettings(), plotSettings);
  for (unsigned iState = 0; iState < nStatesCouple; ++iState) {
    for (unsigned jState = 0; jState < nStatesCouple; ++jState) {
      auto dMatController =
          std::make_shared<DensityMatrixController<Options::SCF_MODES::UNRESTRICTED>>(_densMats[iState][jState]);
      auto matrix = dMatController->getDensityMatrix();
      unsigned iSpin = 0;
      std::vector<std::string> labels{"alpha", "beta"};
      for_spin(matrix) {
        print((std::string) "Printing density matrix to file:" + superSystem->getSystemPath() + "diabDens" + iState +
              jState + "_" + labels[iSpin]);
        cubeWriter.writeMatrixToGrid(superSystem->getSystemPath() + "diabDens" + iState + jState + "_" + labels[iSpin],
                                     superSystem->getGeometry(), superSystem->getBasisController(), matrix_spin);
        ++iSpin;
      };
    }
  }
}

} /* namespace Serenity */
