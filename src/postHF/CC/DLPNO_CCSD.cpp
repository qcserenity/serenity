/**
 * @file DLPNO_CCSD.cpp
 *
 * @date May 16, 2019
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
#include "postHF/CC/DLPNO_CCSD.h" //Associated header file.
/* Include Serenity Internal Headers */
#include "analysis/PAOSelection/PNOConstructor.h"                   //PNO construction.
#include "data/OrbitalPair.h"                                       //Definition of an OrbitalPair.
#include "data/OrbitalPairSet.h"                                    //OrbitalPairSet definition.
#include "data/SingleSubstitution.h"                                //Definition of a SingleSubstitution.
#include "integrals/transformer/Ao2MoExchangeIntegralTransformer.h" //Integral calculation.
#include "integrals/wrappers/Libint.h"                              //Keep engines for Sigma vector construction.
#include "io/FormattedOutput.h"                                     //Captions etc.
#include "io/FormattedOutputStream.h"                               //Filtered output.
#include "io/HDF5.h"                                                //Read from disk.
#include "memory/MemoryManager.h"                                   //Memory handling
#include "misc/WarningTracker.h"                                    //Warnings.
#include "postHF/LocalCorrelation/CouplingOrbitalSet.h"             //Definition of a kSet.
#include "postHF/LocalCorrelation/DomainOverlapMatrixController.h"  //Overlap matrices between domains.
#include "postHF/LocalCorrelation/KLOrbitalSet.h"                   //Definition of a kl-set.
#include "postHF/LocalCorrelation/LocalAOSigmaVectorWrapper.h"      //Sigma vector construction
#include "postHF/LocalCorrelation/LocalCorrelationController.h"     //Various resources.
#include "postHF/LocalCorrelation/OrbitalPairDIISWrapper.h"         //DIIS
#include "system/SystemController.h"
/* Includes for direct testing of canonical vs. local amplitudes. */
#include "data/ElectronicStructure.h"                    //Test purpose only.
#include "data/OrbitalController.h"                      //Test purpose only.
#include "data/PAOController.h"                          //Test purpose only.
#include "integrals/OneElectronIntegralController.h"     //Test purpose only.
#include "integrals/looper/TwoElecFourCenterIntLooper.h" //Test purpose only.
#include "integrals/transformer/Ao2MoTransformer.h"      //Test purpose only.
/* Include Std and External Headers */
#include <iomanip> // std::setprecision/std::fixed

namespace Serenity {

DLPNO_CCSD::DLPNO_CCSD(std::shared_ptr<LocalCorrelationController> localCorrelationController, double maxResidual,
                       unsigned int maxCycles)
  : _localCorrelationController(localCorrelationController),
    _maxResidual(maxResidual),
    _maxCycles(maxCycles),
    _linearScalingSigmaVector(_localCorrelationController->getSettings().linearScalingSigmaVector) {
  if (!_linearScalingSigmaVector) {
    _g_ao_ao = std::make_shared<MatrixInBasis<RESTRICTED>>(
        _localCorrelationController->getActiveSystemController()->getBasisController());
    _g_ao_ao->setZero();
    unsigned int nOcc = _localCorrelationController->getActiveSystemController()->getNOccupiedOrbitals<RESTRICTED>();
    _g_ao_occ = Eigen::MatrixXd::Zero(_g_ao_ao->rows(), nOcc);
    _g_occ_occ = Eigen::MatrixXd::Zero(nOcc, nOcc);
  }
}

Eigen::VectorXd DLPNO_CCSD::calculateElectronicEnergyCorrections() {
  prepareOrbitalPairs();
  calculateCCSDIntegrals();
  if (testRun)
    switchIntegrals();
  optimizeAmplitudes();
  Eigen::VectorXd energies = calculateEnergyCorrection();
  deleteIntegralFiles();
  deleteIntegrals();
  return energies;
}

void DLPNO_CCSD::deleteIntegralFiles() {
  if (!_localCorrelationController->getSettings().dumpIntegrals &&
      _localCorrelationController->getCloseOrbitalPairSets().size() > 1) {
    std::remove(_localCorrelationController->getPairIntegralFileName().c_str());
  }
}

void DLPNO_CCSD::deleteIntegrals() {
  if (!dlpnoCCSDSettings.keepPairIntegrals) {
    OutputControl::dOut << "Deleting pair integrals" << std::endl;
    auto orbitalPairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
    for (auto& pair : orbitalPairs) {
      pair->cleanUp();
    }
  }
}

inline void DLPNO_CCSD::calculate_tau_ij(std::shared_ptr<OrbitalPair> pair) {
  pair->tau_ij = pair->t_ij + pair->getS_ij_i() * pair->singles_i->t_i * pair->singles_j->t_i.transpose() *
                                  pair->getS_ij_j().transpose();
}

inline void DLPNO_CCSD::calculate_Z(std::shared_ptr<OrbitalPair> pair) {
  Eigen::MatrixXd z_ij = Eigen::MatrixXd::Zero(pair->t_ij.cols(), pair->t_ij.cols());
  Eigen::MatrixXd z_ji = Eigen::MatrixXd::Zero(pair->t_ij.cols(), pair->t_ij.cols());

  const Eigen::VectorXd& t_j = pair->singles_j->t_i;
  const Eigen::VectorXd& t_i = pair->singles_i->t_i;
  for (unsigned int cPNO = 0; cPNO < t_j.size(); ++cPNO)
    z_ij += t_j[cPNO] * pair->ic_ab[cPNO];
  for (unsigned int cPNO = 0; cPNO < t_i.size(); ++cPNO)
    z_ji += t_i[cPNO] * pair->jc_ab[cPNO];

  for (const auto& kSet : pair->coupledPairs) {
    const auto& ikPair = kSet->getIKPair();
    const auto& kjPair = kSet->getKJPair();
    const Eigen::MatrixXd k_ki = kSet->getS_ij_ik() * ((pair->i == ikPair->j) ? ikPair->k_ij : ikPair->k_ij.transpose());
    const Eigen::MatrixXd k_kj = kSet->getS_ij_kj() * ((pair->j == kjPair->j) ? kjPair->k_ij : kjPair->k_ij.transpose());
    const Eigen::MatrixXd t_ik = (pair->i == ikPair->i) ? ikPair->t_ij : ikPair->t_ij.transpose();
    const Eigen::MatrixXd t_jk = (pair->j == kjPair->i) ? kjPair->t_ij : kjPair->t_ij.transpose();
    const Eigen::MatrixXd& s_ik_kj = kSet->getS_ik_kj();
    const Eigen::MatrixXd& s_kj_ik = kSet->getS_ik_kj().transpose();
    z_ij -= 0.5 * k_ki * s_ik_kj * t_jk * kSet->getS_ij_kj().transpose();
    z_ji -= 0.5 * k_kj * s_kj_ik * t_ik * kSet->getS_ij_ik().transpose();

    const Eigen::VectorXd& t_k = kSet->getKSingles()->t_i;
    const Eigen::MatrixXd& s_ik_j = kSet->getS_ik_j();
    const Eigen::MatrixXd& s_kj_i = kSet->getS_kj_i();
    const Eigen::MatrixXd& s_ij_k = kSet->getS_ij_k();
    z_ij -= (kSet->ij_ak + k_ki * s_ik_j * t_j) * t_k.transpose() * s_ij_k.transpose();
    z_ji -= (kSet->ij_ak + k_kj * s_kj_i * t_i) * t_k.transpose() * s_ij_k.transpose();
  } // for kSet
  pair->z_ij = z_ij;
  pair->z_ji = z_ji;
}

inline void DLPNO_CCSD::calculate_Y(std::shared_ptr<OrbitalPair> pair) {
  Eigen::MatrixXd y_ij = Eigen::MatrixXd::Zero(pair->t_ij.cols(), pair->t_ij.cols());
  Eigen::MatrixXd y_ji = Eigen::MatrixXd::Zero(pair->t_ij.cols(), pair->t_ij.cols());

  const Eigen::VectorXd& t_j = pair->singles_j->t_i;
  const Eigen::VectorXd& t_i = pair->singles_i->t_i;
  for (unsigned int cPNO = 0; cPNO < t_j.size(); ++cPNO)
    y_ij += t_j[cPNO] * (pair->ia_bc[cPNO] - 0.5 * pair->ic_ab[cPNO]);
  for (unsigned int cPNO = 0; cPNO < t_i.size(); ++cPNO)
    y_ji += t_i[cPNO] * (pair->ja_bc[cPNO] - 0.5 * pair->jc_ab[cPNO]);

  for (const auto& kSet : pair->coupledPairs) {
    const auto& ikPair = kSet->getIKPair();
    const auto& kjPair = kSet->getKJPair();
    const Eigen::MatrixXd k_ik = (pair->i == ikPair->i) ? ikPair->k_ij : ikPair->k_ij.transpose();
    const Eigen::MatrixXd k_jk = (pair->j == kjPair->i) ? kjPair->k_ij : kjPair->k_ij.transpose();
    const Eigen::MatrixXd t_ki = (pair->i == ikPair->j) ? ikPair->t_ij : ikPair->t_ij.transpose();
    const Eigen::MatrixXd t_kj = (pair->j == kjPair->j) ? kjPair->t_ij : kjPair->t_ij.transpose();
    const Eigen::MatrixXd& s_ik_kj = kSet->getS_ik_kj();
    const Eigen::MatrixXd s_kj_ik = kSet->getS_ik_kj().transpose();
    const Eigen::MatrixXd L_ik = kSet->getS_ij_ik() * (2.0 * k_ik - k_ik.transpose());
    const Eigen::MatrixXd L_jk = kSet->getS_ij_kj() * (2.0 * k_jk - k_jk.transpose());
    y_ij += 0.25 * L_ik * s_ik_kj * (2.0 * t_kj - t_kj.transpose()) * kSet->getS_ij_kj().transpose();
    y_ji += 0.25 * L_jk * s_kj_ik * (2.0 * t_ki - t_ki.transpose()) * kSet->getS_ij_ik().transpose();

    const Eigen::VectorXd& t_k = kSet->getKSingles()->t_i;
    const Eigen::MatrixXd& s_ik_j = kSet->getS_ik_j();
    const Eigen::MatrixXd& s_kj_i = kSet->getS_kj_i();
    const Eigen::MatrixXd& s_ij_k = kSet->getS_ij_k();
    y_ij -= 0.5 * ((2.0 * kSet->ia_jk - kSet->ij_ak) + (L_ik * s_ik_j * t_j)) * t_k.transpose() * s_ij_k.transpose();
    y_ji -= 0.5 * ((2.0 * kSet->ja_ik - kSet->ij_ak) + (L_jk * s_kj_i * t_i)) * t_k.transpose() * s_ij_k.transpose();
  } // for kSet
  pair->y_ij = y_ij;
  pair->y_ji = y_ji;
}

inline void DLPNO_CCSD::calculate_tilde_ikjl(std::shared_ptr<OrbitalPair> pair) {
  for (auto& klPairSet : pair->klPairSets) {
    const auto klPair = klPairSet->getKLPair();
    const Eigen::MatrixXd& s_ij_kl = klPairSet->getS_ij_kl();
    klPairSet->tilde_ik_jl =
        klPairSet->ik_jl + (pair->tau_ij * s_ij_kl * klPair->k_ij.transpose() * s_ij_kl.transpose()).trace() +
        klPairSet->ki_la.transpose() * pair->singles_j->t_i + klPairSet->lj_ka.transpose() * pair->singles_i->t_i;
    klPairSet->tilde_il_jk = klPairSet->il_jk + (pair->tau_ij * s_ij_kl * klPair->k_ij * s_ij_kl.transpose()).trace() +
                             klPairSet->li_ka.transpose() * pair->singles_j->t_i +
                             klPairSet->kj_la.transpose() * pair->singles_i->t_i;
  } // for klPairSet
}

inline Eigen::MatrixXd DLPNO_CCSD::calculate_K_tau(std::shared_ptr<OrbitalPair> pair) {
  assert(pair->tau_ij.rows() != 0);
  const Matrix<Eigen::MatrixXd>& ac_bd = *pair->ac_bd;
  unsigned int nPNOs = ac_bd.rows();
  Eigen::MatrixXd resultTest = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
  for (unsigned int a = 0; a < nPNOs; ++a) {
    for (unsigned int b = 0; b < nPNOs; ++b) {
      Eigen::MatrixXd tmp = (a <= b) ? ac_bd(a, b) : ac_bd(b, a).transpose();
      for (const auto& kSet : pair->coupledPairs) {
        const std::shared_ptr<SingleSubstitution>& kSingle = kSet->getKSingles();
        const Eigen::MatrixXd kd_ac = kSet->ka_bc[a].transpose();
        const Eigen::MatrixXd kc_bd = kSet->ka_bc[b];
        const Eigen::VectorXd t_kinIJ = kSet->getS_ij_k() * kSingle->t_i;
        tmp -= kd_ac * t_kinIJ(b);
        tmp -= kc_bd * t_kinIJ(a);
      }
      resultTest(a, b) += (tmp.array() * pair->tau_ij.array()).sum();
    }
  }
  return resultTest;
}

inline void DLPNO_CCSD::calculate_tildeF_ij(std::shared_ptr<OrbitalPair> pair) {
  double tildeF_ij = pair->f_ij;
  double tildeF_ji = pair->f_ij;
  unsigned int i = pair->i;
  unsigned int j = pair->j;
  for (const auto& coupledPair : pair->coupledPairs) {
    auto ikPair = coupledPair->getIKPair();
    auto kjPair = coupledPair->getKJPair();
    const Eigen::MatrixXd tau_ik = (ikPair->i == i) ? ikPair->tau_ij : ikPair->tau_ij.transpose();
    const Eigen::MatrixXd tau_jk = (kjPair->i == j) ? kjPair->tau_ij : kjPair->tau_ij.transpose();
    const Eigen::MatrixXd t_ik = (ikPair->i == i) ? ikPair->t_ij : ikPair->t_ij.transpose();
    const Eigen::MatrixXd t_jk = (kjPair->i == j) ? kjPair->t_ij : kjPair->t_ij.transpose();
    const Eigen::VectorXd& t_k = coupledPair->getKSingles()->t_i;
    const Eigen::VectorXd& t_i = pair->singles_i->t_i;
    const Eigen::VectorXd& t_j = pair->singles_j->t_i;
    const Eigen::MatrixXd& s_kj_i = coupledPair->getS_kj_i();
    const Eigen::MatrixXd& s_ik_j = coupledPair->getS_ik_j();
    const Eigen::MatrixXd& s_kj_k = coupledPair->getS_kj_k();
    const Eigen::MatrixXd& s_ik_k = coupledPair->getS_ik_k();

    const Eigen::MatrixXd& s_kj_ik = coupledPair->getS_kj_ik();
    const Eigen::MatrixXd k_jk = (kjPair->i == j) ? kjPair->k_ij : kjPair->k_ij.transpose();
    const Eigen::MatrixXd k_ik = (ikPair->i == i) ? ikPair->k_ij : ikPair->k_ij.transpose();
    // Introduced projection to save computational time. This term contributes in second order
    // in the amplitudes.
    const Eigen::MatrixXd L_ki = 2.0 * k_ik.transpose() - k_ik;
    const Eigen::MatrixXd L_kj = 2.0 * k_jk.transpose() - k_jk;
    tildeF_ij +=
        (t_ik * s_kj_ik.transpose() * L_kj * s_kj_ik).trace() + t_k.transpose() * s_kj_k.transpose() * L_kj * s_kj_i * t_i;
    tildeF_ji +=
        (t_jk * s_kj_ik * L_ki * s_kj_ik.transpose()).trace() + t_k.transpose() * s_ik_k.transpose() * L_ki * s_ik_j * t_j;
  } // for kIndex
  pair->tf_ij = tildeF_ij;
  pair->tf_ji = tildeF_ji;
}

inline void DLPNO_CCSD::calculate_tildeF_ab(std::shared_ptr<OrbitalPair> pair) {
  Eigen::MatrixXd tildeF_ab = pair->f_ab.asDiagonal();
  Eigen::MatrixXd sumOverKL = Eigen::MatrixXd::Zero(tildeF_ab.cols(), tildeF_ab.cols());
  for (const auto& klPairSet : pair->klPairSets) {
    const auto& klPair = klPairSet->getKLPair();
    const Eigen::MatrixXd& s_ij_kl = klPairSet->getS_ij_kl();
    const Eigen::MatrixXd& tau_kl = klPair->tau_ij;
    Eigen::MatrixXd A = (2.0 * klPair->k_ij.transpose() - klPair->k_ij) * tau_kl;
    if (klPair->i != klPair->j)
      A += (2.0 * klPair->k_ij - klPair->k_ij.transpose()) * tau_kl.transpose();
    sumOverKL -= s_ij_kl * A * s_ij_kl.transpose();
  } // for kIndex
  tildeF_ab += sumOverKL;
  if (pair->i == pair->j) {
    auto& singles = pair->singles_i;
    singles->tf_ab =
        (Eigen::MatrixXd)singles->f_ab.asDiagonal() + pair->getS_ij_i().transpose() * sumOverKL * pair->getS_ij_i();
  }
  pair->tf_ab = tildeF_ab;
}

inline void DLPNO_CCSD::calculate_tildeF_ai(std::shared_ptr<SingleSubstitution> single) {
  unsigned int i = single->i;
  Eigen::VectorXd tildeF_ai = single->f_ai;
  for (const auto& ikPair_ptr : single->orbitalPairs) {
    const auto ikPair = ikPair_ptr.lock();
    bool i_is_i = ikPair->i == i;
    const std::shared_ptr<SingleSubstitution> single_k = (i_is_i) ? ikPair->singles_j : ikPair->singles_i;
    const Eigen::MatrixXd& s_ik_k = (i_is_i) ? ikPair->getS_ij_j() : ikPair->getS_ij_i();
    const Eigen::MatrixXd& s_ik_i = (i_is_i) ? ikPair->getS_ij_i() : ikPair->getS_ij_j();
    const Eigen::MatrixXd k_ik = (i_is_i) ? ikPair->k_ij : ikPair->k_ij.transpose();
    tildeF_ai += s_ik_i.transpose() * (2.0 * k_ik - k_ik.transpose()) * s_ik_k * single_k->t_i;
  } // for k
  single->tf_ai = tildeF_ai;
}

inline std::pair<double, double> DLPNO_CCSD::calculate_G_t1_ij(std::shared_ptr<OrbitalPair> pair) {
  double G_t1_ij = 0.0;
  double G_t1_ji = 0.0;
  if (!_linearScalingSigmaVector) {
    G_t1_ij = _g_occ_occ(pair->i, pair->j);
    G_t1_ji = _g_occ_occ(pair->j, pair->i);
  }
  else if (not testRun) {
    for (const auto& kSet : pair->coupledPairs) {
      const auto& singles_k = kSet->getKSingles();
      G_t1_ij += kSet->ij_akX2_M_ia_jk.transpose() * singles_k->t_i;
      G_t1_ji += kSet->ij_akX2_M_ja_ik.transpose() * singles_k->t_i;
    } // for coupledSet
  }
  // The following lines are for extremely strict tests only!
  if (testRun) {
    double tmpGt1_ij = G_t1_ij;
    double tmpGt1_ji = G_t1_ji;
    G_t1_ij = 0.0;
    G_t1_ji = 0.0;
    for (const auto& coupledSet : pair->coupledPairs) {
      const auto& singles_k = coupledSet->getKSingles();
      const Eigen::VectorXd tmp = coupledSet->getS_ij_k() * singles_k->t_i;
      G_t1_ij += (2.0 * coupledSet->ij_ak - coupledSet->ja_ik).transpose() * tmp;
      G_t1_ji += (2.0 * coupledSet->ij_ak - coupledSet->ia_jk).transpose() * tmp;
    } // for coupledSet
    if (std::fabs(tmpGt1_ij - G_t1_ij) > 1e-9)
      OutputControl::dOut << " G(t1)_ij Check: " << tmpGt1_ij - G_t1_ij << std::endl;
    if (std::fabs(tmpGt1_ji - G_t1_ji) > 1e-9)
      OutputControl::dOut << " G(t1)_ji Check: " << tmpGt1_ji - G_t1_ji << std::endl;
  } // if testRun
  return std::pair<double, double>(G_t1_ij, G_t1_ji);
}

inline Eigen::VectorXd DLPNO_CCSD::calculate_G_t1_ia(std::shared_ptr<SingleSubstitution> single) {
  Eigen::VectorXd G_t1_i;
  if (!_linearScalingSigmaVector) {
    Eigen::MatrixXd p_ii = (Eigen::MatrixXd)_localCorrelationController->getPAOController()->getAllPAOs() *
                           single->getDiagonalPair()->domainProjection.transpose() * single->toPAODomain;
    if (testRun) {
      p_ii = _localCorrelationController->getActiveSystemController()
                 ->getActiveOrbitalController<RESTRICTED>()
                 ->getCoefficients()
                 .rightCols(single->t_i.rows())
                 .eval();
    }
    G_t1_i = (p_ii.transpose() * _g_ao_occ.col(single->i)).eval();
    // The following lines are for extremely strict tests only!
  }
  else {
    G_t1_i = Eigen::VectorXd::Zero(single->t_i.size());
    unsigned int i = single->i;
    for (const auto& ijPair_ptr : single->orbitalPairs) {
      const auto ijPair = ijPair_ptr.lock();
      bool i_is_i = ijPair->i == i;
      std::shared_ptr<SingleSubstitution> single_j = (i_is_i) ? ijPair->singles_j : ijPair->singles_i;
      G_t1_i += ((i_is_i) ? ijPair->iaS_jbSX2_M_ij_aSbS : ijPair->iaS_jbSX2_M_ij_aSbS.transpose()) * single_j->t_i;
    } // for ijPair
  }
  if (testRun) {
    auto tmpTest = G_t1_i;
    G_t1_i = Eigen::VectorXd::Zero(single->t_i.size());
    unsigned int i = single->i;
    for (const auto& ijPair_ptr : single->orbitalPairs) {
      const auto ijPair = ijPair_ptr.lock();
      bool i_is_i = ijPair->i == i;
      std::shared_ptr<SingleSubstitution> single_j = (i_is_i) ? ijPair->singles_j : ijPair->singles_i;
      G_t1_i += ((i_is_i) ? ijPair->iaS_jbSX2_M_ij_aSbS : ijPair->iaS_jbSX2_M_ij_aSbS.transpose()) * single_j->t_i;
    } // for ijPair
    if ((tmpTest - G_t1_i).array().abs().sum() > 1e-9)
      OutputControl::dOut << "G(t1)_ia Check: " << (tmpTest - G_t1_i).array().abs().sum() << std::endl;
  } // if testRun
  return G_t1_i;
}

inline Eigen::MatrixXd DLPNO_CCSD::calculate_G_t1_ab(std::shared_ptr<OrbitalPair> pair) {
  Eigen::MatrixXd G_t1_ab;
  if (!_linearScalingSigmaVector) {
    Eigen::MatrixXd p_ij = (Eigen::MatrixXd)_localCorrelationController->getPAOController()->getAllPAOs() *
                           pair->domainProjection.transpose() * pair->toPAODomain;
    if (testRun) {
      p_ij = _localCorrelationController->getActiveSystemController()
                 ->getActiveOrbitalController<RESTRICTED>()
                 ->getCoefficients()
                 .rightCols(pair->t_ij.rows())
                 .eval();
    }
    G_t1_ab = (p_ij.transpose() * *_g_ao_ao * p_ij).eval();
  }
  else {
    unsigned int nPNOs = pair->t_ij.cols();
    G_t1_ab = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
    for (const auto& coupledSet : pair->coupledPairs) {
      const Eigen::VectorXd& t_k = coupledSet->getKSingles()->t_i;
      const std::vector<Eigen::MatrixXd>& ints = coupledSet->ab_kcX2_M_ak_bc;
      for (unsigned int c = 0; c < t_k.size(); ++c) {
        G_t1_ab += t_k(c) * ints[c];
      } // for c
    }   // for k
  }
  // The following lines are for extremely strict tests only!
  if (testRun) {
    auto tmpTest = G_t1_ab;
    unsigned int nPNOs = pair->t_ij.cols();
    G_t1_ab = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
    for (const auto& coupledSet : pair->coupledPairs) {
      const Eigen::VectorXd& t_k = coupledSet->getKSingles()->t_i;
      const std::vector<Eigen::MatrixXd>& ints = coupledSet->ab_kcX2_M_ak_bc;
      for (unsigned int c = 0; c < t_k.size(); ++c) {
        G_t1_ab += t_k(c) * ints[c];
      } // for c
    }   // for k
    if ((tmpTest - G_t1_ab).array().abs().sum() > 1e-9)
      OutputControl::dOut << " G(t1)_ab Check: " << (tmpTest - G_t1_ab).array().abs().sum() << std::endl;
  } // if testRun
  return G_t1_ab;
}

inline void DLPNO_CCSD::calculate_tildetildeF_ij(std::shared_ptr<OrbitalPair> pair) {
  double ttF_ij = pair->tf_ij;
  double ttF_ji = pair->tf_ji;
  const auto G_t1_ijANDji = calculate_G_t1_ij(pair);
  ttF_ij += G_t1_ijANDji.first;
  ttF_ji += G_t1_ijANDji.second;
  ttF_ij += pair->singles_i->f_ai.transpose() * pair->getS_i_j() * pair->singles_j->t_i;
  ttF_ji += pair->singles_j->f_ai.transpose() * pair->getS_i_j().transpose() * pair->singles_i->t_i;
  pair->ttf_ij = ttF_ij;
  pair->ttf_ji = ttF_ji;
}

inline Eigen::MatrixXd DLPNO_CCSD::calculate_tildetildeF_ab(std::shared_ptr<OrbitalPair> pair) {
  Eigen::MatrixXd ttF_ab = pair->tf_ab;
  ttF_ab += calculate_G_t1_ab(pair);
  for (const auto& coupledSet : pair->coupledPairs) {
    const Eigen::MatrixXd& s_ij_k = coupledSet->getS_ij_k();
    const std::shared_ptr<SingleSubstitution>& single_k = coupledSet->getKSingles();
    ttF_ab -= s_ij_k * single_k->t_i * single_k->f_ai.transpose() * s_ij_k.transpose();
  } // for k
  return ttF_ab;
}

inline Eigen::VectorXd DLPNO_CCSD::calculateSinglesResidual(std::shared_ptr<SingleSubstitution> single) {
  Eigen::VectorXd residual = single->f_ai;
  // sum_b tF_ba*t_b^i
  residual += single->tf_ab.transpose() * single->t_i;
  // G(t_1)_ia
  residual += calculate_G_t1_ia(single);

  for (const auto& ijPair_ptr : single->orbitalPairs) {
    const auto ijPair = ijPair_ptr.lock();
    bool i_is_i = ijPair->i == single->i;
    const Eigen::MatrixXd& s_ij_j = (i_is_i) ? ijPair->getS_ij_j() : ijPair->getS_ij_i();
    const Eigen::MatrixXd s_i_j = (i_is_i) ? ijPair->getS_i_j() : ijPair->getS_i_j().transpose();
    const std::shared_ptr<SingleSubstitution>& single_j = (i_is_i) ? ijPair->singles_j : ijPair->singles_i;
    // tF_ij*t_a^j
    residual -= ((i_is_i) ? ijPair->tf_ij : ijPair->tf_ji) * s_i_j * single_j->t_i;
    const Eigen::MatrixXd t_ij = (i_is_i) ? ijPair->t_ij : ijPair->t_ij.transpose();
    const Eigen::MatrixXd t_ji = (i_is_i) ? ijPair->t_ij.transpose() : ijPair->t_ij;
    const Eigen::MatrixXd s_i_ij = (i_is_i) ? ijPair->getS_ij_i().transpose() : ijPair->getS_ij_j().transpose();
    // sum_b (2t^ji-t^ji^T)*tF_bj
    residual += s_i_ij * (2.0 * t_ij - t_ij.transpose()) * s_ij_j * single_j->tf_ai;
    // sum_kb (2(ik|jb)-(ij|kb))tau_ab^kj
    for (const auto& kSet : ijPair->coupledPairs) {
      const unsigned int k = kSet->getK();
      const std::shared_ptr<OrbitalPair>& kjPair = (i_is_i) ? kSet->getKJPair() : kSet->getIKPair();
      const Eigen::MatrixXd tau_kj = (kjPair->i == k) ? kjPair->tau_ij : kjPair->tau_ij.transpose();
      const Eigen::MatrixXd& s_ij_kj = (i_is_i) ? kSet->getS_ij_kj() : kSet->getS_ij_ik();
      const Eigen::VectorXd& ik_jb = (i_is_i) ? kSet->ja_ik : kSet->ia_jk;
      const Eigen::MatrixXd& s_kj_i = (i_is_i) ? kSet->getS_kj_i() : kSet->getS_ik_j();
      residual -= s_kj_i.transpose() * tau_kj * s_ij_kj.transpose() * (2.0 * ik_jb - kSet->ij_ak);
    } // for kSet

    const Eigen::MatrixXd tau_ij = (i_is_i) ? ijPair->tau_ij : ijPair->tau_ij.transpose();
    const std::vector<Eigen::MatrixXd>& ja_bc = (i_is_i) ? ijPair->ja_bc : ijPair->ia_bc;
    for (unsigned int a = 0; a < residual.size(); ++a) {
      // sum_jbc (jc|ab)(2tau_ij-tau_ji)
      //[a] determines one of the indices for coordinate 2.
      //-->(ib|ac) is given by ia_bc[a], since b and c define the two remaining indices.
      residual(a) += (ja_bc[a].transpose().array() * (2.0 * tau_ij - tau_ij.transpose()).array()).sum();
    } // for a
    // p = sum_b (tF_jb - 2F_jb)t_b^i
    double prefactor = ((single_j->tf_ai - 2.0 * single_j->f_ai).array() * (s_i_j.transpose() * single->t_i).array()).sum();
    // p * t_a^j
    residual += prefactor * s_i_j * single_j->t_i;
  } // for ijPair
  return residual;
}

inline Eigen::MatrixXd DLPNO_CCSD::calculateDoublesResidual(std::shared_ptr<OrbitalPair> pair) {
  Eigen::MatrixXd residual = pair->k_ij;
  residual += calculate_K_tau(pair);
  const Eigen::MatrixXd ttF_ab = calculate_tildetildeF_ab(pair);
  residual += ttF_ab.transpose() * pair->t_ij + pair->t_ij * ttF_ab;
  for (const auto& kSet : pair->coupledPairs) {
    const std::shared_ptr<OrbitalPair>& ikPair = kSet->getIKPair();
    const std::shared_ptr<OrbitalPair>& kjPair = kSet->getKJPair();
    const std::shared_ptr<SingleSubstitution>& singles_k = kSet->getKSingles();
    const Eigen::MatrixXd& s_ij_ik = kSet->getS_ij_ik();
    const Eigen::MatrixXd& s_ij_kj = kSet->getS_ij_kj();
    const Eigen::MatrixXd s_ik_ij = s_ij_ik.transpose();
    const Eigen::MatrixXd s_kj_ij = s_ij_kj.transpose();
    const Eigen::MatrixXd s_ik_kj = kSet->getS_ik_kj();
    const Eigen::MatrixXd s_jk_ik = s_ik_kj.transpose();
    const Eigen::MatrixXd s_k_ij = kSet->getS_ij_k().transpose();
    const Eigen::MatrixXd s_j_k = kSet->getS_j_k();
    const Eigen::MatrixXd s_i_k = kSet->getS_i_k();
    const Eigen::MatrixXd t_ik = (pair->i == ikPair->i) ? ikPair->t_ij : ikPair->t_ij.transpose();
    const Eigen::MatrixXd t_kj = (pair->j == kjPair->j) ? kjPair->t_ij : kjPair->t_ij.transpose();
    double ttF_jk = (kjPair->i == pair->j) ? kjPair->ttf_ij : kjPair->ttf_ji;
    double ttF_ik = (ikPair->i == pair->i) ? ikPair->ttf_ij : ikPair->ttf_ji;
    // sum_k (ttF_jk * t^kj^T + ttF_ik + t^kj)_ab
    residual -= ttF_jk * s_ij_ik * t_ik * s_ij_ik.transpose();
    residual -= ttF_ik * s_ij_kj * t_kj * s_ij_kj.transpose();
    // Semi-joint pair--pair interaction term.
    const Eigen::MatrixXd y_kj = (kjPair->j == pair->j) ? kjPair->y_ij : kjPair->y_ji;
    const Eigen::MatrixXd y_ki = (ikPair->j == pair->i) ? ikPair->y_ij : ikPair->y_ji;
    const Eigen::MatrixXd fullY_kj = kSet->ja_kc.transpose() - 0.5 * kSet->jk_ca + s_ik_kj * y_kj * s_ij_kj.transpose(); //[ik]x[ij]
    const Eigen::MatrixXd fullY_ki =
        kSet->ia_kc.transpose() - 0.5 * kSet->ik_ca + s_ik_kj.transpose() * y_ki * s_ij_ik.transpose(); //[kj]x[ij]
    const Eigen::MatrixXd z_kj = (kjPair->j == pair->j) ? kjPair->z_ij : kjPair->z_ji;
    const Eigen::MatrixXd z_ki = (ikPair->j == pair->i) ? ikPair->z_ij : ikPair->z_ji;
    const Eigen::MatrixXd fullZ_kj = kSet->jk_ca + s_ik_kj * z_kj * s_ij_kj.transpose();             //[ik]x[ij]
    const Eigen::MatrixXd fullZ_ki = kSet->ik_ca + s_ik_kj.transpose() * z_ki * s_ij_ik.transpose(); //[kj]x[ij]
    residual += s_ij_ik * (2.0 * t_ik - t_ik.transpose()) * fullY_kj +
                fullY_ki.transpose() * (2.0 * t_kj - t_kj.transpose()) * s_ij_kj.transpose();
    residual +=
        -0.5 * (s_ij_ik * t_ik.transpose() * fullZ_kj + fullZ_ki.transpose() * t_kj.transpose() * s_ij_kj.transpose());
    residual += -(s_ij_kj * t_kj * fullZ_ki + fullZ_kj.transpose() * t_ik * s_ij_ik.transpose());

    //(jk|ia)t^k_b+(ik|jb)t_a^k
    const Eigen::MatrixXd t_k_T = singles_k->t_i.transpose();
    residual -= kSet->ia_jk * t_k_T * s_k_ij;
    residual -= kSet->getS_ij_k() * singles_k->t_i * kSet->ja_ik.transpose();
    // Last sum over k
    //(K^ik*t_j)_a * t^k_b
    const Eigen::MatrixXd k_ik = kSet->ia_kc;
    // This should be symmetric!
    const Eigen::MatrixXd& J_jk = kSet->jk_ca.transpose();
    const Eigen::VectorXd Kik_tj = k_ik * kSet->getS_kj_j() * pair->singles_j->t_i;
    const Eigen::VectorXd Jjk_ti = J_jk * kSet->getS_ik_i() * pair->singles_i->t_i;
    residual -= (Kik_tj + Jjk_ti) * t_k_T * s_k_ij;
    //(K^jk t^i)_b t^k_a + (J^ik t^j)_b t^k_a
    const Eigen::MatrixXd k_kj = kSet->ja_kc.transpose();
    const Eigen::MatrixXd& J_ik = kSet->ik_ca; //[kj] x [ij]
    // These are row vectors.
    const Eigen::MatrixXd ti_Kjk = pair->singles_i->t_i.transpose() * kSet->getS_ik_i().transpose() * k_kj;
    const Eigen::MatrixXd tj_Jik = pair->singles_j->t_i.transpose() * kSet->getS_kj_j().transpose() * J_ik;
    residual -= kSet->getS_ij_k() * singles_k->t_i * (ti_Kjk + tj_Jik);
  } // for kSet
  // sum_kl tilde(ik|jl)*tau^kl_ab
  calculate_tilde_ikjl(pair);
  for (const auto& klPairSet : pair->klPairSets) {
    const auto& klPair = klPairSet->getKLPair();
    const Eigen::MatrixXd& s_ij_kl = klPairSet->getS_ij_kl();
    const Eigen::MatrixXd& s_ij_k = klPairSet->getS_ij_k();
    const Eigen::MatrixXd& s_ij_l = klPairSet->getS_ij_l();
    const Eigen::VectorXd& t_k = klPair->singles_i->t_i;
    const Eigen::VectorXd& t_l = klPair->singles_j->t_i;
    residual += klPairSet->tilde_ik_jl * s_ij_kl * klPair->t_ij * s_ij_kl.transpose();
    residual += klPairSet->tilde_ik_jl * s_ij_k * t_k * t_l.transpose() * s_ij_l.transpose();
    if (klPair->i != klPair->j) {
      residual += klPairSet->tilde_il_jk * s_ij_kl * klPair->t_ij.transpose() * s_ij_kl.transpose();
      residual += klPairSet->tilde_il_jk * s_ij_l * t_l * t_k.transpose() * s_ij_k.transpose();
    }
  } // for klPairSet
  // sum_c (jb|ac)t_c^i
  //[c] defines of the indices for coordinate 2.
  //--> The set  (jb|ac) is given by ja_bc[c] where c is in [i]
  for (unsigned int c = 0; c < pair->ja_bc.size(); ++c) {
    residual += pair->ja_bc[c].transpose() * pair->singles_i->t_i[c];
  } // for c
  // sum_c (ia|bc)t_c^j
  for (unsigned int c = 0; c < pair->ia_bc.size(); ++c) {
    residual += pair->ia_bc[c] * pair->singles_j->t_i[c];
  } // for c
  return residual;
}

inline void DLPNO_CCSD::dressSingles() {
  std::vector<std::shared_ptr<SingleSubstitution>> singles = _localCorrelationController->getSingles();
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iSingle = 0; iSingle < singles.size(); ++iSingle) {
    std::shared_ptr<SingleSubstitution> single = singles[iSingle];
    calculate_tildeF_ai(single);
  } // for single
}

inline void DLPNO_CCSD::dressPairs() {
  std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs =
      _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    std::shared_ptr<OrbitalPair> pair = orbitalPairs[iPair];
    calculate_tau_ij(pair);
  } // for pair
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
    std::shared_ptr<OrbitalPair> pair = orbitalPairs[iPair];
    calculate_Y(pair);
    calculate_Z(pair);
    calculate_tildeF_ij(pair);
    calculate_tildetildeF_ij(pair);
    calculate_tildeF_ab(pair);
  } // for pair
}

inline void DLPNO_CCSD::prepareOrbitalPairs() {
  std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs =
      _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  unsigned int nStrongPlusWeakPairs = orbitalPairs.size();
  unsigned int nVeryWeakPairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::VERY_DISTANT).size();
  const auto activeSystem = _localCorrelationController->getActiveSystemController();
  const auto paoController = _localCorrelationController->getPAOController();

  auto qcPAOConstructor = _localCorrelationController->produceQCPAOConstructor();

  // Ensure that at least some pair energies were set before to avoid removing all orbital pairs
  // because they have a pair energy of 0.0 .
  bool pairEnergiesWereSetBefore = false;
  if (dlpnoCCSDSettings.skipCrudePrescreening) {
    for (const auto& pair : orbitalPairs) {
      if (pair->scMP2PairEnergy != 0.0 || pair->dlpnoCCSDPairEnergy != 0.0) {
        pairEnergiesWereSetBefore = true;
        break;
      }
    }
  }

  if (not dlpnoCCSDSettings.skipCrudePrescreening || not pairEnergiesWereSetBefore) {
    OutputControl::nOut << " Crude SC-MP2 guess" << std::endl;
    Ao2MoExchangeIntegralTransformer::transformExchangeIntegrals(
        activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL),
        _localCorrelationController->getApproximateMO3CenterIntegralController(), orbitalPairs, qcPAOConstructor);
    _localCorrelationController->removeApproximateMO3CenterIntegralController();
  }
  _localCorrelationController->selectDistantOrbitalPairs();

  std::shared_ptr<PNOConstructor> pnoConstructor = _localCorrelationController->producePNOConstructor();
  orbitalPairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  OutputControl::vOut << "Omitted pairs after crude screening: " << nStrongPlusWeakPairs - orbitalPairs.size() << std::endl;
  _localCorrelationController->buildSingles();
  OutputControl::nOut << " Accurate SC-MP2 guess" << std::endl;
  std::vector<std::shared_ptr<OrbitalPair>> toTransform = orbitalPairs;
  std::vector<std::shared_ptr<OrbitalPair>> distantTriplePairs =
      _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::DISTANT_TRIPLES);
  toTransform.insert(toTransform.end(), distantTriplePairs.begin(), distantTriplePairs.end());
  Ao2MoExchangeIntegralTransformer::transformExchangeIntegrals(
      activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL),
      _localCorrelationController->getMO3CenterIntegralController(false), toTransform, pnoConstructor);
  _localCorrelationController->selectDistantOrbitalPairs();
  _localCorrelationController->initializeSingles();
  orbitalPairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  auto singles = _localCorrelationController->getSingles();
  _localCorrelationController->buildOrbitalPairCouplingMap();
  _localCorrelationController->buildKLOrbitalPairs();
  OutputControl::nOut << "  Calculating overlap matrices                           ...";
  OutputControl::nOut.flush();
  auto domainOverlapController = _localCorrelationController->getDomainOverlapMatrixController();
  OutputControl::nOut << " done" << std::endl;
  for (auto pair : orbitalPairs) {
    pair->setOverlapMatrixController(domainOverlapController);
    for (auto kSet : pair->coupledPairs)
      kSet->setOverlapMatrixController(domainOverlapController);
    for (auto klSet : pair->klPairSets)
      klSet->setOverlapMatrixController(domainOverlapController);
  } // for pair
  if (testRun)
    return;
  // Print some information about the actual correlated system.
  unsigned int sumOfPNOs = 0;
  Eigen::VectorXi pnoSpread = Eigen::VectorXi::Zero(11);
  Eigen::VectorXi pnoTraceRecovery = Eigen::VectorXi::Zero(4);
  unsigned int maxNPNOs = 0;
  for (const auto& pair : orbitalPairs) {
    unsigned int nPNOs = pair->k_ij.cols();
    if (nPNOs > maxNPNOs)
      maxNPNOs = nPNOs;
    int index = (nPNOs - 1) / 5;
    if (index > 10)
      index = 10;
    pnoSpread[index] += 1;
    sumOfPNOs += pair->k_ij.cols();
    if (pair->pnoNormError > 0.999)
      pnoTraceRecovery[0] += 1;
    if (pair->pnoNormError <= 0.999 && pair->pnoNormError > 0.99)
      pnoTraceRecovery[1] += 1;
    if (pair->pnoNormError <= 0.99 && pair->pnoNormError > 0.9)
      pnoTraceRecovery[2] += 1;
    if (pair->pnoNormError <= 0.9)
      pnoTraceRecovery[3] += 1;
  }
  OutputControl::nOut << std::fixed;
  Eigen::VectorXd energies = calculateEnergyCorrection();
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  OutputControl::nOut << "  Number of strong pairs (CCSD)           " << orbitalPairs.size() << std::endl;
  OutputControl::nOut << "  Number of weak pairs (SC-MP2)           " << nStrongPlusWeakPairs - orbitalPairs.size()
                      << std::endl;
  OutputControl::nOut << "  Number of very distant pairs (dipole)   " << nVeryWeakPairs << std::endl;
  OutputControl::nOut << "  Total number of PNOs for strong pairs:  " << sumOfPNOs << std::endl;
  OutputControl::nOut << "  Average number of PNOs per strong pair: " << sumOfPNOs / orbitalPairs.size() << std::endl;
  OutputControl::nOut << "  Maximum number of PNOs                  " << maxNPNOs << std::endl;
  OutputControl::nOut << "  Semi-Canonical MP2 energy               " << std::setw(10) << energies.sum() << " E_h"
                      << std::endl;
  OutputControl::nOut << "    SC-MP2 (close pairs)                  " << std::setw(10) << energies[0] << " E_h" << std::endl;
  OutputControl::nOut << "    SC-MP2 (distant pairs)                " << std::setw(10) << energies[2] << " E_h" << std::endl;
  OutputControl::nOut << "    dipole approx. (very distant pairs)   " << std::setw(10) << energies[3] << " E_h" << std::endl;
  OutputControl::nOut << "    PNO truncation                        " << std::setw(10) << energies[4] << " E_h" << std::endl;
  OutputControl::nOut << "-----------------------------------------------------" << std::endl;
  OutputControl::nOut << std::scientific;
  OutputControl::vOut << "Pairs with 1-5   PNOs " << pnoSpread[0] << std::endl;
  OutputControl::vOut << "Pairs with 6-10  PNOs " << pnoSpread[1] << std::endl;
  OutputControl::vOut << "Pairs with 11-15 PNOs " << pnoSpread[2] << std::endl;
  OutputControl::vOut << "Pairs with 16-20 PNOs " << pnoSpread[3] << std::endl;
  OutputControl::vOut << "Pairs with 21-25 PNOs " << pnoSpread[4] << std::endl;
  OutputControl::vOut << "Pairs with 26-30 PNOs " << pnoSpread[5] << std::endl;
  OutputControl::vOut << "Pairs with 31-35 PNOs " << pnoSpread[6] << std::endl;
  OutputControl::vOut << "Pairs with 36-40 PNOs " << pnoSpread[7] << std::endl;
  OutputControl::vOut << "Pairs with 41-45 PNOs " << pnoSpread[8] << std::endl;
  OutputControl::vOut << "Pairs with 46-50 PNOs " << pnoSpread[9] << std::endl;
  OutputControl::vOut << "Pairs with  > 50 PNOs " << pnoSpread[10] << std::endl;
  OutputControl::vOut << "----------------------------------------------" << std::endl;
  OutputControl::vOut << "D_ij trace recovery" << std::endl;
  OutputControl::vOut << "Pairs with tr > 0.999          " << pnoTraceRecovery[0] << std::endl;
  OutputControl::vOut << "Pairs with tr in ]0.99,0.999]  " << pnoTraceRecovery[1] << std::endl;
  OutputControl::vOut << "Pairs with tr in ]0.9,0.99]    " << pnoTraceRecovery[2] << std::endl;
  OutputControl::vOut << "Pairs with tr <= 0.9           " << pnoTraceRecovery[3] << std::endl;
  OutputControl::vOut << "----------------------------------------------" << std::endl;
}

std::string DLPNO_CCSD::getTimeString() {
  timespec now;
  clock_gettime(CLOCK_REALTIME, &now);
  double sec = (double)(now.tv_sec - _time.tv_sec) + (now.tv_nsec - _time.tv_nsec) * 0.000000001;
  int ms = (int)(sec * 1000) - 1000 * (int)sec;
  std::string timeString =
      std::to_string((int)(sec / 60)) + ":" + std::to_string((int)(sec) % 60) + "." + std::to_string(ms) + " min:sec.ms";
  clock_gettime(CLOCK_REALTIME, &_time);
  return timeString;
}

inline void DLPNO_CCSD::optimizeAmplitudes() {
  std::vector<std::shared_ptr<OrbitalPair>> orbitalPairs =
      _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  std::vector<std::shared_ptr<SingleSubstitution>> singles = _localCorrelationController->getSingles();

  std::cout << std::endl;
  takeTime("Local CCSD Amplitude Optimization");
  Timings::takeTime("Local Cor. -    Amplitude Opt.");
  printSmallCaption("Local CCSD Amplitude Optimization");
  std::printf("%6s %12s %12s %12s %s\n", "Cycle", "abs. max. Res.", "Corr. Energy", "Delta E_corr", "Time");
  double largestResidual = 10;
  double correlationEnergy = 0.0;
  double deltaE_corr = 1e+9;
  unsigned int cycle = 0;

  auto& libint = Libint::getInstance();
  if (!_linearScalingSigmaVector) {
    OutputControl::nOut << "    Using integral direct sigma vector construction" << std::endl;
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  }
  else {
    OutputControl::nOut << "    Using linear scaling sigma vector construction" << std::endl;
  }

  const auto& settings = _localCorrelationController->getSettings();
  OrbitalPairDIISWrapper diis(settings.diisMaxStore);
  double damping = settings.dampingFactor;
  clock_gettime(CLOCK_REALTIME, &_time);
  unsigned int nThreads = 1;
  // Read integrals from file if necessary.
  std::vector<std::shared_ptr<OrbitalPairSet>> closePairSets = _localCorrelationController->getCloseOrbitalPairSets();
  std::shared_ptr<HDF5::H5File> file = tryLoadingIntegrals();
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
#endif
  while (largestResidual > _maxResidual || std::fabs(deltaE_corr) > _maxResidual) {
    largestResidual = 0.0;
    // Residual calculations
    std::vector<double> largestResidualVector = {0.0};
#ifdef _OPENMP
    largestResidualVector = std::vector<double>(nThreads, 0.0);
    Eigen::setNbThreads(1);
#endif
    // Pair/single dressing
    if (!_linearScalingSigmaVector)
      updateSigmaVector();
    dressSingles();
    dressPairs();
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iSingle = 0; iSingle < singles.size(); ++iSingle) {
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif
      std::shared_ptr<SingleSubstitution> single = singles[iSingle];
      single->residual = calculateSinglesResidual(single);
      double maxCoeff = single->residual.array().abs().maxCoeff();
      if (maxCoeff > largestResidualVector[threadId])
        largestResidualVector[threadId] = maxCoeff;
    } // for singles
    for (int i = closePairSets.size() - 1; i >= 0; --i) {
      OrbitalPairSet& orbitalPairSet = *closePairSets[i];
      if (!orbitalPairSet.integralsReady()) {
        orbitalPairSet.fromHDF5(*file);
      }
#pragma omp parallel for schedule(dynamic)
      for (unsigned int iPair = 0; iPair < orbitalPairSet.size(); ++iPair) {
#ifdef _OPENMP
        const unsigned int threadId = omp_get_thread_num();
#else
        const unsigned int threadId = 0;
#endif
        std::shared_ptr<OrbitalPair> pair = orbitalPairSet[iPair];
        if (pair->type != OrbitalPairTypes::CLOSE)
          continue;
        pair->residual = calculateDoublesResidual(pair);
        double maxCoeff = pair->residual.array().abs().maxCoeff();
        if (maxCoeff > largestResidualVector[threadId])
          largestResidualVector[threadId] = maxCoeff;
      } // for pair
      if (file)
        orbitalPairSet.removeInteralsFromMemory();
    } // for orbitalPairSet
    for (const auto& resEntry : largestResidualVector)
      if (resEntry > largestResidual)
        largestResidual = resEntry;
    // Amplitude update.
    if (damping > settings.finalDamping) {
      damping = (damping - settings.dampingChange < settings.finalDamping) ? settings.finalDamping
                                                                           : damping - settings.dampingChange;
    }
    double update = 1.0 - damping;
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iSingle = 0; iSingle < singles.size(); ++iSingle) {
      std::shared_ptr<SingleSubstitution> single = singles[iSingle];
      single->t_i.array() -= update * single->residual.array() / single->epsMinusF.array();
    }
#pragma omp parallel for schedule(dynamic)
    for (unsigned int iPair = 0; iPair < orbitalPairs.size(); ++iPair) {
      std::shared_ptr<OrbitalPair> pair = orbitalPairs[iPair];
      pair->t_ij.array() -= update * pair->residual.array() / pair->uncoupledTerm.array();
    }

    ++cycle;
    Eigen::VectorXd energies = calculateEnergyCorrection();
    double newEnergy = energies.sum();
    if (settings.diisStartResidual > largestResidual && cycle > 1)
      diis.optimize(orbitalPairs, singles);
    energies = calculateEnergyCorrection();
    newEnergy = energies.sum();
    deltaE_corr = energies.sum() - correlationEnergy;
    correlationEnergy = newEnergy;
    std::printf("%6d %12f %12f %12f %s\n", cycle, largestResidual, newEnergy, deltaE_corr, getTimeString().c_str());
    std::cout.flush();
    if (cycle > _maxCycles - 1) {
      throw SerenityError((std::string) "Canceling amplitude optimization after " + cycle + " cycles. NOT CONVERGED!!!");
    } // if cycle > _maxCycles-1
    if (largestResidual > 1e+3) {
      throw SerenityError("The local CCSD amplitude optimization encountered a very large residual!\n\
                           The amplitude optimization will fail and will be stopped here.\n\
                           Typical error sources are:\n\
                             1. Wrong Fock matrix (use HF orbitals, check embedding settings)\n\
                             2. Inconsistent use of PAO normalization with respect to canonical orthogonalization\n\
                             3. To loose integral threshold.");
    }
  } // while largestResidual > _maxResidual
  OutputControl::mOut << " Amplitude Optimization Converged!" << std::endl;
  runDiagnostics();
#ifdef _OPENMP
  Eigen::setNbThreads(nThreads);
#endif
  if (file)
    file->close();
  if (!_linearScalingSigmaVector) {
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
    libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  }
  timeTaken(1, "Local CCSD Amplitude Optimization");
  Timings::timeTaken("Local Cor. -    Amplitude Opt.");
}

void DLPNO_CCSD::runDiagnostics() {
  auto singles = _localCorrelationController->getSingles();
  unsigned int nThreads = 1;
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
#endif
  std::vector<double> norms(nThreads, 0.0);
  std::vector<double> largestT1(nThreads, 0.0);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iSingle = 0; iSingle < singles.size(); ++iSingle) {
#ifdef _OPENMP
    const unsigned int threadId = omp_get_thread_num();
#else
    const unsigned int threadId = 0;
#endif
    const std::shared_ptr<SingleSubstitution>& single = singles[iSingle];
    norms[threadId] += single->t_i.squaredNorm();
    double maxT1 = single->t_i.array().abs().maxCoeff();
    largestT1[threadId] = std::max(largestT1[threadId], std::fabs(maxT1));
  }
  double singlesNorm = 0.0;
  double maxT1 = 0.0;
  for (unsigned int threadID = 0; threadID < norms.size(); ++threadID) {
    singlesNorm += norms[threadID];
    maxT1 = std::max(largestT1[threadID], maxT1);
  }
  singlesNorm = sqrt(singlesNorm);
  OutputControl::nOut << std::fixed;
  OutputControl::nOut << " Singles Norm sqrt(<S|S>):  " << singlesNorm << std::endl;
  unsigned int nElectrons = _localCorrelationController->getActiveSystemController()->getNElectrons<RESTRICTED>();
  double t_1Diagnostic = singlesNorm / sqrt((double)nElectrons);
  OutputControl::nOut << " T_1 diagnostic             " << t_1Diagnostic << std::endl;
  OutputControl::nOut << " Abs. max T_1 amplitude     " << maxT1 << std::endl;
  if (t_1Diagnostic > 0.02)
    WarningTracker::printWarning("WARNING: The T_1 diagnostic exceeds the threshold (0.02) commonly used for detecting\n\
         multireference character. Note that for transition metal species the threshold is\n\
         less strict:\n\
           3d species: 0.05\n\
           4d species: 0.045\n\
         according to JCTC 11, 5865-5872 (11).",
                                 true);
  OutputControl::nOut << std::scientific;
  OutputControl::nOut.flush();
}

inline void DLPNO_CCSD::calculateCCSDIntegrals() {
  auto closePairSets = _localCorrelationController->getCloseOrbitalPairSets();
  const auto activeSystem = _localCorrelationController->getActiveSystemController();
  Ao2MoExchangeIntegralTransformer::transformAllIntegrals(
      activeSystem->getBasisController(Options::BASIS_PURPOSES::AUX_CORREL),
      _localCorrelationController->getMO3CenterIntegralController(true), closePairSets,
      _localCorrelationController->getSettings().dumpIntegrals || closePairSets.size() > 1,
      _localCorrelationController->getPairIntegralFileName(), _linearScalingSigmaVector,
      _localCorrelationController->getSettings().lowMemory,
      _localCorrelationController->getSettings().ignoreMemoryConstraints);
  if (!dlpnoCCSDSettings.keepMO3CenterIntegrals)
    _localCorrelationController->removeMO3CenterIntegralController();
}

Eigen::VectorXd DLPNO_CCSD::calculateEnergyCorrection() {
  auto closePairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  auto distantPairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::DISTANT);
  auto veryDistantPairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::VERY_DISTANT);
  auto singles = _localCorrelationController->getSingles();

  // Singles contribution
  double ccSinglesCorrection = 0.0;
  for (const auto& single : singles)
    ccSinglesCorrection += 2.0 * single->f_ai.transpose() * single->t_i;
  // Doubles contribution
  double ccDoublesCorrection = 0.0;
  double pnoTruncation = 0.0;

  Eigen::setNbThreads(1);
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iPair = 0; iPair < closePairs.size(); ++iPair) {
    std::shared_ptr<OrbitalPair> pair = closePairs[iPair];
    calculate_tau_ij(pair);
    pair->dlpnoCCSDPairEnergy = ((pair->i == pair->j) ? 1.0 : 2.0) *
                                ((2.0 * pair->k_ij - pair->k_ij.transpose()).array() * pair->tau_ij.array()).sum();
  } // for pair
  Eigen::setNbThreads(0);
  for (const auto& pair : closePairs) {
    ccDoublesCorrection += pair->dlpnoCCSDPairEnergy;
    pnoTruncation += pair->deltaPNO;
  }
  // distant pairs (Semi-canonical pair energies)
  double scMP2Correction = 0.0;
  for (const auto& pair : distantPairs)
    scMP2Correction += pair->scMP2PairEnergy;
  // very distant pairs (dipole approximation)
  double dipoleCorrection = 0.0;
  for (const auto& pair : veryDistantPairs) {
    if (pair->scMP2PairEnergy != 0.0) {
      scMP2Correction += pair->scMP2PairEnergy;
    }
    else {
      dipoleCorrection += pair->dipolePairEnergy;
    }
  }
  Eigen::VectorXd energies(5);
  energies << ccDoublesCorrection, ccSinglesCorrection, scMP2Correction, dipoleCorrection, pnoTruncation;
  return energies;
}

void DLPNO_CCSD::updateSigmaVector() {
  auto activeSystem = _localCorrelationController->getActiveSystemController();
  _g_ao_ao = std::make_shared<MatrixInBasis<RESTRICTED>>(LocalAOSigmaVectorWrapper::getSigmaVector_AO_AO(
      activeSystem, _localCorrelationController->getSingles(), _localCorrelationController->getPAOController(), testRun));
  unsigned int nOcc = activeSystem->getNOccupiedOrbitals<RESTRICTED>();
  const Eigen::MatrixXd occCoefficients =
      activeSystem->getActiveOrbitalController<RESTRICTED>()->getCoefficients().leftCols(nOcc);
  _g_ao_occ = _g_ao_ao->transpose() * occCoefficients;
  _g_occ_occ = _g_ao_occ.transpose() * occCoefficients;
}

std::shared_ptr<HDF5::H5File> DLPNO_CCSD::tryLoadingIntegrals() {
  std::vector<std::shared_ptr<OrbitalPairSet>> closePairSets = _localCorrelationController->getCloseOrbitalPairSets();
  const double totalMemoryAvailable = _localCorrelationController->getSettings().maximumMemoryRatio *
                                      MemoryManager::getInstance()->getAvailableSystemMemory();
  double memoryUsed = 0.0;
  bool allLoaded = true;
  std::shared_ptr<HDF5::H5File> file = nullptr;
  for (int i = closePairSets.size() - 1; i >= 0; --i) {
    auto orbitalPairSet = closePairSets[i];
    if (!orbitalPairSet->integralsReady()) {
      if (!file) {
        std::string fileName = _localCorrelationController->getPairIntegralFileName();
        HDF5::Filepath name(fileName);
        file = std::make_shared<HDF5::H5File>(name.c_str(), H5F_ACC_RDONLY);
      } // if !file
      if (memoryUsed + orbitalPairSet->memoryDemand(_linearScalingSigmaVector) < totalMemoryAvailable) {
        orbitalPairSet->fromHDF5(*file);
        memoryUsed += orbitalPairSet->memoryDemand(_linearScalingSigmaVector);
      }
      else {
        OutputControl::dOut << "Unable to hold all integrals in memory! Integrals will be read from disk!" << std::endl;
        allLoaded = false;
        break;
      }
    } // if !orbitalPairSet->integralsReady()
  }
  if (allLoaded && file) {
    file->close();
    return nullptr;
  }
  return file;
}

void DLPNO_CCSD::switchIntegrals() {
  auto activeSystem = _localCorrelationController->getActiveSystemController();
  const Eigen::MatrixXd s = activeSystem->getOneElectronIntegralController()->getOverlapIntegrals();
  const auto paoController = _localCorrelationController->getPAOController();
  auto closePairs = _localCorrelationController->getOrbitalPairs(OrbitalPairTypes::CLOSE);
  const unsigned int nOcc = activeSystem->getNOccupiedOrbitals<RESTRICTED>();
  const unsigned int nBasisFunc = activeSystem->getBasisController()->getNBasisFunctions();
  const unsigned int nCanVir = nBasisFunc - nOcc;
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeffs =
      activeSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getMolecularOrbitals()->getCoefficients();
  const Eigen::VectorXd eigenvalues =
      activeSystem->getElectronicStructure<Options::SCF_MODES::RESTRICTED>()->getMolecularOrbitals()->getEigenvalues();
  // Calculate fully transformed integrals
  TwoElecFourCenterIntLooper looper(LIBINT_OPERATOR::coulomb, 0, activeSystem->getBasisController(), 1E-10);
  RegularRankFourTensor<double> eris(nBasisFunc, 0.0);
  Ao2MoTransformer aoToMo(activeSystem->getBasisController());
  auto const storeERIS = [&eris](const unsigned int& a, const unsigned int& b, const unsigned int& i,
                                 const unsigned int& j, const Eigen::VectorXd& integral, const unsigned int threadId) {
    (void)threadId; // no warning, please
    eris(b, a, i, j) = integral(0);
    eris(b, a, j, i) = integral(0);
    eris(a, b, j, i) = integral(0);
    eris(a, b, i, j) = integral(0);
    eris(i, j, b, a) = integral(0);
    eris(i, j, a, b) = integral(0);
    eris(j, i, b, a) = integral(0);
    eris(j, i, a, b) = integral(0);
  };
  looper.loop(storeERIS);

  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nCanVir, nCanVir);
  _localCorrelationController->getDomainOverlapMatrixController()->setIdentity(std::make_shared<Eigen::MatrixXd>(I));
  aoToMo.transformTwoElectronIntegrals(eris, eris, coeffs, nBasisFunc);
  const auto result = eris;
  for (auto& pair : closePairs) {
    const unsigned int i = pair->i;
    const unsigned int j = pair->j;
    pair->singles_i->t_i = Eigen::VectorXd::Zero(nCanVir);
    pair->singles_j->t_i = Eigen::VectorXd::Zero(nCanVir);
    pair->singles_i->f_ai = Eigen::VectorXd::Zero(nCanVir);
    pair->singles_j->f_ai = Eigen::VectorXd::Zero(nCanVir);
    pair->singles_i->epsMinusF.array() = eigenvalues.segment(nOcc, nCanVir).array() - eigenvalues(i);
    pair->singles_j->epsMinusF.array() = eigenvalues.segment(nOcc, nCanVir).array() - eigenvalues(j);
    pair->f_ab = eigenvalues.segment(nOcc, nCanVir);
    pair->singles_i->f_ab = pair->f_ab;
    pair->singles_j->f_ab = pair->f_ab;
    // (ac|bd)
    unsigned int nPNOs = nCanVir;
    pair->ac_bd = std::make_unique<Matrix<Eigen::MatrixXd>>(nPNOs, nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));
    pair->ia_bc = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));
    pair->ja_bc = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));
    pair->jc_ab = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));
    pair->ic_ab = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));
    for (unsigned int k = 0; k < nOcc; ++k) {
      pair->coupledPairs[k]->ka_bc = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));
      pair->coupledPairs[k]->ab_kcX2_M_ak_bc = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nPNOs, nPNOs));
      pair->coupledPairs[k]->ia_kc = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
      pair->coupledPairs[k]->ja_kc = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
      pair->coupledPairs[k]->ik_ca = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
      pair->coupledPairs[k]->jk_ca = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
      pair->coupledPairs[k]->ij_ak = Eigen::VectorXd::Zero(nPNOs);
      pair->coupledPairs[k]->ja_ik = Eigen::VectorXd::Zero(nPNOs);
      pair->coupledPairs[k]->ia_jk = Eigen::VectorXd::Zero(nPNOs);
    }
    pair->k_ij = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
    pair->ij_ab = Eigen::MatrixXd::Zero(nPNOs, nPNOs);
    pair->ki_la = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nOcc, nOcc));
    pair->kj_la = std::vector<Eigen::MatrixXd>(nPNOs, Eigen::MatrixXd::Zero(nOcc, nOcc));
    // Check symmetry of the integrals and compare to AO integrals if possible.
    for (unsigned int a = 0; a < nPNOs; ++a) {
      for (unsigned int b = 0; b < nPNOs; ++b) {
        for (unsigned int c = 0; c < nPNOs; ++c) {
          for (unsigned int d = 0; d < nPNOs; ++d) {
            (*pair->ac_bd)(a, b)(c, d) = result(a + nOcc, c + nOcc, b + nOcc, d + nOcc);
          } // for d
          pair->ia_bc[c](a, b) = result(i, a + nOcc, b + nOcc, c + nOcc);
          pair->ja_bc[c](a, b) = result(j, a + nOcc, b + nOcc, c + nOcc);
          pair->jc_ab[c](a, b) = result(j, c + nOcc, a + nOcc, b + nOcc);
          pair->ic_ab[c](a, b) = result(i, c + nOcc, a + nOcc, b + nOcc);
          for (unsigned int k = 0; k < nOcc; ++k) {
            pair->coupledPairs[k]->ka_bc[c](a, b) = result(k, a + nOcc, b + nOcc, c + nOcc);
            pair->coupledPairs[k]->ab_kcX2_M_ak_bc[c](a, b) =
                2.0 * result(a + nOcc, b + nOcc, k, c + nOcc) - result(a + nOcc, k, b + nOcc, c + nOcc);
          }
        } // for c
        pair->k_ij(a, b) = result(i, a + nOcc, j, b + nOcc);
        pair->ij_ab(a, b) = result(i, j, a + nOcc, b + nOcc);
        for (unsigned int k = 0; k < nOcc; ++k) {
          pair->coupledPairs[k]->ia_kc(a, b) = result(i, a + nOcc, k, b + nOcc);
          pair->coupledPairs[k]->ja_kc(a, b) = result(j, a + nOcc, k, b + nOcc);
          pair->coupledPairs[k]->ik_ca(a, b) = result(i, k, a + nOcc, b + nOcc);
          pair->coupledPairs[k]->jk_ca(a, b) = result(j, k, a + nOcc, b + nOcc);
        }
      } // for b
      for (unsigned int k = 0; k < nOcc; ++k) {
        pair->coupledPairs[k]->ij_ak(a) = result(i, j, a + nOcc, k);
        pair->coupledPairs[k]->ja_ik(a) = result(j, a + nOcc, i, k);
        pair->coupledPairs[k]->ia_jk(a) = result(i, a + nOcc, j, k);
      }
    } // for a
    pair->iaS_jbSX2_M_ij_aSbS = 2.0 * pair->k_ij - pair->ij_ab;
    pair->ik_jl = Eigen::MatrixXd::Zero(nOcc, nOcc);
    for (auto& klPairSet : pair->klPairSets) {
      klPairSet->ki_la = Eigen::VectorXd(nPNOs);
      klPairSet->kj_la = Eigen::VectorXd(nPNOs);
      klPairSet->li_ka = Eigen::VectorXd(nPNOs);
      klPairSet->lj_ka = Eigen::VectorXd(nPNOs);
      unsigned int k = klPairSet->getKLPair()->i;
      unsigned int l = klPairSet->getKLPair()->j;
      klPairSet->ik_jl = result(i, k, j, l);
      klPairSet->il_jk = result(i, l, j, k);
      pair->ik_jl(k, l) = result(i, k, j, l);
      pair->ik_jl(l, k) = result(i, l, j, k);
      for (unsigned int a = 0; a < nPNOs; ++a) {
        klPairSet->ki_la(a) = result(k, i, l, a + nOcc);
        klPairSet->kj_la(a) = result(k, j, l, a + nOcc);
        klPairSet->li_ka(a) = result(l, i, k, a + nOcc);
        klPairSet->lj_ka(a) = result(l, j, k, a + nOcc);
        pair->ki_la[a](k, l) = result(k, i, l, a + nOcc);
        pair->ki_la[a](l, k) = result(l, i, k, a + nOcc);
        pair->kj_la[a](k, l) = result(k, j, l, a + nOcc);
        pair->kj_la[a](l, k) = result(l, j, k, a + nOcc);
      }
    }
    Eigen::MatrixXd uncoupledTerm = Eigen::MatrixXd::Zero(nCanVir, nCanVir);
    uncoupledTerm.rowwise() += eigenvalues.segment(nOcc, nCanVir).transpose();
    uncoupledTerm.colwise() += eigenvalues.segment(nOcc, nCanVir);
    uncoupledTerm.array() -= eigenvalues(i) + eigenvalues(j);
    pair->uncoupledTerm = uncoupledTerm;
    pair->t_ij.array() = -pair->k_ij.array() / pair->uncoupledTerm.array();
  }
  std::cout << "Over writing done!" << std::endl;
}

} /* namespace Serenity */
