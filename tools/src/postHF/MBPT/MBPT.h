/**
 * @file MBPT.h
 *
 * @date Apr 14, 2020
 * @author Johannes Toelle
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
 */

#ifndef POSTHF_MBPT_MBPT_H_
#define POSTHF_MBPT_MBPT_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "tasks/GWTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>

namespace Serenity {

/* Forward declarations */
class SystemController;
class DIIS;

/**
 * @class  MBPT MBPT.h
 * @brief  Calculates MBPT correlation energy (dRPA) and helper functions used for GW/dRPA
 *
 * Coded in reference to:\n
 * New J. of Phys. 14 053020 2012\n
 * JCTC 9 4856–4869 2018\n
 * https://publikationen.bibliothek.kit.edu/1000095752
 */
template<Options::SCF_MODES SCFMode>
class MBPT {
 public:
  /**
   * @brief Constructor
   * @param systemController
   * @param settings The GW task settings
   * @param envSystemController
   * @param riInts The RI integrals
   * @param startOrb The first orbital index included in a GW calculation
   * @param endOrb The last orbital index included in a GW calculation
   */
  MBPT(std::shared_ptr<LRSCFController<SCFMode>> lrscf, GWTaskSettings settings,
       std::vector<std::shared_ptr<SystemController>> envSystemController,
       std::shared_ptr<RIIntegrals<SCFMode>> riInts = nullptr, int startOrb = 0, int endOrb = 0);
  /**
   * @brief Default destructor;
   */
  virtual ~MBPT() = default;
  /**
   * @brief Calculates the static environmental response in auxiliary basis
   * @return environmental response
   */
  Eigen::MatrixXd environmentRespose();
  /**
   * @brief Calculates virtual-occupied orbital energy difference
   * @param orbEigenValues The orbital eigenvalues
   * @param nOcc Number of occupied orbitals
   * @param nVirt Number of virtual orbitals
   * @return virtual-occupied orbital energy difference
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> calculateEia(SpinPolarizedData<SCFMode, Eigen::VectorXd>& orbEigenValues,
                                                           SpinPolarizedData<SCFMode, unsigned int>& nOcc,
                                                           SpinPolarizedData<SCFMode, unsigned int>& nVirt);
  /**
   * @brief Calculates the transformation matrix from the environmental screening
   * @param transformation The transformation matrix
   * @param environmentResponse The environment response
   * @return sparse projection matrix
   */
  Eigen::SparseMatrix<double> calculateTransformation(Eigen::MatrixXd& transformation, Eigen::MatrixXd& environmentResponse);
  /**
   * @brief Calculates the screened Coulomb interaction for complex frequencies
   * @return screened Coulomb interaction
   */
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> calculateWnmComplex();
  /**
   * @brief Calculates the screened Coulomb interaction for complex frequencies using Laplace transform
   * @return screened Coulomb interaction
   */
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> calculateWnmComplexLT();
  /**
   * @brief Calculates the fermi level
   * @return fermi level
   */
  SpinPolarizedData<SCFMode, double> calculateFermiLevel();
  /**
   * @brief Calculates Pi(w) in auxiliary basis
   * @param frequency The complex fequency
   * @param jia The JiaQ integrals
   * @param e_ia The occupied virtual orbital energy difference
   * @return Pi(w) in auxiliary basis
   */
  Eigen::MatrixXd calculatePiOmega(std::complex<double>& frequency, Eigen::MatrixXd& jia, Eigen::VectorXd& e_ia);
  /**
   * @brief Calculates dRPA correlation energy
   * @return dRPA correlation energy
   */
  double calculateRPACorrelationEnergy();
  /**
   * @brief Checks for convergence in G0W0 quasi-particle iterations/ evGW
   * @param new_qp The new quasi-particle energy
   * @param old_qp The old quasi-particle energy
   * @param diis The DIIS convergence accelerator
   * @param evGW Whether the calcualation is a evGW calculation
   * @param iteration The current iteration
   * @return whether the calculation is converged
   */
  bool convergenCheck(SpinPolarizedData<SCFMode, Eigen::VectorXd>& new_qp, SpinPolarizedData<SCFMode, Eigen::VectorXd>& old_qp,
                      std::shared_ptr<DIIS> diis, bool evGW, unsigned int iteration);
  /**
   * @brief Prints the quasi-particle iteration information
   * @param new_qp The new quasi-particle energy
   * @param old_qp The old quasi-particle energy
   */
  void printQPInfo(SpinPolarizedData<SCFMode, Eigen::VectorXd>& new_qp, SpinPolarizedData<SCFMode, Eigen::VectorXd>& old_qp);
  /**
   * @brief Shifts the occupied and virtual orbitals not included in the GW calculation \n
   *        by the change of the highest occupied/lowest virtual orbital relative to the underlying mean-field
   * calculation
   * @param new_qp The new quasi-particle energy
   * @param old_qp The old quasi-particle energy
   */
  void shiftByGap(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEner, SpinPolarizedData<SCFMode, Eigen::VectorXd>& startOrbEner);
  /// @brief Sets virtual occupied orbital energy difference
  void setEia(SpinPolarizedData<SCFMode, Eigen::VectorXd>& eia) {
    _eia = eia;
  };
  /// @brief Sets the quasi-particle energies
  void setQPEner(SpinPolarizedData<SCFMode, Eigen::VectorXd>& qpEner) {
    _orbEig = qpEner;
  };

 protected:
  /// @brief systemcontroller
  std::shared_ptr<LRSCFController<SCFMode>> _lrscfController;
  /// @brief GWTask settings
  GWTaskSettings _settings;
  /// @brief environmental systems
  std::vector<std::shared_ptr<SystemController>> _env;
  /// @brief RI integrals
  std::shared_ptr<RIIntegrals<SCFMode>> _riInts;
  /// @brief orbital eigenvalues
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _orbEig;
  /// @brief start orbital index
  int _startOrb;
  /// @brief end orbital index
  int _endOrb;
  /// @brief the ri integrals needed for the calculation
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jia, _Jpq, _Jia_transformed, _Jpq_transformed;
  /// @brief the occupied orbitals
  SpinPolarizedData<SCFMode, unsigned int> _nOcc;
  /// @brief the virtual orbitals
  SpinPolarizedData<SCFMode, unsigned int> _nVirt;
  /// @brief the virtual-occupied orbital energy difference
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _eia;
  /// @brief the numerical integration nodes
  Eigen::VectorXd _nodes;
  /// @brief the numerical integration weights
  Eigen::VectorXd _weights;
  /// @brief the number of auxiliary functions
  unsigned int _nAux;
  /// @brief the sparse projection matrix
  Eigen::SparseMatrix<double> _proj;

  Eigen::MatrixXd _envResponse;

  Eigen::VectorXd _eValues;
};

} /* namespace Serenity */

#endif /* POSTHF_MBPT_MBPT_H_ */