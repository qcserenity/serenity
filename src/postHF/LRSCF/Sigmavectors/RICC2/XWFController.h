/**
 * @file XWFController.h
 *
 * @date Apr 20, 2020
 * @author Niklas Niemeyer
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

#ifndef LRSCF_XWFCONTROLLER
#define LRSCF_XWFCONTROLLER

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/LRSCFController.h"

/* Include Std and External Headers */

namespace Serenity {
class TwoElecThreeCenterCalculator;
/**
 * @class XWFController XWFController.h
 * @brief Base class for ADC(2), CIS(D) and CC2 sigma vectors.\n
 *
 * The whole CC2 framework loosely follows the algorithms and working
 * equations devised by Hättig et al., at least for CC2.\n
 * The main references along with naming conventions for intermediates etc.
 * can be found below.\n
 *
 * Please note that there are some typos in there which should be
 * corrected in the cheat sheet that I have added to the Serenity manual.
 *   -- Niklas.
 *
 * [1] C. Hättig, F. Weigend, J. Chem. Phys. 2000, 113, 5154-5161.
 * [2] C. Hättig, A. Köhn, J. Chem. Phys. 2002, 117, 6939-6951.
 *
 */
template<Options::SCF_MODES SCFMode>
class XWFController {
  friend class LRSCFController<SCFMode>;

 public:
  /**
   * @brief Constructor.
   * @param lrscf LRSCFController holding all the necessary information.
   */
  XWFController<SCFMode>(std::shared_ptr<LRSCFController<SCFMode>> lrscf);

  /**
   * @brief Destructor.
   */
  virtual ~XWFController() = default;

  /**
   * @brief Initiliazes this controller, calculates the E intermediate and the ground state.
   */
  void initialize();

  /**
   * @brief Return the right Jacobian transformation of the incoming guess vector.\n
   *         Must be overriden by derived class.
   * @param guessVector The guess vector.
   * @param eigenvalue The current eigenvalue (non-linearity of the eigenvalue problem).
   * @return The sigma vector.
   */
  virtual Eigen::VectorXd getRightXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) = 0;

  /**
   * @brief Return the left Jacobian transformation of the incoming guess vector.\n
   *        Must be overriden by derived class.
   * @param guessVector The guess vector.
   * @param eigenvalue The current eigenvalue (non-linearity of the eigenvalue problem).
   * @return The sigma vector.
   */
  virtual Eigen::VectorXd getLeftXWFSigma(Eigen::Ref<Eigen::VectorXd> guessVector, double eigenvalue) = 0;

  /**
   * @brief Normalizes eigenvectors, either <R|R> = 1 (one set of eigenvectors) or the <L|R> = 1 criterion (two).
   * @param eigenvectors The eigenvectors to be normalized.
   * @param eigenvectors The associated eigenvalues.
   */
  void normalizeEigenvectors(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues);

  /**
   * @brief Calculates the ground-state Lagrange multiplier.
   */
  virtual void calcGroundStateLagrangeMultiplier() = 0;

  /**
   * @brief Calculates the transition-moment Lagrange multiplier.
   * @param eigenvectors The eigenvectors.
   * @param eigenvectors The associated eigenvalues.
   */
  virtual void calcTransitionMomentLagrangeMultiplier(std::vector<Eigen::MatrixXd>& eigenvectors,
                                                      Eigen::VectorXd eigenvalues) = 0;

  /**
   * @brief Calculates the one-particle density matrices.
   * @param eigenvectors The eigenvectors.
   * @param eigenvectors The associated eigenvalues.
   * @param densityMatrices The densityMatrices to be calcualted.
   */
  virtual void calcDensityMatrices(std::vector<Eigen::MatrixXd>& eigenvectors, Eigen::VectorXd eigenvalues,
                                   std::vector<Eigen::MatrixXd>& densityMatrices) = 0;

 protected:
  ///@brief The LRSCFController.
  std::weak_ptr<LRSCFController<SCFMode>> _lrscf;

  ///@brief HF energy.
  double _hfEnergy;

  ///@brief The MP2 energy.
  double _mp2Energy;

  ///@brief The correlation energy.
  double _corrEnergy;

  ///@brief Orbital energies.
  SpinPolarizedData<SCFMode, Eigen::VectorXd> _e;

  ///@brief SCF coefficient matrix.
  CoefficientMatrix<SCFMode> _C;

  ///@brief [V^{-1/2}] [* naf-matrix].
  Eigen::MatrixXd _riTrafo;

  ///@brief Left singles coefficient matrix.
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _P;

  ///@brief Right singles coefficient matrix.
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _H;

  ///@brief Number occupied orbitals.
  SpinPolarizedData<SCFMode, unsigned> _no;

  ///@brief Number virtual orbitals.
  SpinPolarizedData<SCFMode, unsigned> _nv;

  ///@brief Simple helper ints.
  unsigned _nv_a, _nv_b, _no_a, _no_b, _alpha, _beta;

  ///@brief The dimension of the single excitation space.
  unsigned _nDim;

  ///@brief Number of basis functions.
  unsigned _nb;

  ///@brief Number of auxiliary functions (might be < nxb if naf is used).
  unsigned _nx;

  ///@brief Number of actual auxiliary functions.
  unsigned _nxb;

  ///@brief Same-spin and opposite-spin scaling (and, for convenience, the sum of both).
  double _sss, _oss, _soss;

  ///@brief Laplace Transformation roots.
  Eigen::VectorXd _roots;

  ///@brief Laplace Transformation weights.
  Eigen::VectorXd _weights;

  ///@brief A transformation vector.
  Eigen::VectorXd _trafoVector;

  ///@brief Ground-state singles amplitudes (zero unless CC2).
  Eigen::VectorXd _singles;

  ///@brief Residual of the vectorfunction (zero unless CC2).
  Eigen::VectorXd _residual;

  ///@brief Iteration counter for ground state calculation.
  unsigned _iter;

  ///@brief Integral lists (ia|Q).
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jia;

  ///@brief Integral lists (ij|Q).
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Jij;

  ///@brief Singles-transformed integral lists (ia|Q).
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Bia;

  ///@brief Singles-transformed integral lists (ij|Q).
  std::shared_ptr<SpinPolarizedData<SCFMode, Eigen::MatrixXd>> _Bij;

  ///@brief Intermediates (only live inside the sigma vector object).
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _Xia, _Yia, _Zia;

  ///@brief Occ-occ part of the E intermediate.
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _Eij;

  ///@brief virt-virt part of the E intermediate.
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _Eab;

  ///@brief Used to store transformed Fock matrices.
  Eigen::VectorXd _Fai, _Fia;

  ///@brief Ground-state Lagrange multiplier.
  Eigen::VectorXd _gsLagrange;

  ///@brief Occ part of the amplitude denominator.
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _eOcc;

  ///@brief Virt part of the amplitude denominator.
  SpinPolarizedData<SCFMode, Eigen::MatrixXd> _eVirt;

  ///@brief Occ part of the amplitude denominator (mixed spin).
  Eigen::MatrixXd _eOccab;

  ///@brief Virt part of the amplitude denominator (mixed spin).
  Eigen::MatrixXd _eVirtab;

  ///@brief Looper for RI integrals.
  std::shared_ptr<TwoElecThreeCenterCalculator> _integrals;

  ///@brief Enum for underlying XWF method.
  Options::LR_METHOD _xwfModel;

  ///@brief The LRSCFTask settings.
  const LRSCFTaskSettings& _settings;

  ///@brief Prescreening threshold for integral-direct algorithms.
  double _prescreeningThreshold;

  /**
   * @brief Calculates the E-intermediate (must be overriden).
   */
  virtual void calculateE() = 0;

  /**
   * @brief Calculate the ground state (MP2 or CC2).
   */
  void calculateGroundstate();

  /**
   * @brief Calculates the CC2 residual function (must be overriden).
   */
  virtual void calculateResidual() = 0;

  /**
   * @brief Applies the singles similarity transformation to the (ia|Q) and (ij|Q) integrals.
   */
  void transformIntegrals();

  /**
   * @brief Calculates left- or right-transformed integrals for amplitudes.
   * @param Xia The tensor to hold the intermediate to be calculated.
   * @param trafoVector The singles part of the transformation.
   * @param fromLeft Perform transformation from the left or not.
   * @param iSpin Spin identifier.
   * @param singlesTransformed Use similarity transformed integrals (Bia vs Jia) or not.
   */
  void performTransformation(Eigen::MatrixXd& Xia, Eigen::Ref<Eigen::VectorXd> trafoVector, bool fromLeft,
                             int iSpin = 1, bool singlesTransformed = true);

  /**
   * @brief Return J.2 and G like contributions in an integral-direct fashion.
   * @param YiaPtr Yia intermediate to be contracted.
   * @param Jij Jij intermediate to be contracted.
   * @param guessVector The guess vector to be contracted.
   * @param fromLeft Perform transformation from the left or not.
   * @param iSpin Spin identifier.
   * @return The sigma vector-like contributions.
   */
  Eigen::VectorXd getJ2GContribution(double* YiaPtr, double* JijPtr, Eigen::Ref<Eigen::VectorXd> guessVector,
                                     bool fromLeft, int iSpin = 1);

  /**
   * @brief Performs a transformation with the sqrt of the inverse Coulomb metric.
   * @param Yia Yia intermediate to be transformed.
   * @param no Number of occupied orbitals (for block-wise transformation).
   * @param transposeRITrafo Q -> P transformation?
   */
  void performRITransformation(Eigen::MatrixXd& Yia, unsigned no, bool transposeRITrafo = false);

  /**
   * @brief Reorders Yia-type tensor in occ and virt index (can also be a single vector).
   * @param Yia The tensor to be reordered.
   * @param n The old major index in each column (to be swapped with m).
   * @param m The old minor index in each column (to be swapped with n).
   * @param thirdIndex How many columns there are.
   */
  void reorderTensor(Eigen::Ref<Eigen::MatrixXd> Yia, unsigned n, unsigned m, unsigned thirdIndex = 1);

  /**
   * doubles amplitudes/multiplier to be
   * calculated and contracted on-the-fly
   * these can return all amplitudes for
   * a given ij-pair or ab-pair.
   *
   * affixes:
   *   A: not antisymmetrized
   *   V: ab-wise
   */

  // Ground-state doubles amplitudes.
  Eigen::MatrixXd getAmplitudes(unsigned i, unsigned j, int iSpin = 1);
  Eigen::MatrixXd getAmplitudesA(unsigned i, unsigned j, int iSpin = 1);
  Eigen::MatrixXd getAmplitudesV(unsigned a, unsigned b, int iSpin = 1);
  Eigen::MatrixXd getAmplitudesAV(unsigned a, unsigned b, int iSpin = 1);

  // Right excited-state doubles amplitudes.
  Eigen::MatrixXd getRightAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin = 1);
  Eigen::MatrixXd getRightAmplitudesA(unsigned i, unsigned j, double eigenvalue, int iSpin = 1);
  Eigen::MatrixXd getRightAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin = 1);
  Eigen::MatrixXd getRightAmplitudesAV(unsigned a, unsigned b, double eigenvalue, int iSpin = 1);

  // Left excited-state doubles amplitudes.
  Eigen::MatrixXd getLeftAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin = 1);
  Eigen::MatrixXd getLeftAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin = 1);

  // Fock doubles multiplier.
  Eigen::MatrixXd getFockAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin = 1);

  // Ground-state Lagrangian doubles multiplier.
  Eigen::MatrixXd getGLagrangeAmplitudes(unsigned i, unsigned j, int iSpin = 1);
  Eigen::MatrixXd getGLagrangeAmplitudesV(unsigned a, unsigned b, int iSpin = 1);

  // Transition-moment Lagrangian doubles multiplier.
  Eigen::MatrixXd getELagrangeAmplitudes(unsigned i, unsigned j, double eigenvalue, int iSpin = 1);
  Eigen::MatrixXd getELagrangeAmplitudesV(unsigned a, unsigned b, double eigenvalue, int iSpin = 1);
}; // class XWFController
} // namespace Serenity

#endif /* LRSCF_XWFCONTROLLER */
