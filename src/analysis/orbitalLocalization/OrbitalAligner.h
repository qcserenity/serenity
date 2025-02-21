/**
 * @file OrbitalAligner.h
 *
 * @date Mar 14, 2019
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

#ifndef ANALYSIS_ORBITALLOCALIZATION_ORBITALALIGNER_H_
#define ANALYSIS_ORBITALLOCALIZATION_ORBITALALIGNER_H_

/* Include Serenity Internal Headers */
#include "analysis/orbitalLocalization/Localization.h" //Base class.
#include "data/matrices/SPMatrix.h"                    //spin-polarized matrices.
/* Include Std and External Headers */
#include <Eigen/Dense>  //Dense matrices.
#include <Eigen/Sparse> //Sparse matrices.
#include <memory>       //smrt_ptr

namespace Serenity {

/* Forward declarations */
template<Options::SCF_MODES SCFMode>
class OrbitalController;
class SystemController;

/**
 * @class OrbitalAligner OrbitalAligner.h
 * @brief This class implements orbital alignment strategies that are supposed to make a given set
 *        of orbitals as similar to a given template set as possible.
 *        The orbitals can be aligned using partial charges and the orbital kinetic energy.
 *        For more information see the associated pdf/tex document in the manual-directory.
 */
template<Options::SCF_MODES SCFMode>
class OrbitalAligner : public Localization<SCFMode> {
 public:
  /**
   * @brief Constructor. Performs sanity checks.
   * @param systemController  The system controller for the orbitals that should be aligned.
   * @param templateSystem    The system controller for the template orbitals.
   * @param exponent          The exponent used in the alignment. Has to be even.
   * @param kineticAlign      Use the kinetic energy for the orbital alignment.
   * @param replaceVirtuals   If true, the virtual valence orbitals are reconstructed to fit the IAO basis.
   */
  OrbitalAligner(std::shared_ptr<SystemController> systemController, std::shared_ptr<SystemController> templateSystem,
                 unsigned int exponent, bool kineticAlign, bool replaceVirtuals);
  virtual ~OrbitalAligner() = default;

  /**
   * @brief Align the given set of orbitals to the template. All orbitals will be aligned.
   * @param orbitals      The orbitals.
   * @param maxSweeps     Maximum number of alignment cycles before aborting.
   * @param orbitalRange  The range of orbitals to be rotated in.
   */
  virtual void localizeOrbitals(OrbitalController<SCFMode>& orbitals, unsigned int maxSweeps,
                                SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange) override final;

  /**
   * @brief Aligns two sets of orbitals.
   *
   *  The orbitals that will be rotated can be adjusted by the orbital range. This means that orbitals not included
   *  in this range will not be changed in the procedure even if they are included in the unpairedOrbitals vector.
   *
   * @param orbitals            The orbitals to be aligned.
   * @param unpairedOrbitals    The index list that defines which orbitals should be aligned. Entry of 1 will lead to
   * skipping the orbital.
   * @param unpairedOrbitalsRef The index list for the reference orbitals.
   * @param orbitalRange        Orbital subset to be rotated in.
   * @param maxSweeps           Maximum number of align iterations before aborting.
   * @param angleThreshold      Convergence threshold for the largest rotation angle.
   */
  void alignOrbitals(OrbitalController<SCFMode>& orbitals, const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitals,
                     const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitalsRef,
                     SpinPolarizedData<SCFMode, std::vector<unsigned int>> orbitalRange, unsigned int maxSweeps,
                     double angleThreshold = 1e-4);
  /**
   * @brief Aligns two sets of orbitals. Without restrictions to the set of orbitals in which is rotated.
   * @param orbitals            The orbitals to be aligned.
   * @param unpairedOrbitals    The index list that defines which orbitals should be aligned. Entry of 1 will lead to
   * skipping the orbital.
   * @param unpairedOrbitalsRef The index list for the reference orbitals.
   * @param maxSweeps           Maximum number of align iterations before aborting.
   * @param angleThreshold      Convergence threshold for the largest rotation angle.
   */
  void alignOrbitals(OrbitalController<SCFMode>& orbitals, const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitals,
                     const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitalsRef, unsigned int maxSweeps,
                     double angleThreshold = 1e-4);

 private:
  std::shared_ptr<SystemController> _system;
  std::shared_ptr<SystemController> _templateSystem;

  /**
   * @brief Create a vector that contains the indices of the reference orbitals for every orbital.
   * @param unpairedOrbitals    A list of 1|0 for the actual orbitals. If 1, the orbital will not be aligned.
   * @param unpairedOrbitalsRef A list of 1|0 for template orbitals. If 1, the orbital will be used as a template
   * orbital.
   * @return The mapping vector. The orbital is given by the vector row, the reference orbital index by the value.
   *         The value will be -1 if not paired.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXi>
  getReferenceOrbitals(const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitals,
                       const SpinPolarizedData<SCFMode, Eigen::VectorXi>& unpairedOrbitalsRef);

  /**
   * @brief The exponent used in the orbital alignment. This has the be even. Reasonable values are 2 or 4.
   */
  unsigned int _exponent;
  ///@brief Flag for usage of the orbital-wise kinetic energy in the Lagrangian
  bool _kineticAlign;
  /**
   * @brief Add an ij term to the Lagrangian and its derivatives.
   * @param Qii      Diagonal ii term.
   * @param Qij      Cross term.
   * @param Qjj      Diagonal jj term.
   * @param qi       Reference ii value.
   * @param qj       Reference jj value.
   * @param x        Current angle.
   * @param jRef     Reference index of j. Skipped if negative.
   * @param value    Current Lagrangian value.
   * @param gradient Current Lagrangian gradient.
   * @param hessian  Current Lagrangian hessian.
   */
  void addToLagrangian(const double Qii, const double Qjj, const double Qij, const double qi, const double qj,
                       const double x, const int jRef, double& value, double& gradient, double& hessian);
  /**
   * @brief Calculate orbital kinetic energies of the template system.
   * @return The orbital kinetic energies of the template system.
   */
  SpinPolarizedData<SCFMode, Eigen::VectorXd> getReferenceKineticEnergies();
  /**
   * If true, the virtual orbitals are localized. False by default.
   */
  bool _localizeVirtuals = false;
  bool _replaceVirtuals = false;
};

} /* namespace Serenity */

#endif /* ANALYSIS_ORBITALLOCALIZATION_ORBITALALIGNER_H_ */
