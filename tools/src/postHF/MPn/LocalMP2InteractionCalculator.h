/**
 * @file LocalMP2InteractionCalculator.h
 *
 * @date 18 Apr 2020
 * @author Moritz Bensberg
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

#ifndef POSTHF_MPN_LOCALMP2INTERACTIONCALCULATOR_H_
#define POSTHF_MPN_LOCALMP2INTERACTIONCALCULATOR_H_
/* Include Serenity Internal Headers */
#include "postHF/LocalCorrelation/LocalCorrelationController.h" //Default constructor of LocalCorrelationSettings
/* Include Std and External Headers */
#include <memory> //smrt_ptr
#include <vector> //std::vector

namespace Serenity {

/* Forward Declarations */
class SystemController;
class OrbitalPair;
class LocalCorrelationController;
/**
 * @class LocalMP2InteractionCalculator LocalMP2InteractionCalculator.h
 * @brief A class to calculate the MP2-correlation energy from intersubsystem orbital pairs.
 */
class LocalMP2InteractionCalculator {
 public:
  /**
   * @brief Constructor.
   * @param activeSystem The active system.
   * @param environmentSystems The environment systems.
   * @param lcSettings The local correlation settings.
   * @param maxCycles The max. number of iterations for the local-MP2 amplitude optimization.
   * @param maxResidual The convergence threshold for the local-MP2 amplitude optimization.
   * @param fullCoupling If true, the LMP2-correlation energy for the active orbital pairs is calculated in a
   * supersystem calculation to capture the effect of environment orbital pairs on the active-pair amplitudes.
   * @param supersystem The supersystem (optional, may be constructed on-the-fly)
   * @param environmentOrbitalIndices The environment orbital indices in the supersystem.
   */
  LocalMP2InteractionCalculator(std::shared_ptr<SystemController> activeSystem,
                                std::vector<std::shared_ptr<SystemController>> environmentSystems,
                                LocalCorrelationSettings lcSettings, unsigned int maxCycles, double maxResidual,
                                bool fullCoupling = false, std::shared_ptr<SystemController> supersystem = nullptr,
                                std::vector<unsigned int> environmentOrbitalIndices = {});
  /**
   * @brief Getter for the subsystem correlation-interaction (intersystem orbital-pairs).
   * @return The subsystem interaction.
   */
  double getInteractionEnergy();
  /**
   * @brief Getter for the correlation energy of the environment orbitals.
   * @return The correlation energy of the environment orbital pairs.
   */
  double getEnvironmentEnergy();
  /**
   * @brief Getter for the correlation energy of the active-system orbital pairs.
   * @return The correlation energy of the active system.
   */
  double getActiveEnergy();
  /**
   * @brief Getter for the coupling correction/influence of the environment orbital-pairs on the active-system orbital
   * pairs.
   * @return The coupling correction.
   */
  double getCouplingEnergyCorrection();
  /**
   * @brief Non-default Destructor. Cleans up tmp-system directories.
   */
  ~LocalMP2InteractionCalculator();

 private:
  std::shared_ptr<SystemController> _activeSystem;
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  LocalCorrelationSettings _lcSettings;
  unsigned int _maxCycles;
  double _maxResidual;
  bool _fullCoupling;
  std::shared_ptr<SystemController> _supersystem;
  void calculateLocalMP2Energies();
  void setUpSupersystem();
  void buildOrbitalPairs();
  std::shared_ptr<LocalCorrelationController> runLocalMP2();
  double calculateActiveOnlyLocalMP2Energy(std::shared_ptr<LocalCorrelationController> lcController);
  double getPairEnergy(std::shared_ptr<OrbitalPair> pair);
  bool _energiesAvailable = false;
  double _environmentEnergy;
  double _interactionEnergy;
  double _couplingCorrectionEnergy;
  double _activeEnergy;
  std::vector<std::shared_ptr<OrbitalPair>> _env_envPairs;
  std::vector<std::shared_ptr<OrbitalPair>> _env_actPairs;
  std::vector<std::shared_ptr<OrbitalPair>> _act_actPairs = {};
  std::vector<unsigned int> _environmentOrbitalIndices;
  void restrictCouplings(std::vector<std::shared_ptr<OrbitalPair>> pairs);
};

} /* namespace Serenity */

#endif /* POSTHF_MPN_LOCALMP2INTERACTIONCALCULATOR_H_ */
