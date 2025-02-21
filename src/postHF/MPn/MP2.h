/**
 * @file   MP2.h
 *
 * @date   Jul 14, 2014
 * @author Jan Unsleber
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
#ifndef MP2_H_
#define MP2_H_
/* Include Std and External Headers */
#include <memory>

namespace Serenity {
/* Forward declarations */
class SystemController;
namespace Options {
enum class SCF_MODES;
}
/**
 * @class  MP2EnergyCorrector MP2.h
 * @brief  Calculates the MP2 correction of the electronic energy.
 */
template<Options::SCF_MODES SCFMode>
class MP2EnergyCorrector {
 public:
  /**
   * @brief Constructor.
   * @param systemController from which the electronic structure is / the orbitals are taken
   */
  MP2EnergyCorrector(std::shared_ptr<SystemController> systemController, const double ssScaling = 1.0,
                     const double osScaling = 1.0);
  /**
   * @brief Default destructor.
   */
  virtual ~MP2EnergyCorrector() = default;
  /**
   * @brief Calculates the MP2 energy.
   * @return Reurns the  MP2 energy.
   */
  double calculateElectronicEnergy();

 private:
  std::shared_ptr<SystemController> _systemController;
  const double _ssScaling;
  const double _osScaling;
};

} /* namespace Serenity */

#endif /* MP2_H_ */
