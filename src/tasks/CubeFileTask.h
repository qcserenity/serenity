/**
* @file   CubeFileTask.h
 *
 * @date   Nov 24, 2015
 * @author Jan Unsleber
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
#ifndef CUBEFILETASK_H_
#define CUBEFILETASK_H_

/* Include Serenity Internal Headers */
#include "settings/Options.h"
#include "settings/Reflection.h"
#include "data/SpinPolarizedData.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <map>
#include <memory>
#include <vector>

namespace Serenity {

class SystemController;
using namespace Serenity::Reflection;

struct CubeFileTaskSettings {
  CubeFileTaskSettings():
    cubeSpacing(0.1),
    cubeBorder(5.0),
    density(false),
    allOrbitals(false),
    occOrbitals(false),
    electrostaticPot(false),
    ntos(false),
    sedd(false),
    dori(false),
    elf(false),
    elfts(false),
    signedDensity(false),
    ntoPlotThreshold(0.1),
    orbitals(std::vector<unsigned int>(0)){
  };
  REFLECTABLE(
    (double) cubeSpacing,
    (double) cubeBorder,
    (bool) density,
    (bool) allOrbitals,
    (bool) occOrbitals,
    (bool) electrostaticPot,
    (bool) ntos,
    (bool) sedd,
    (bool) dori,
    (bool) elf,
    (bool) elfts,
    (bool) signedDensity,
    (double) ntoPlotThreshold,
    (std::vector<unsigned int>) orbitals
  )
};

/**
 * @class  CubeFileTask CubeFileTask.h
 * @brief  Write out the density on a cube file
 */
template<Options::SCF_MODES SCFMode>
class CubeFileTask : public Serenity::Task {
public:
  /**
   * @brief Constructor.
   * @param systems The system to work with, will be combined into one supersystem, its properties are plotted.
   * @param environmentSystems Environment systems to be added to the geometry.
   */
  CubeFileTask(const std::vector<std::shared_ptr<SystemController> >& systems,
               const std::vector<std::shared_ptr<SystemController> >& environmentSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~CubeFileTask() = default;

  void run();

  /**
   * @brief The settings/keywords for CubeFileTask. Set one of the following booleans
   *        to true to write the respective quantity to a cube file. Per default, all
   *        booleans are set to false. \n
   *        - density
   *        - allOrbitals
   *        - occOrbitals
   *        - electrostaticPot
   *        - NTOs (Natural Transition Orbitals)
   *        - sedd (Electron localization method)
   * 	    - dori (Electron localization method)
   *	    - signedDensity (density multiplied by the sign of the density Laplacian, for more SEDD.h)
   *	    - orbital (specify, which orbital to print)
   *	    - ntoPlotThreshold (threshold for which eigenvalues in the NTO evaluation are taken into account)
   *	    - cubeSpacing (Sets the step sizes for the three vectors spanning the cube mesh)
   * 	    - cubeBorder(Sets the border width around the geometry for which the cube grid is created)
   */
  CubeFileTaskSettings settings;

private:
  const std::vector<std::shared_ptr<SystemController> >  _systems;
  const std::vector<std::shared_ptr<SystemController> > _environmentSystems;
  SpinPolarizedData<SCFMode,std::string> _fnameSuffix;
};

} /* namespace Serenity */

#endif /* CUBEFILETASK_H_ */
