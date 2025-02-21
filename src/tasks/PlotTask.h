/**
 * @file   PlotTask.h
 *
 * @date   Nov 24, 2015
 * @author Jan Unsleber, Anja Massolle
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
#ifndef PLOTTASK_H_
#define PLOTTASK_H_

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h" // Output suffix
#include "tasks/Task.h"             // Task
/* Include Std and External Headers */
#include <limits> // std::numeric_limits
namespace Serenity {

class SystemController;
using namespace Serenity::Reflection;

struct PlotTaskSettings {
  PlotTaskSettings()
    : p1({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity()}),
      p2({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity()}),
      p3({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity()}),
      p4({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity()}),
      atom1(std::numeric_limits<int>::infinity()),
      atom2(std::numeric_limits<int>::infinity()),
      atom3(std::numeric_limits<int>::infinity()),
      atom4(std::numeric_limits<int>::infinity()),
      gridSpacing({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()}),
      borderWidth(2.0),
      xUnitVector({1.0, 0.0, 0.0}),
      yUnitVector({0.0, 1.0, 0.0}),
      zUnitVector({0.0, 0.0, 1.0}),
      projectCutOffRadius(3.0),
      xyHeatmap(false),
      density(false),
      allOrbitals(false),
      occOrbitals(false),
      electrostaticPot(false),
      sedd(false),
      dori(false),
      elf(false),
      elfts(false),
      signedDensity(false),
      orbitals(std::vector<unsigned int>(0)),
      maxGridPoints(1e+8),
      cavity(false),
      gridCoordinates(false),
      ntos(false),
      ntoPlotThreshold(0.1),
      excitations({0}),
      nros(false),
      nrominimum(0.75),
      cctrdens(false),
      ccexdens(false) {
  }
  REFLECTABLE((std::vector<double>)p1, (std::vector<double>)p2, (std::vector<double>)p3, (std::vector<double>)p4,
              (int)atom1, (int)atom2, (int)atom3, (int)atom4, (std::vector<double>)gridSpacing, (double)borderWidth,
              (std::vector<double>)xUnitVector, (std::vector<double>)yUnitVector, (std::vector<double>)zUnitVector,
              (double)projectCutOffRadius, (bool)xyHeatmap, (bool)density, (bool)allOrbitals, (bool)occOrbitals,
              (bool)electrostaticPot, (bool)sedd, (bool)dori, (bool)elf, (bool)elfts, (bool)signedDensity,
              (std::vector<unsigned int>)orbitals, (double)maxGridPoints, (bool)cavity, (bool)gridCoordinates,
              (bool)ntos, (double)ntoPlotThreshold, (std::vector<unsigned int>)excitations, (bool)nros,
              (double)nrominimum, (bool)cctrdens, (bool)ccexdens)
};

/**
 * @class  PlotTask PlotTask.h
 * @brief  Write out the density on a file
 */
template<Options::SCF_MODES SCFMode>
class PlotTask : public Serenity::Task {
 public:
  /**
   * @brief Constructor.
   * @param systems The system to work with, will be combined into one
   * supersystem, its properties are plotted.
   * @param environmentSystems Environment systems to be added to the geometry.
   */
  PlotTask(const std::vector<std::shared_ptr<SystemController>>& systems,
           const std::vector<std::shared_ptr<SystemController>>& environmentSystems, std::string filename = "");
  /**
   * @brief Default destructor.
   */
  virtual ~PlotTask() = default;

  /**
   * @see Task
   */
  void run();

  /**
   * @brief The settings/keywords for PlotTask. Set one of the following
   *booleans to true to write the respective quantity to a file. Per
   *default, all booleans are set to false. \n
   *        - density
   *        - allOrbitals
   *        - occOrbitals
   *	      - orbital (specify, which orbitals to print)
   *        - electrostaticPot
   *        - sedd  (Electron localization method)
   * 	      - dori  (Electron localization method)
   *	      - signedDensity (density multiplied by the sign of the density
   *Laplacian, for more SEDD.h)
   *        - ELF   (Electron localization method)
   *        - ELFTS (Electron localization method)
   *        - xyHeatmap \n\n
   * In addition the following settings can be set: \n
   * 	      - borderWidth (default 2.0) Sets the border width in A around the geometry for which
   *the grid is created
   *        - p1, p2, p3 The points which define the plane
   *        - atom1, atom2, atom3 The atoms which defines the plane
   *        - gridSpacing (default {0.25 0.25 0.25} 3D, {0.1 0.1 0.1} 2D) Sets the step size in A in x, y and z direction
   *        - xUnitVector, yUnitVector, zUnitVector The unit vectors of the cube grid
   *        - projectCutOffRadius (default 3.0) The maximal distance in A an atom in a 2D plot is
   * allowed to have to be projected on the plane for the evaluation of min / max values for the generation of the grid
   *        - maxGridPoints (default 1e+8) Maximum number of grid points.
   *        - cavity (default false) plot values on the molecular cavity.
   *        - gridCoordinates plot the gird coordinates and weights to file.
   *        - NTOs (Natural Transition Orbitals)
   *	    - ntoPlotThreshold (threshold for which eigenvalues in the NTO evaluation are taken into account)
   */
  PlotTaskSettings settings;

 private:
  ///@brief true if a cube grid should be generated
  bool _cubeTask = false;
  ///@brief true if a planar grid should be generated
  bool _planeTask = false;
  ///@brief true if a linear grid should be generated
  bool _lineTask = false;
  ///@brief true if a grid with a single point should be generated
  bool _pointTask = false;
  ///@brief true if the molecular surface grid should be used for plotting.
  bool _cavityGrid = false;
  ///@brief systems for which the property should be evaluated
  const std::vector<std::shared_ptr<SystemController>> _systems;
  ///@brief systems for which the property should NOT be evaluated, they
  /// are only considered for the generation of the grid
  const std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  ///@brief name of the cube/dat file
  std::string _filename;
  ///@brief alpha / beta suffix of the filename for an unrestricted calculation
  SpinPolarizedData<SCFMode, std::string> _fnameSuffix;
};

} /* namespace Serenity */

#endif /* PLOTTASK_H_ */
