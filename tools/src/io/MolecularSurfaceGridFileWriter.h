/**
 * @file MolecularSurfaceGridFileWriter.h
 *
 * @author Moritz Bensberg
 * @date May 19, 2020
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

#ifndef IO_MOLECULARSURFACEGRIDFILEWRITER_H_
#define IO_MOLECULARSURFACEGRIDFILEWRITER_H_
/* Include Serenity Internal Headers */
#include "io/GeneralGridFileWriter.h" //DataOnGridWriter (base class)
/* Include Std and External Headers */
#include <memory>
#include <string>
#include <vector>

namespace Serenity {

class MolecularSurfaceController;
struct Settings;
/**
 * @class
 * @brief A class that writes values on the molecular surface to a xyz file.
 *        The element type in the xyz file will correspond to the element
 *        the surface point is associated to.
 */
class MolecularSurfaceGridFileWriter : public GeneralGridFileWriter {
 public:
  /**
   * @brief Constructor.
   * @param settings          The settings.
   * @param surfaceController The molecular surface grid controller.
   */
  MolecularSurfaceGridFileWriter(const Settings& settings, std::shared_ptr<MolecularSurfaceController> surfaceController);
  /**
   * @brief Default destructor.
   */
  virtual ~MolecularSurfaceGridFileWriter() = default;

 private:
  /**
   * @brief Writes the data to the cube file
   *
   * @param filename The base name of the cube file.
   * @param data The data which shall be printed to the cube file.
   */
  void writeData(std::string filename, const Eigen::VectorXd& data) override final;
  ///@brief The (surface) grid controller.
  std::shared_ptr<MolecularSurfaceController> _surfaceGridController;
};

} /* namespace Serenity */

#endif /* IO_MOLECULARSURFACEGRIDFILEWRITER_H_ */
