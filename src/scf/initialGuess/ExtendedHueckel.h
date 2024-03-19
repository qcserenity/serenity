/**
 * @file   ExtendedHueckel.h
 *
 * @date   Nov 23, 2013
 * @author Thomas Dresselhaus
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
#ifndef EXTENDEDHUECKEL_H_
#define EXTENDEDHUECKEL_H_
/* Include Serenity Internal Headers */
#include "scf/initialGuess/InitialGuessCalculator.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <map>
#include <string>
#include <vector>

namespace Serenity {
/* Forward declarations */
class AtomCenteredBasisController;
template<Options::SCF_MODES>
class MatrixInBasis;
/**
 * @class ExtendedHueckel ExtendedHueckel.h
 * @brief Implementation of an extended Hueckel calculator.
 *
 * Should be according to:
 * Roald Hoffmann, J. Chem. Phys. 39, 1397 (1963); http://dx.doi.org/10.1063/1.1734456
 */
class ExtendedHueckel : public InitialGuessCalculator<Options::SCF_MODES::RESTRICTED> {
 public:
  /*
   * Although all methods and fields are static, instantiation is allowed to suit the interface.
   */
  ExtendedHueckel() = default;
  virtual ~ExtendedHueckel() = default;

  std::unique_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>
  calculateInitialGuess(std::shared_ptr<SystemController> systemController) override final;

 private:
  /**
   * @param minimalBasisController the (shared) ownership will be transferred to the resulting orbitals.
   * @param minimalBasisOverlaps
   * @returns Hueckel orbitals corresponding to the geometry. They can only be defined in a basis.
   *          for which parameters are actually available, typically a minimal basis.
   */
  static std::unique_ptr<OrbitalController<Options::SCF_MODES::RESTRICTED>>
  calculateHueckelOrbitals(std::shared_ptr<SystemController> systemController,
                           std::shared_ptr<AtomCenteredBasisController> minimalBasisController,
                           const MatrixInBasis<RESTRICTED>& minimalBasisOverlaps);

  static constexpr double _k = 1.75;

  /*
   * The parameters for the diagonal elements of the Hückel matrix.
   * Each parameter shall be used for one shell of basis functions.
   */
  static const std::map<std::string, std::vector<double>> _parameters;

  /*
   * Function to set the parameters for the elements.
   */
  static std::map<std::string, std::vector<double>> assignParameters() {
    std::map<std::string, std::vector<double>> parameters;
    std::vector<double> paramsForElement;

    // EHT parameters
    // using AO energies from http://www.few.vu.nl/~visscher/FiniteNuclei/Table3.html
    //   (Nov. 8, 2017) author Lucas Visscher

    // First period
    paramsForElement.resize(1);
    // Hydrogen
    paramsForElement[0] = -5.0000666E-01;
    parameters["H"] = paramsForElement;
    // Helium
    paramsForElement[0] = -9.1799069E-01;
    parameters["He"] = paramsForElement;

    // Second period
    paramsForElement.resize(3);
    // Lithium
    paramsForElement[0] = -2.4779791E+00;
    paramsForElement[1] = -1.9633887E-01;
    paramsForElement[2] = 0.0;
    parameters["Li"] = paramsForElement;
    // Berylium
    paramsForElement[0] = -4.7334980E+00;
    paramsForElement[1] = -3.0932208E-01;
    paramsForElement[2] = 0.0;
    parameters["Be"] = paramsForElement;
    // Bor
    paramsForElement[0] = -7.6976454E+00;
    paramsForElement[1] = -4.9491914E-01;
    paramsForElement[2] = -3.0972112E-01;
    parameters["B"] = paramsForElement;
    // Carbon
    paramsForElement[0] = -1.1343601E+01;
    paramsForElement[1] = -7.1260267E-01;
    paramsForElement[2] = -4.0664220E-01;
    parameters["C"] = paramsForElement;
    // Nitrogen
    paramsForElement[0] = -1.5676414E+01;
    paramsForElement[1] = -9.6478811E-01;
    paramsForElement[2] = -5.0818404E-01;
    parameters["N"] = paramsForElement;
    // Oxygen
    paramsForElement[0] = -2.0698561E+01;
    paramsForElement[1] = -1.2524124E+00;
    paramsForElement[2] = -6.1537304E-01;
    parameters["O"] = paramsForElement;
    // Fluorine
    paramsForElement[0] = -2.6411757E+01;
    paramsForElement[1] = -1.5759831E+00;
    paramsForElement[2] = -7.2866290E-01;
    parameters["F"] = paramsForElement;
    // Neon
    paramsForElement[0] = -3.2817472E+01;
    paramsForElement[1] = -1.9358461E+00;
    paramsForElement[2] = -8.4826678E-01;
    parameters["Ne"] = paramsForElement;
    return parameters;
  }
};

} /* namespace Serenity */
#endif /* EXTENDEDHUECKEL_H_ */
