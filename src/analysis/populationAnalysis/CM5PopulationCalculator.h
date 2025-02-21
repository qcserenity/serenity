/**
 * @file CM5PopulationCalculator.h
 *
 * @date October 1, 2024
 * @author: Thorben Wiegmann
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
#ifndef ANALYSIS_CM5POPULATIONCALCULATOR_H_
#define ANALYSIS_CM5POPULATIONCALCULATOR_H_

/* Include Serenity Internal Headers */
#include "analysis/populationAnalysis/HirshfeldPopulationCalculator.h"
#include "data/SpinPolarizedData.h"
#include "parameters/Constants.h"
#include "settings/ElectronicStructureOptions.h"
#include "settings/Options.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace Serenity {
/* Forward declarations */
class SystemController;
class Atom;
/**
 * @class CM5PopulationCalculator CM5PopulationCalculator.h
 * @brief Performs a CM5 Population Analysis for a given system.
 *        [1] J. Chem. Theory Comput. 2012, 8, 2, 527–541; https://doi.org/10.1021/ct200866d
 */
template<Options::SCF_MODES SCFMode>
class CM5PopulationCalculator {
 public:
  /**
   * @brief Constructor
   * @param system The SystemController of the calculated system.
   * @param hirshFeld HirshfeldPopulationCalculator to calculate Hirshfeld charges.
   */
  CM5PopulationCalculator(std::shared_ptr<SystemController> system,
                          std::shared_ptr<HirshfeldPopulationCalculator<SCFMode>> hirshFeld);
  /**
   * @brief Default destructor
   */
  virtual ~CM5PopulationCalculator() = default;
  /**
   * @brief Calculates the atom populations (if necessary) and returns them.
   * @return The atom populations.
   */
  const SpinPolarizedData<SCFMode, Eigen::VectorXd>& getAtomPopulations() {
    if (!_atomPopulations) {
      calculateCM5Populations();
    }
    return *_atomPopulations;
  }

 private:
  // calculation of CM5 charges
  void calculateCM5Populations();
  Eigen::MatrixXd getPaulingBondOrderMatrix();
  Eigen::MatrixXd getParameterMatrix();
  std::shared_ptr<SystemController> _system;
  std::shared_ptr<HirshfeldPopulationCalculator<SCFMode>> _hirshFeld;
  std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>> _atomPopulations;
  std::vector<std::shared_ptr<Atom>> _atoms;
  unsigned int _nAtoms;
  std::vector<std::string> _atSymbols;
  // parameters from [1]
  double _alpha = 2.474;
  // D as given in the SI of [1]
  std::map<std::string, double> _atomicParameters = {
      {"H", 0.0056},   {"Li", 0},       {"Na", 0.0184},  {"K", 0.0130},   {"Rb", 0.0092},  {"Cs", 0.0065},
      {"Fr", 0.0046},  {"Be", 0.0333},  {"Mg", 0},       {"Ca", 0},       {"Sr", 0},       {"Ba", 0},
      {"Ra", 0},       {"Sc", 0},       {"Y", 0},        {"La", 0},       {"Ac", 0},       {"Ce", 0},
      {"Th", 0},       {"Pr", 0},       {"Pa", 0},       {"Nd", 0},       {"U", 0},        {"Pm", 0},
      {"Np", 0},       {"Sm", 0},       {"Pu", 0},       {"Eu", 0},       {"Am", 0},       {"Gd", 0},
      {"Cm", 0},       {"Tb", 0},       {"Bk", 0},       {"Dy", 0},       {"Cf", 0},       {"Ho", 0},
      {"Es", 0},       {"Er", 0},       {"Fm", 0},       {"Tm", 0},       {"Md", 0},       {"Yb", 0},
      {"No", 0},       {"Lu", 0},       {"Lr", 0},       {"Ti", 0},       {"Zr", 0},       {"Hf", 0},
      {"Rf", 0},       {"V", 0},        {"Nb", 0},       {"Ta", 0},       {"Db", 0},       {"Cr", 0},
      {"Mo", 0},       {"W", 0},        {"Sg", 0},       {"Mn", 0},       {"Tc", 0},       {"Re", 0},
      {"Bh", 0},       {"Fe", 0},       {"Ru", 0},       {"Os", 0},       {"Hs", 0},       {"Co", 0},
      {"Rh", 0},       {"Ir", 0},       {"Mt", 0},       {"Ni", 0},       {"Pd", 0},       {"Pt", 0},
      {"Ds", 0},       {"Cu", 0},       {"Ag", 0},       {"Au", 0},       {"Rg", 0},       {"Zn", 0},
      {"Cd", 0},       {"Hg", 0},       {"Cn", 0},       {"B", -0.1030},  {"Al", -0.0726}, {"Ga", -0.0512},
      {"In", -0.0361}, {"Tl", -0.0255}, {"Nh", -0.0179}, {"C", -0.0446},  {"Si", -0.0790}, {"Ge", -0.0557},
      {"Sn", -0.0393}, {"Pb", -0.0277}, {"Fl", -0.0195}, {"N", -0.1072},  {"P", -0.0756},  {"As", -0.0533},
      {"Sb", -0.0376}, {"Bi", -0.0265}, {"Mc", -0.0187}, {"O", -0.0802},  {"S", -0.0565},  {"Se", -0.0399},
      {"Te", -0.0281}, {"Po", -0.0198}, {"Lv", -0.0140}, {"F", -0.0629},  {"Cl", -0.0444}, {"Br", -0.0313},
      {"I", -0.0220},  {"At", -0.0155}, {"Tn", -0.0110}, {"He", -0.1543}, {"Ne", -0.1088}, {"Ar", -0.0767},
      {"Kr", -0.0541}, {"Xe", -0.0381}, {"Rn", -0.0269}, {"Og", -0.0189}};
  // special bond parameters as given in [1]
  std::map<std::string, double> _bondParameters = {{"HC", 0.0502},  {"HN", 0.1747},  {"HO", 0.1671},  {"CN", 0.0556},
                                                   {"CO", 0.0234},  {"NO", -0.0346}, {"CH", -0.0502}, {"NH", -0.1747},
                                                   {"OH", -0.1671}, {"NC", -0.0556}, {"OC", -0.0234}, {"ON", 0.0346}};
};
} // namespace Serenity

#endif /* ANALYSIS_CM5POPULATIONCALCULATOR_H_ */
