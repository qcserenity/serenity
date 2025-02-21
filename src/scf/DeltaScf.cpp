/**
 * @file DeltaScf.cpp
 *
 * @date Oct. 11, 2022
 * @author Niklas Goellmann
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
/* Include Class Header*/
#include "scf/DeltaScf.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "io/FormattedOutputStream.h"
#include "misc/SerenityError.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
std::shared_ptr<SPMatrix<SCFMode>> DeltaScf<SCFMode>::prepareMOMMatrix(const std::vector<int>& exca,
                                                                       const std::vector<int>& excb,
                                                                       std::shared_ptr<ElectronicStructure<SCFMode>> es) {
  std::shared_ptr<SPMatrix<SCFMode>> momMatrix = nullptr;
  if (!exca.empty() || !excb.empty()) {
    auto orbitalController = es->getMolecularOrbitals();
    momMatrix = std::make_shared<SPMatrix<SCFMode>>();
    auto& mom = (*momMatrix);
    const auto& C = orbitalController->getCoefficients();
    const auto nocc = es->getNOccupiedOrbitals();
    DensityMatrix<SCFMode> D(C.getBasisController());

    bool isAlpha = true;
    for_spin(D, C, mom, nocc) {
      const auto& exc = isAlpha ? exca : excb;
      mom_spin = C_spin;
      // do not swap occupations if exca or excb is just {0}
      if (!(exc.size() == 1 && exc[0] == 0)) {
        if (exc.size() % 2 != 0) {
          throw SerenityError("You have given a wrong number of orbital replacements in the exca/b vectors (SCFTask).");
        }
        if (isAlpha) {
          if (SCFMode == RESTRICTED) {
            OutputControl::nOut << " Swapped orbitals: " << std::endl;
          }
          else {
            OutputControl::nOut << " Swapped alpha orbitals: " << std::endl;
          }
        }
        else {
          OutputControl::nOut << " Swapped beta orbitals: " << std::endl;
        }
        for (unsigned i = 0; i < exc.size(); i += 2) {
          mom_spin.col(nocc_spin - 1 - abs(exc[i])).swap(mom_spin.col(nocc_spin + abs(exc[i + 1])));
          OutputControl::n.printf(" HOMO - %-2i <--> LUMO + %-2i\n", abs(exc[i]), abs(exc[i + 1]));
        }
      }
      mom_spin = mom_spin.leftCols(nocc_spin).eval();
      D_spin = (SCFMode == Options::SCF_MODES::RESTRICTED ? 2.0 : 1.0) * mom_spin * mom_spin.transpose();
      isAlpha = false;
    };
    es->getDensityMatrixController()->setDensityMatrix(D);
  }
  return momMatrix;
}

template class DeltaScf<Options::SCF_MODES::RESTRICTED>;
template class DeltaScf<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
