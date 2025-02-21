/**
 * @file LRSCFAnalysis.cpp
 * @author: Michael Boeckers
 *
 * @date Dec. 18, 2018
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
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"
/* Include Serenity Internal Headers */
#include "parameters/Constants.h"
#include "postHF/LRSCF/LRSCFController.h"
#include "settings/LRSCFOptions.h"
#include "tasks/LRSCFTask.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void LRSCFAnalysis<SCFMode>::printDominantContributions(const std::vector<std::shared_ptr<LRSCFController<SCFMode>>>& lrscf,
                                                        const std::vector<Eigen::MatrixXd>& eigenvectors,
                                                        const Eigen::VectorXd& eigenvalues, const double contribThresh) {
  auto settings = lrscf[0]->getLRSCFSettings();
  bool isNotCC2 = settings.method == Options::LR_METHOD::TDA || settings.method == Options::LR_METHOD::TDDFT;

  Eigen::MatrixXd contributions = eigenvectors[0].cwiseProduct(eigenvectors[1]);

  // This is a simple helper struct to identify orbital transitions.
  struct Excitation {
    // The subsystem.
    unsigned I;
    // The occupied orbital.
    unsigned i;
    // The virtual orbital.
    unsigned a;
    // The spin.
    char spin;
  };

  std::function<Excitation(unsigned)> getIndices = [&](unsigned ia) {
    Excitation exc;

    unsigned iCount = 0;
    for (unsigned I = 0; I < lrscf.size(); ++I) {
      auto no = lrscf[I]->getNOccupied();
      auto nv = lrscf[I]->getNVirtual();

      bool isAlpha = true;
      for_spin(no, nv) {
        for (unsigned i = 0; i < no_spin; ++i) {
          for (unsigned a = no_spin; a < no_spin + nv_spin; ++a, ++iCount) {
            if (iCount == ia) {
              exc.I = I;
              exc.i = i;
              exc.a = a;
              exc.spin = isAlpha ? 'a' : 'b';
            }
          }
        }
        isAlpha = false;
      };
    }
    return exc;
  };

  printf("---------------------------------------------------------------------------------------\n");
  if (SCFMode == Options::SCF_MODES::RESTRICTED) {
    if (settings.method == Options::LR_METHOD::TDA) {
      if (settings.scfstab == Options::STABILITY_ANALYSIS::NONE) {
        printf("                                    TDA Summary                                      \n\n");
      }
      else {
        printf("                           SCF Stability Analysis Summary                            \n\n");
      }
    }
    else if (settings.method == Options::LR_METHOD::TDDFT) {
      printf("                                   TDDFT Summary                                     \n\n");
    }
    else if (settings.method == Options::LR_METHOD::CC2) {
      printf("                                    CC2 Summary                                      \n\n");
    }
    else if (settings.method == Options::LR_METHOD::CISDINF) {
      printf("                                 CIS(Dinf) Summary                                   \n\n");
    }
    else if (settings.method == Options::LR_METHOD::ADC2) {
      printf("                                   ADC(2) Summary                                    \n\n");
    }
    else if (settings.method == Options::LR_METHOD::CISD) {
      printf("                                   CIS(D) Summary                                    \n\n");
    }
  }
  else {
    if (settings.method == Options::LR_METHOD::TDA) {
      if (settings.scfstab == Options::STABILITY_ANALYSIS::NONE) {
        printf("                                  uTDA Summary                                       \n\n");
      }
      else {
        printf("                           SCF Stability Analysis Summary                            \n\n");
      }
    }
    else if (settings.method == Options::LR_METHOD::TDDFT) {
      printf("                                 uTDDFT Summary                                      \n\n");
    }
    else if (settings.method == Options::LR_METHOD::CC2) {
      printf("                                  uCC2 Summary                                       \n\n");
    }
    else if (settings.method == Options::LR_METHOD::CISDINF) {
      printf("                               uCIS(Dinf) Summary                                    \n\n");
    }
    else if (settings.method == Options::LR_METHOD::ADC2) {
      printf("                                 uADC(2) Summary                                     \n\n");
    }
    else if (settings.method == Options::LR_METHOD::CISD) {
      printf("                                   uCIS(D) Summary                                   \n\n");
    }
  }
  printf("---------------------------------------------------------------------------------------\n");
  if (SCFMode == Options::SCF_MODES::RESTRICTED && settings.scfstab == Options::STABILITY_ANALYSIS::NONE) {
    if (settings.triplet) {
      printf("                                Triplet Excitations.                                 \n\n");
    }
    else {
      printf("                                Singlet Excitations.                                 \n\n");
    }
  }
  printf("    Convergence Threshold            : %16.1e\n", settings.conv);
  if (!isNotCC2) {
    if (settings.ltconv != 0) {
      printf("    Laplace Transformation Threshold : %16.1e\n", settings.ltconv);
    }
    if (settings.sss != 1.0 || settings.oss != 1.0) {
      printf("    Same-Spin Scaling                : %16.3f\n", settings.sss);
      printf("    Opposite-Spin Scaling            : %16.3f\n", settings.oss);
    }
  }

  printf("---------------------------------------------------------------------------------------\n");
  printf("                               Dominant Contributions                                  \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" %5s %27s %24s %5s %9s %11s\n", "state", "energy", "sys", "i", "a", "|c|^2*100");
  printf(" %15s %10s %10s %12s\n", "(a.u.)", "(eV)", "(nm)", "%singles");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned iEigen = 0; iEigen < eigenvalues.size(); ++iEigen) {
    if (isNotCC2 && contributions.col(iEigen).sum() < 0.95) {
      printf("     Normalization of the following eigenvector is probably faulty. Be careful.\n");
    }
    // sort orbital transition by their squared coefficients
    std::vector<std::pair<double, unsigned>> coeffs(0);
    for (unsigned ia = 0; ia < contributions.rows(); ++ia) {
      coeffs.push_back(std::make_pair(contributions(ia, iEigen), ia));
    }
    std::stable_sort(coeffs.begin(), coeffs.end());
    std::reverse(coeffs.begin(), coeffs.end());

    printf(" %4i %11.6f %10.5f %10.2f %10.2f ", iEigen + 1, eigenvalues(iEigen), eigenvalues(iEigen) * HARTREE_TO_EV,
           HARTREE_TO_NM / eigenvalues(iEigen), contributions.col(iEigen).sum() * 100);
    double sum = 0.0;
    unsigned index = 0;
    do {
      unsigned ia = coeffs[index].second;
      auto exc = getIndices(ia);
      double contribution = coeffs[index].first;
      const char* str = (index == 0) ? " %6i %5i %2c %5i %2c %8.2f \n" : "%58i %5i %2c %5i %2c %8.2f \n";
      printf(str, exc.I + 1, exc.i + 1, exc.spin, exc.a + 1, exc.spin, 100 * contribution);
      sum += contribution / contributions.col(iEigen).sum();
      ++index;
    } while (sum < contribThresh && index < coeffs.size());
  }
  printf("---------------------------------------------------------------------------------------\n");
}

template class LRSCFAnalysis<Options::SCF_MODES::RESTRICTED>;
template class LRSCFAnalysis<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
