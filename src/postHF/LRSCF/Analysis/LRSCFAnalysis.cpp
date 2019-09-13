/**
 * @file LRSCFAnalysis.cpp
 * @author: Michael Boeckers
 *
 * @date Dec. 18, 2018
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
/* Include Class Header*/
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
void LRSCFAnalysis<SCFMode>::printDominantContributions(
            const std::vector<std::shared_ptr<LRSCFController<SCFMode> > >& lrscf,
            const std::vector<Eigen::MatrixXd >& eigenvectors,
            const Eigen::VectorXd& eigenvalues,
            const double th){
  
  Eigen::MatrixXd x,y;
  x = eigenvectors[0];
  if(eigenvectors.size() == 2) {
    y = eigenvectors[1];
  } else {
    //For TDA, y is zero and not calculated.
    y = Eigen::MatrixXd::Zero(eigenvectors[0].rows(),eigenvectors[0].cols());
  }

  Eigen::MatrixXd c = 100 * (x.array().square() + y.array().square());

  printf("---------------------------------------------------------------------------------------\n");
  printf("                               Dominant contributions                                  \n");
  printf("---------------------------------------------------------------------------------------\n");
  printf(" %5s %27s %24s %5s %9s %11s\n","state","energy","sys","i","a","|c|^2*100");
  printf(" %15s %10s %10s %12s\n","(a.u.)","(eV)","(nm)","(cm^-1)");
  printf("---------------------------------------------------------------------------------------\n");
  for (unsigned int iState = 0; iState < eigenvalues.rows(); ++iState) {
    printf(" %5i %10.6f %10.4f %10.2f %12.2f ",
        iState + 1,
        eigenvalues(iState),
        eigenvalues(iState) * HARTREE_TO_EV,
        HARTREE_TO_NM / eigenvalues(iState),
        eigenvalues(iState) * HARTREE_TO_OOCM);
    unsigned int iCount = 0;
    bool first = true;
    for (unsigned int I = 0; I < lrscf.size(); ++I) {
      auto nOcc = lrscf[I]->getNOccupied();
      auto nVirt = lrscf[I]->getNVirtual();
      unsigned int iSpin = 0;
      for_spin(nOcc,nVirt) {
        std::string spin = (iSpin == 0) ? "a" : "b";
        for (unsigned int j = 0, jb = iCount; j <  nOcc_spin; ++j) {
          for (unsigned int b = nOcc_spin; b < nOcc_spin + nVirt_spin; ++b, ++jb) {
            if (c(jb,iState) > th) {
              if (first) {
                printf("%5i %5i %2s %5i %2s %7.2f \n",I+1,j+1,spin.c_str(),b+1,spin.c_str(),c(jb,iState));
                first = false;
              } else {
                printf("%58i %5i %2s %5i %2s %7.2f \n",I+1,j+1,spin.c_str(),b+1,spin.c_str(),c(jb,iState));
              }
            }
            iCount += 1;
          }
        }
        iSpin += 1;
      };
      //Break after first cycle in uncoupled case
      if (lrscf.size() == 1) break;
    }
    if (first) printf("\n");
  }
}

template class LRSCFAnalysis<Options::SCF_MODES::RESTRICTED>;
template class LRSCFAnalysis<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
