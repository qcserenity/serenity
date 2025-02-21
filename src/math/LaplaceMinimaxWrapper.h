/**
 * @file   LaplaceMinimaxWrapper.h
 *
 * @date   Dec 26, 2021
 * @author Niklas Niemeyer
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

#ifndef LAPLACEMINIMAXWRAPPER_H
#define LAPLACEMINIMAXWRAPPER_H

/* Include Serenity Internal Headers */
#include "misc/SerenityError.h"
/* Include Std and External Headers */
#include <stdlib.h>
#include <Eigen/Dense>
#include <fstream>

namespace Serenity {

#ifdef SERENITY_USE_LAPLACE_MINIMAX

extern "C" void laplace_minimax_(double* errmax, double* xpnts, double* wghts, int* nlap, double* ymin, double* ymax,
                                 int* mxiter = 0, int* iprint = 0, double* stepmx = 0, double* tolrng = 0,
                                 double* tolpar = 0, double* tolerr = 0, double* delta = 0, double* afact = 0,
                                 double* do_rmsd = 0, double* do_init = 0, double* do_nlap = 0);

#endif /* SERENITY_USE_LAPLACE_MINIMAX */

static void getMinimaxRoots(Eigen::VectorXd& roots, Eigen::VectorXd& weights, double ymin, double ymax, double conv = 1e-6) {
#ifdef SERENITY_USE_LAPLACE_MINIMAX
  printBigCaption("Laplace-Minimax");
  printf("   Lower Bound       : %-6.2f\n", ymin);
  printf("   Upper Bound       : %-6.2f\n", ymax);

  std::string default_loc(std::getenv("SERENITY_HOME") ? std::getenv("SERENITY_HOME") : "");
  default_loc += "/build/lib/laplace-minimax/";
  std::ifstream file((default_loc + "/data/init_para.txt").c_str());
  if (file.good()) {
    setenv("LAPLACE_ROOT", default_loc.c_str(), true);
  }
  else {
    std::string lroot(std::getenv("LAPLACE_ROOT") ? std::getenv("LAPLACE_ROOT") : "");
    std::ifstream file((lroot + "/data/init_para.txt").c_str());
    if (!file.good()) {
      throw SerenityError("Cannot find laplace-minimax file 'init_para.txt'."
                          " Please set $LAPLACE_ROOT so that Serenity finds it in $LAPLACE_ROOT/data/init_para.txt.");
    }
  }

  int nPoints;
  double errmax = 0;
  bool converged = false;
  for (nPoints = 3; nPoints < 100; ++nPoints) {
    double xpnts[nPoints];
    double wghts[nPoints];

    // Make actual call to the library.
    laplace_minimax_(&errmax, xpnts, wghts, &nPoints, &ymin, &ymax);

    if (std::abs(errmax) < conv) {
      converged = true;
      roots = Eigen::Map<Eigen::VectorXd>(xpnts, nPoints);
      weights = Eigen::Map<Eigen::VectorXd>(wghts, nPoints);
      break;
    }
  }

  if (!converged) {
    throw SerenityError("laplace-minimax not converged!");
  }

  printf("   Quadrature Points : %-6i\n", nPoints);
  printf("   Abs(Max. Error)   : %-6.2e\n\n", std::abs(errmax));
#else
  // Avoid warnings if compiled without laplace-minimax.
  (void)roots;
  (void)weights;
  (void)ymin;
  (void)ymax;
  (void)conv;
  throw SerenityError("Link laplace-minimax with cmake -DSERENITY_USE_LAPLACE_MINIMAX=ON.");
#endif /* SERENITY_USE_LAPLACE_MINIMAX */
}
} /* namespace Serenity */

#endif /* LAPLACEMINIMAXWRAPPER_H */