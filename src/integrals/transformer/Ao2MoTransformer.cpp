/**
 * @file Ao2MoTransformer.cpp
 *
 * @date Apr 29, 2016
 * @author Michael Boeckers
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
#include "integrals/transformer/Ao2MoTransformer.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"


namespace Serenity {

Ao2MoTransformer::Ao2MoTransformer(std::shared_ptr<BasisController> basisController):
  _nBasisFunc(basisController->getNBasisFunctions()){

}


//ToDo: This routine works for Cholesky vectors calculated with eigen. However, since we
//      very soon have our own decomposition function (already present locally on my machine),
//      it needs to be revised.
//      (The Transpositions are no longer needed since the Cholesky vectors are stored in
//      the correct oreder)
//template<> void Ao2MoTransformer<Options::SCF_MODES::RESTRICTED>::transformCholeskyVectors(
//    Matrix<double>& choleskyVectors,
//    Transpositions<-1, -1, int>& P) {
//  /*
//   * Transforms cholesky vectors from AO to MO basis.
//   *
//   * ToDo:
//   * The transformation is currently parallelized over the columns of the cholesky
//   * vectors. Since the dimension of the cholesky vectors is usually smaller than the
//   * number of basis functions, it could be more efficient to parallelize over
//   * the second loop. This should be tested.
//   */
//  //Permutate rows, so that cholesky vectors are multiplied with correct coefficients
//  choleskyVectors = P.transpose() * choleskyVectors ;
//  //Transformation 1: L_{\mu j}^K = \sum_nu c_{\nu j} L_{\mu\nu}^K
//  Matrix<double> tmpCholeskyVectors(choleskyVectors.rows(), choleskyVectors.cols());
//  tmpCholeskyVectors.setZero();
//#pragma omp parallel shared(tmpCholeskyVectors,choleskyVectors)
//  {
//#pragma omp for schedule(dynamic)
//  for (unsigned int K = 0; K < choleskyVectors.cols(); ++K) {
//    for (unsigned int j = 0, muj = 0; j < _nBasisFunc; ++j) {
//      for (unsigned int mu = 0, munu = 0; mu < _nBasisFunc; ++mu, ++muj) {
//        for (unsigned int nu = 0; nu < _nBasisFunc; ++nu, ++munu) {
//          tmpCholeskyVectors(muj, K) += _coeffMat(nu, j) * choleskyVectors(munu, K);
//        }/* AO index nu */
//      }/* AO index mu */
//    }/* MO index j */
//  }/* cholesky vector column K */
//  //Transformation 2: L_{ij}^K = \sum_nu c_{\mu i} L_{\mu j}^K
//  choleskyVectors.setZero();
//#pragma omp for schedule(dynamic)
//  for (unsigned int K = 0; K < choleskyVectors.cols(); ++K) {
//    for (unsigned int i = 0, ij = 0; i < _nBasisFunc; ++i) {
//      for (unsigned int j = 0, muj = 0; j < _nBasisFunc; ++j, ++ij) {
//        for (unsigned int mu = 0; mu < _nBasisFunc; ++mu, ++muj) {
//          choleskyVectors(ij, K) += _coeffMat(mu, i) * tmpCholeskyVectors(muj, K);
//        }/* AO index mu */
//      }/* MO index j */
//    }/* MO index i */
//  }/* cholesky vector column K */
//  }/* omp parrallel */
//  //Permute back to old order so that CD_Integral calculator can work with this transformed matrix
//  choleskyVectors = P * choleskyVectors;
//}

void Ao2MoTransformer::transformTwoElectronIntegrals(
    RegularRankFourTensor<double>& twoElectronIntegrals,
    RegularRankFourTensor<double>& result,
    const Eigen::MatrixXd coefficients,
    const unsigned int& nOrbitals) {
  assert(coefficients.rows()==_nBasisFunc);
  assert(result.getKOffset()==nOrbitals);
  /*
   * Function copied and adapted from MP2.cpp
   *
   * Transforms two electron integrals from AO to MO basis. Scales with N^5
   *
   * Transformation 1 (a,b||i,j) -> (a,b||i,kappa)
   */
  RegularRankFourTensor<double> tempRegularRankFourTensor(_nBasisFunc, 0.0);
#pragma omp parallel shared(tempRegularRankFourTensor,twoElectronIntegrals)
  {
#pragma omp for schedule(dynamic)
    for (unsigned int a = 0; a < _nBasisFunc; ++a) {
      for (unsigned int b = 0; b < _nBasisFunc; ++b) {
        for (unsigned int i = 0; i < _nBasisFunc; ++i) {

          for (unsigned int kappa = 0; kappa < nOrbitals; ++kappa) {

            tempRegularRankFourTensor(a, b, i, kappa) = 0.0;
            for (unsigned int j = 0; j < _nBasisFunc; ++j) {
              tempRegularRankFourTensor(a, b, i, kappa) += coefficients(j, kappa)
                  * twoElectronIntegrals(a, b, i, j);
            }
          }
        }
      }
    }
    /*
     * Transformation 2 (a,b||i,kappa) -> (a,b||iota,kappa)
     */
#pragma omp for schedule(dynamic)
    for (unsigned int a = 0; a < _nBasisFunc; ++a) {
      for (unsigned int b = 0; b < _nBasisFunc; ++b) {

        for (unsigned int kappa = 0; kappa < nOrbitals; ++kappa) {
          for (unsigned int iota = 0; iota < nOrbitals; ++iota) {

            twoElectronIntegrals(a, b, iota, kappa) = 0.0;
            for (unsigned int i = 0; i < _nBasisFunc; ++i) {
              twoElectronIntegrals(a, b, iota, kappa) += coefficients(i, iota)
                  * tempRegularRankFourTensor(a, b, i, kappa);
            }
          }
        }
      }
    }

    /*
     * Transformation 3 (a,b||iota,kappa) -> (a,beta||iota,kappa)
     */
#pragma omp for schedule(dynamic)
    for (unsigned int a = 0; a < _nBasisFunc; ++a) {

      for (unsigned int kappa = 0; kappa < nOrbitals; ++kappa) {
        for (unsigned int iota = 0; iota < nOrbitals; ++iota) {
          for (unsigned int beta = 0; beta < nOrbitals; ++beta) {

            tempRegularRankFourTensor(a, beta, iota, kappa) = 0.0;
            for (unsigned int b = 0; b < _nBasisFunc; ++b) {
              tempRegularRankFourTensor(a, beta, iota, kappa) += coefficients(b, beta)
                  * twoElectronIntegrals(a, b, iota, kappa);
            }
          }
        }
      }
    }

    /*
     * Transformation 4 (a,beta||iota,kappa) -> (alpha,beta||iota,kappa)
     */
#pragma omp for schedule(dynamic)
    for (unsigned int kappa = 0; kappa < nOrbitals; ++kappa) {
      for (unsigned int iota = 0; iota < nOrbitals; ++iota) {
        for (unsigned int beta = 0; beta < nOrbitals; ++beta) {
          for (unsigned int alpha = 0; alpha < nOrbitals; ++alpha) {

            result(alpha, beta, iota, kappa) = 0.0;
            for (unsigned int a = 0; a < _nBasisFunc; ++a) {
              result(alpha, beta, iota, kappa) += coefficients(a, alpha)
                  * tempRegularRankFourTensor(a, beta, iota, kappa);
            }
          }
        }
      }
    }
  }
}

}
