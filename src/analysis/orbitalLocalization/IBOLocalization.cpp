/**
 * @file IBOLocalization.cpp
 *
 * @date Jun 15, 2016
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

/* Include Class Header*/
#include "analysis/orbitalLocalization/IBOLocalization.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisController.h"
#include "basis/AtomCenteredBasisControllerFactory.h"
#include "basis/Basis.h"
#include "basis/BasisFunctionProvider.h"
#include "data/matrices/CoefficientMatrix.h"
#include "integrals/wrappers/Libint.h"
#include "math/linearAlgebra/JacobiRotation.h"
#include "data/matrices/MatrixInBasis.h"
#include "integrals/OneElectronIntegralController.h"
#include "settings/Options.h"
#include "data/OrbitalController.h"
#include "system/SystemController.h"
#include "math/linearAlgebra/MatrixFunctions.h"


namespace Serenity {

template<Options::SCF_MODES SCFMode>
IBOLocalization<SCFMode>::IBOLocalization(std::shared_ptr<SystemController> systemController,
                                          bool IAOsOnly) :
  _system(systemController),
  _IAOsOnly(IAOsOnly){
};


template<Options::SCF_MODES SCFMode>
void IBOLocalization<SCFMode>::localizeOrbitals(OrbitalController<SCFMode>& orbitals,
    unsigned int maxSweeps) {
  /*
   * Most variables are named after their names in the
   *  paper describing the initial implementation.
   * Please refere to:
   * Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts
   * Gerald Knizia, J. Chem. Theory Comput., 2013, 9 (11), pp 4834–4843
   *
   * this link contains the reprint with corrected formulas in the appendix:
   * http://www.theochem.uni-stuttgart.de/~knizia/bin/iao_preprint.pdf
   */
  auto nOccOrbs = _system->getNOccupiedOrbitals<SCFMode>();
  // Create new basis
  std::shared_ptr<AtomCenteredBasisController> minaoBasis =
    AtomCenteredBasisControllerFactory::produce(
          _system->getGeometry(),
          _system->getSettings().basis.basisLibPath,
          _system->getSettings().basis.makeSphericalBasis,
          false,
          _system->getSettings().basis.firstECP,
          "MINAO");

  // store the basis
  _system->setBasisController(minaoBasis,Options::BASIS_PURPOSES::IAO_LOCALIZATION);

  // bases
  auto B1 = _system->getBasisController(Options::BASIS_PURPOSES::DEFAULT);
  auto B2 = _system->getBasisController(Options::BASIS_PURPOSES::IAO_LOCALIZATION);

  // coefficent matrix
  CoefficientMatrix<SCFMode> C = orbitals.getCoefficients();

  /*
   * Gather and calculate overlap integrals
   */
  auto& libint = Libint::getInstance();


  const auto& S1 = _system->getOneElectronIntegralController()->getOverlapIntegrals();
  Eigen::MatrixXd S2  = libint.compute1eInts(libint2::Operator::overlap,B2,B2);
  Eigen::MatrixXd S12 = libint.compute1eInts(libint2::Operator::overlap,B2,B1);

  /*
   * IAO calculation
   */
  // Coefficients in IAO basis
  Eigen::MatrixXd CIAO;
  Eigen::MatrixXd othoA;

  //  projection from basis 1 to basis 2 (P12)
  Eigen::MatrixXd P12 = S1.ldlt().solve(S12);
  Eigen::MatrixXd P21 = S2.ldlt().solve(S12.transpose());
  for_spin(C,nOccOrbs){

    /*
     * EQ 1
     */
    Eigen::MatrixXd Ct;
    Eigen::MatrixXd tmp = P12*P21*C_spin.leftCols(nOccOrbs_spin);

    // Symmetric orthogonalization
    {
      Eigen::MatrixXd tmpStmp = tmp.transpose()*S1*tmp;
//      Ct = tmpStmp.completeOrthogonalDecomposition().pseudoInverse().sqrt().eval();
      Ct = pseudoInversSqrt_Sym(tmpStmp);
      Ct = (tmp * Ct).eval();
    }

    /*
     * EQ 2
     */
    // slow but exact:
    Eigen::MatrixXd CCS = C_spin.leftCols(nOccOrbs_spin)*C_spin.leftCols(nOccOrbs_spin).transpose()*S1;
    Eigen::MatrixXd CtCtS = Ct*Ct.transpose()*S1;
    tmp = CCS*CtCtS*P12 + (Eigen::MatrixXd::Identity(CCS.cols(),CCS.rows())-CCS)
                          *(Eigen::MatrixXd::Identity(CtCtS.cols(),CtCtS.rows())-CtCtS)*P12;
    // alternative (faster, but approximate)
    //tmp = P12 + (C_spin.leftCols(nOccOrbs_spin)*C_spin.leftCols(nOccOrbs_spin).transpose()- Ct*Ct.transpose())*S12;

    // Symmetric orthogonalization
    { // these brackets make sure the solver is delete when its used and not needed
      Eigen::MatrixXd tmpStmp = tmp.transpose()*S1*tmp;
//      othoA = tmpStmp.completeOrthogonalDecomposition().pseudoInverse().sqrt().eval();
      othoA = pseudoInversSqrt_Sym(tmpStmp);
      othoA = (tmp * othoA).eval();
    }

    // transform occupied MOs in C into IAOs
    CIAO = othoA.transpose()*S1*C_spin.leftCols(nOccOrbs_spin);

    if (_IAOsOnly){
      /*
       * Stop here for IAOs
       */
      C_spin.leftCols(nOccOrbs_spin) = othoA*CIAO;
    } else {
      /*
       * Rotation from IAOs to IBOs
       */
      bool converged = false;
      unsigned int cycle = 0;
      while (!converged){
        ++cycle;

        // In the reference implementation of Knizia it is explained that
        //  Bij is effectively the gradient, it is therefore used as convergence crit.
        double gradient = 0.0;

        /*
         * 2x2 rotation loop
         */
        for (unsigned int i=0; i<nOccOrbs_spin; ++i){
          for (unsigned int j=0; j<i ; ++j){
            double Aij = 0;
            double Bij = 0;
            for (unsigned int atomIndex=0;atomIndex<_system->getGeometry()->getNAtoms() ; ++atomIndex){
              auto indices = minaoBasis->getBasisIndices();
              double Qii = 0;
              double Qij = 0;
              double Qjj = 0;
              for (unsigned int mu=indices[atomIndex].first; mu<indices[atomIndex].second; ++mu){
                Qii += CIAO(mu,i)*CIAO(mu,i);
                Qij += CIAO(mu,j)*CIAO(mu,i);
                Qjj += CIAO(mu,j)*CIAO(mu,j);
              }
              // p=2 simmilar to PM
              //Aij += 4.0*Qij*Qij - (Qii - Qjj)*(Qii - Qjj);
              //Bij += 4.0*Qij*(Qii - Qjj);
              // Alternative p=4
              const double Qii3 = Qii*Qii*Qii;
              const double Qjj3 = Qjj*Qjj*Qjj;
              Aij += - Qii3*Qii - Qjj3*Qjj;
              Aij += 6.0*(Qii*Qii + Qjj*Qjj)*Qij*Qij;
              Aij += Qii*Qjj3 + Qii3*Qjj;
              Bij += 4.0*Qij*(Qii3 - Qjj3);
            }

            gradient += Bij*Bij;

            // rotate orbital pair
            const double angle = 0.25*atan2(Bij,-Aij);
            JacobiRotation::rotate(CIAO.col(i),CIAO.col(j),angle);
          }
        }
        // Convergence Check
        if (fabs(gradient)<=1e-8) {
          std::cout << "    Converged after "<<cycle<<" orbital rotation cycles." << std::endl;
          C_spin.leftCols(nOccOrbs_spin) = othoA*CIAO;
          converged=true;
        } else if (cycle == maxSweeps) {
          // The paper claims that 5-10 sweeps should converge the orbitals
          //  therefore the procedure will stop after 200 cycles.
          // TODO: maybe this should throw a real ERROR
          std::cout << "    ERROR: IBO procedure did not converged after "<<cycle<<" orbital rotation cycles." << std::endl;
          std::cout << "           The orbitals will still be stored for error analysis." << std::endl;
          C_spin.leftCols(nOccOrbs_spin) = othoA*CIAO;
          converged=true;
        }
      }
    } /* rotations while loop */
  }; /* spin loop */

  orbitals.updateOrbitals(C,orbitals.getEigenvalues());
}

template class IBOLocalization<Options::SCF_MODES::RESTRICTED>;
template class IBOLocalization<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
