/**
 * @file   XCFun.cpp
 *
 * @date   Feb 3, 2017
 * @author M. Boeckers
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
#include "dft/functionals/wrappers/XCFun.h"
/* Include Serenity Internal Headers */
#include "grid/GridController.h"
#include "misc/Timing.h"
/* Include Std and External Headers */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


namespace Serenity {

template<Options::SCF_MODES T>
XCFun<T>::XCFun(unsigned int maxBlockSize):
  _maxBlockSize(maxBlockSize){
}

template<Options::SCF_MODES T>
char* XCFun<T>::getAlias(BASIC_FUNCTIONALS functional) {
  char* alias = (char*) "null";
  switch(functional) {
      case BASIC_FUNCTIONALS::NONE:
        break;
      case BASIC_FUNCTIONALS::SLATERX:
        alias = (char*) "slaterx";
        break;
      case BASIC_FUNCTIONALS::PW86X:
        alias = (char*) "pw86x";
        break;
      case BASIC_FUNCTIONALS::VWN3C:
        alias = (char*) "vwn3c";
        break;
      case BASIC_FUNCTIONALS::VWN5C:
        alias = (char*) "vwn5c";
        break;
      case BASIC_FUNCTIONALS::PBEC:
        alias = (char*) "pbec";
        break;
      case BASIC_FUNCTIONALS::PBEX:
        alias = (char*) "pbex";
        break;
      case BASIC_FUNCTIONALS::BECKEX:
        alias = (char*) "beckex";
        break;
      case BASIC_FUNCTIONALS::BECKECORRX:
        alias = (char*) "beckecorrx";
        break;
      case BASIC_FUNCTIONALS::BECKESRX:
        alias = (char*) "beckesrx";
        break;
      case BASIC_FUNCTIONALS::BECKECAMX:
        alias = (char*) "beckecamx";
        break;
      case BASIC_FUNCTIONALS::BRX:
        alias = (char*) "brx";
        break;
      case BASIC_FUNCTIONALS::BRC:
        alias = (char*) "brc";
        break;
      case BASIC_FUNCTIONALS::BRXC:
        alias = (char*) "brxc";
        break;
      case BASIC_FUNCTIONALS::LDAERFX:
        alias = (char*) "ldaerfx";
        break;
      case BASIC_FUNCTIONALS::LDAERFC:
        alias = (char*) "ldaerfc";
        break;
      case BASIC_FUNCTIONALS::LDAERFC_JT:
        alias = (char*) "ldaerfc_jt";
        break;
      case BASIC_FUNCTIONALS::LYPC:
        alias = (char*) "lypc";
        break;
      case BASIC_FUNCTIONALS::OPTX:
        alias = (char*) "optx";
        break;
      case BASIC_FUNCTIONALS::OPTXCORR:
        alias = (char*) "optxcorr";
        break;
      case BASIC_FUNCTIONALS::REVPBEX:
        alias = (char*) "revpbex";
        break;
      case BASIC_FUNCTIONALS::RPBEX:
        alias = (char*) "rpbex";
        break;
      case BASIC_FUNCTIONALS::SPBEC:
        alias = (char*) "spbec";
        break;
      case BASIC_FUNCTIONALS::VWN_PBEC:
        alias = (char*) "vwn_pbec";
        break;
      case BASIC_FUNCTIONALS::KTX:
        alias = (char*) "ktx";
        break;
      case BASIC_FUNCTIONALS::TFK:
        alias = (char*) "tfk";
        break;
      case BASIC_FUNCTIONALS::TW:
        alias = (char*) "tw";
        break;
      case BASIC_FUNCTIONALS::PW91X:
        alias = (char*) "pw91x";
        break;
      case BASIC_FUNCTIONALS::PW91K:
        alias = (char*) "pw91k";
        break;
      case BASIC_FUNCTIONALS::PW91C:
        alias = (char*) "pw91c";
        break;
      case BASIC_FUNCTIONALS::M05X:
        alias = (char*) "m05x";
        break;
      case BASIC_FUNCTIONALS::M05X2X:
        alias = (char*) "m05x2x";
        break;
      case BASIC_FUNCTIONALS::M06X:
        alias = (char*) "m06x";
        break;
      case BASIC_FUNCTIONALS::M06X2X:
        alias = (char*) "m06x2x";
        break;
      case BASIC_FUNCTIONALS::M06LX:
        alias = (char*) "m06lx";
        break;
      case BASIC_FUNCTIONALS::M06HFX:
        alias = (char*) "m06hfx";
        break;
      case BASIC_FUNCTIONALS::M05X2C:
        alias = (char*) "m05x2c";
        break;
      case BASIC_FUNCTIONALS::M05C:
        alias = (char*) "m05c";
        break;
      case BASIC_FUNCTIONALS::M06C:
        alias = (char*) "m06c";
        break;
      case BASIC_FUNCTIONALS::M06HFC:
        alias = (char*) "m06hfc";
        break;
      case BASIC_FUNCTIONALS::M06LC:
        alias = (char*) "m06lc";
        break;
      case BASIC_FUNCTIONALS::M06X2C:
        alias = (char*) "m06x2c";
        break;
      case BASIC_FUNCTIONALS::TPSSC:
        alias = (char*) "tpssc";
        break;
      case BASIC_FUNCTIONALS::TPSSX:
        alias = (char*) "tpssx";
        break;
      case BASIC_FUNCTIONALS::REVTPSSC:
        alias = (char*) "revtpssc";
        break;
      case BASIC_FUNCTIONALS::REVTPSSX:
        alias = (char*) "revtpssx";
        break;
      case BASIC_FUNCTIONALS::PZ81C:
        alias = (char*) "pz81c";
        break;
      case BASIC_FUNCTIONALS::P86C:
        alias = (char*) "p86c";
        break;
      case BASIC_FUNCTIONALS::P86CORRC:
        alias = (char*) "p86corrc";
        break;
      case BASIC_FUNCTIONALS::BTK:
        alias = (char*) "btk";
        break;
      case BASIC_FUNCTIONALS::VWK:
        alias = (char*) "vwk";
        break;
      case BASIC_FUNCTIONALS::B97X:
        alias = (char*) "b97x";
        break;
      case BASIC_FUNCTIONALS::B97C:
        alias = (char*) "b97c";
        break;
      case BASIC_FUNCTIONALS::B97_1X:
        alias = (char*) "b97_1x";
        break;
      case BASIC_FUNCTIONALS::B97_1C:
        alias = (char*) "b97_1c";
        break;
      case BASIC_FUNCTIONALS::B97_2X:
        alias = (char*) "b97_2x";
        break;
      case BASIC_FUNCTIONALS::B97_2C:
        alias = (char*) "b97_2c";
        break;
      case BASIC_FUNCTIONALS::CSC:
        alias = (char*) "csc";
        break;
      case BASIC_FUNCTIONALS::APBEC:
        alias = (char*) "apbec";
        break;
      case BASIC_FUNCTIONALS::APBEX:
        alias = (char*) "apbex";
        break;
      case BASIC_FUNCTIONALS::ZVPBESOLC:
        alias = (char*) "zvpbesolc";
        break;
      case BASIC_FUNCTIONALS::BLOCX:
        alias = (char*) "blocx";
        break;
      case BASIC_FUNCTIONALS::PBEINTC:
        alias = (char*) "pbeintc";
        break;
      case BASIC_FUNCTIONALS::PBEINTX:
         alias = (char*) "pbeintx";
         break;
      case BASIC_FUNCTIONALS::PBELOCC:
         alias = (char*) "pbelocc";
         break;
      case BASIC_FUNCTIONALS::PBESOLX:
         alias = (char*) "pbesolx";
         break;
      case BASIC_FUNCTIONALS::TPSSLOCC:
    	  alias = (char*) "tpsslocc";
    	  break;
      case BASIC_FUNCTIONALS::ZVPBEINTC:
    	  alias = (char*) "zvpbeintc";
    	  break;
      case BASIC_FUNCTIONALS::LLP91K:
    	  alias = (char*) "llp91k";
    	  break;
      case BASIC_FUNCTIONALS::PBE2K:
    	  alias = (char*) "PBE2";
    	  break;
      case BASIC_FUNCTIONALS::PBE2KS:
    	  alias = (char*) "PBE2S";
    	  break;
      case BASIC_FUNCTIONALS::LLP91KS:
    	  alias = (char*) "llp91ks";
    	  break;
      case BASIC_FUNCTIONALS::PBE3K:
    	  alias = (char*) "PBE3";
    	  break;
      case BASIC_FUNCTIONALS::PBE4K:
    	  alias = (char*) "PBE4";
    	  break;
      case BASIC_FUNCTIONALS::E2000K:
    	  alias = (char*) "E00";
    	  break;
      default:
        throw SerenityError("Unknown functional to XCFun.cpp");
        break;
    }
  return alias;
}

template<Options::SCF_MODES T>
FunctionalData<T> XCFun<T>::calcData(
    FUNCTIONAL_DATA_TYPE type,
    const Functional functional,
    const std::shared_ptr<DensityOnGridController<T> > densityOnGridController,
    unsigned int order) {
  //Check input
  assert(order <= 2 && "Only implemented up to 2nd order");
  //Build functional from basic functionals
  auto func = getFunctional(functional);
  //GGA?
  const bool gga = (xc_is_gga(func) != 0) ? true : false;
  //Get density
  auto& density = densityOnGridController->getDensityOnGrid();
  //If gga, get gradient of density
  std::shared_ptr<Gradient<DensityOnGrid<T> > > gradient = nullptr;
  //If gga and potential is requested, get also hessian
  std::shared_ptr<Hessian<DensityOnGrid<T> > > hessian = nullptr;
  if (gga) {
    assert(densityOnGridController->getHighestDerivative() >= 1);
    gradient = std::make_shared<Gradient<DensityOnGrid<T> > >(densityOnGridController->getDensityGradientOnGrid());
    if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
      assert(densityOnGridController->getHighestDerivative() >= 2);
      hessian = std::make_shared<Hessian<DensityOnGrid<T> > >(densityOnGridController->getDensityHessianOnGrid());
    }
  }
  Timings::takeTime("Tech. - XCFun Functional Eval.");
  //Setup xcFun
  int iErr;
  xc_vars xcVars = xc_vars::XC_N;
  xc_mode xcMode = (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) ? XC_POTENTIAL : XC_PARTIAL_DERIVATIVES;
  if (!gga){
    //If LDA:
    if (T == Options::SCF_MODES::RESTRICTED) {
      xcVars = xc_vars::XC_N;
    } else {
      xcVars = xc_vars::XC_A_B;
    }
    iErr = xc_eval_setup(func,xcVars,XC_PARTIAL_DERIVATIVES,order);
  } else {
    //If GGA:
    if (T == Options::SCF_MODES::RESTRICTED) {
      //If restricted
      if (type == FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS) {
        xcVars = xc_vars::XC_N_GNN;
      } else if (type == FUNCTIONAL_DATA_TYPE::GRADIENTS) {
        xcVars = xc_vars::XC_N_NX_NY_NZ;
      } else if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
        order = 1;
        xcVars = xc_vars::XC_N_2ND_TAYLOR;
      } else {
        assert(false);
      }
    } else {
      //If unrestricted
      if (type == FUNCTIONAL_DATA_TYPE::GRADIENT_INVARIANTS) {
        xcVars = xc_vars::XC_A_B_GAA_GAB_GBB;
      } else if (type == FUNCTIONAL_DATA_TYPE::GRADIENTS) {
        xcVars = xc_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ;
      }  else if (type == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
        order = 1;
        xcVars = xc_vars::XC_A_B_2ND_TAYLOR;
      } else {
        assert(false && "Problem with functional data type");
      }
    }
    iErr = xc_eval_setup(func,xcVars,xcMode,order);
  }
  if (iErr != 0) throw SerenityError("Failed to set vars, mode and order for xcFun.");
  //In and output dimension
  int nOut = xc_output_length(func);
  int nIn = xc_input_length(func);

  //Number of grid points and blocks
  const unsigned int nPoints = density.getNGridPoints();
  const unsigned int nBlocks = (unsigned int)ceil((double)nPoints/_maxBlockSize);
  //Build FunctionaData object
  FunctionalData<T> funcData(
      order,
      type,
      functional,
      densityOnGridController->getGridController());
  //Loop over blocks of grid points
#pragma omp parallel for schedule(dynamic)
  for (unsigned int iBlock = 0; iBlock < nBlocks; ++iBlock) {
    //first index of block
    const unsigned int firstIndex = iBlock * _maxBlockSize;
    //size of this block
    const unsigned int blockSize= determineBlockSize(iBlock, nPoints, nBlocks);
    bool skip = true;
    for_spin(density){
     skip *= (density_spin.segment(firstIndex,blockSize).array().abs().sum()<blockSize*1e-12);
    };
    if(skip) continue;
    //Input for this block
    Eigen::MatrixXd input(nIn,blockSize);
    prepareInput(firstIndex, xcVars,density, gradient, hessian, input);
    //Output object for this block
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(nOut,blockSize);
    double* ptr = input.data();
    double* outpt = output.data();
    for (unsigned int p=0;p<blockSize;p++){
        xc_eval(func,ptr+p*nIn,outpt+p*nOut);
    }
    //parse data for this block
    parseOutput(order, xcVars, firstIndex, output, funcData);
  }/* Loop over blocks */
  xc_free_functional(func);
  funcData.energy = calcEnergy(funcData.epuv,densityOnGridController->getGridController()->getWeights());
  Timings::timeTaken("Tech. - XCFun Functional Eval.");
  return funcData;
}

template<Options::SCF_MODES T>
unsigned int XCFun<T>::determineBlockSize(
    unsigned int blockIndex, unsigned int nPoints, unsigned int nBlocks) {
  unsigned int blockSize;
    if (blockIndex == nBlocks-1) {
      blockSize = nPoints % _maxBlockSize;
      if (blockSize == 0) blockSize = _maxBlockSize;
    } else {
      blockSize = _maxBlockSize;
    }
    return blockSize;
}

template<>Eigen::MatrixXd XCFun<Options::SCF_MODES::RESTRICTED>::calculateSigma(
    const Gradient<DensityOnGrid<Options::SCF_MODES::RESTRICTED> >& gradient,
    const unsigned int& iGridStart,
    const unsigned int& blocksize) {
  Eigen::MatrixXd sigma(1,blocksize);
  sigma.row(0).array()=
      gradient.x.segment(iGridStart,blocksize).array() * gradient.x.segment(iGridStart,blocksize).array() +
      gradient.y.segment(iGridStart,blocksize).array() * gradient.y.segment(iGridStart,blocksize).array() +
      gradient.z.segment(iGridStart,blocksize).array() * gradient.z.segment(iGridStart,blocksize).array();
  return sigma;
}

template<>Eigen::MatrixXd XCFun<Options::SCF_MODES::UNRESTRICTED>::calculateSigma(
    const Gradient<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED> >& gradient,
    const unsigned int& iGridStart,
    const unsigned int& blocksize) {
  Eigen::MatrixXd sigma(3,blocksize);
  sigma.row(0).array()=
      gradient.x.alpha.segment(iGridStart,blocksize).array() * gradient.x.alpha.segment(iGridStart,blocksize).array() +
      gradient.y.alpha.segment(iGridStart,blocksize).array() * gradient.y.alpha.segment(iGridStart,blocksize).array() +
      gradient.z.alpha.segment(iGridStart,blocksize).array() * gradient.z.alpha.segment(iGridStart,blocksize).array();
  sigma.row(1).array()=
      gradient.x.alpha.segment(iGridStart,blocksize).array() * gradient.x.beta.segment(iGridStart,blocksize).array() +
      gradient.y.alpha.segment(iGridStart,blocksize).array() * gradient.y.beta.segment(iGridStart,blocksize).array() +
      gradient.z.alpha.segment(iGridStart,blocksize).array() * gradient.z.beta.segment(iGridStart,blocksize).array();
  sigma.row(2).array()=
      gradient.x.beta.segment(iGridStart,blocksize).array() * gradient.x.beta.segment(iGridStart,blocksize).array() +
      gradient.y.beta.segment(iGridStart,blocksize).array() * gradient.y.beta.segment(iGridStart,blocksize).array() +
      gradient.z.beta.segment(iGridStart,blocksize).array() * gradient.z.beta.segment(iGridStart,blocksize).array();
  return sigma;
}

template<>void XCFun<Options::SCF_MODES::RESTRICTED>::prepareInput(
    const unsigned int iGridStart,
    const xc_vars xcVars,
    const DensityOnGrid<Options::SCF_MODES::RESTRICTED>& density,
    const std::shared_ptr<Gradient<DensityOnGrid<Options::SCF_MODES::RESTRICTED> > > gradient,
    const std::shared_ptr<Hessian<DensityOnGrid<Options::SCF_MODES::RESTRICTED> > > hessian,
    Eigen::MatrixXd& input){
  unsigned int n = input.cols();
  if (xcVars == xc_vars::XC_N) {
    input.row(0) = density.segment(iGridStart,n);
  } else if (xcVars == xc_vars::XC_N_GNN) {
    input.row(0) = density.segment(iGridStart,n);
    input.row(1) = calculateSigma(*gradient,iGridStart,n);
  } else if (xcVars == xc_vars::XC_N_2ND_TAYLOR) {
    input.row(0) = density.segment(iGridStart,n);
    input.row(1) = gradient->x.segment(iGridStart,n);
    input.row(2) = gradient->y.segment(iGridStart,n);
    input.row(3) = gradient->z.segment(iGridStart,n);
    input.row(4) = hessian->xx.segment(iGridStart,n);
    input.row(5) = hessian->xy.segment(iGridStart,n);
    input.row(6) = hessian->xz.segment(iGridStart,n);
    input.row(7) = hessian->yy.segment(iGridStart,n);
    input.row(8) = hessian->yz.segment(iGridStart,n);
    input.row(9) = hessian->zz.segment(iGridStart,n);
  } else if (xcVars == xc_vars::XC_N_NX_NY_NZ) {
    input.row(0) = density.segment(iGridStart,n);
    input.row(1) = gradient->x.segment(iGridStart,n);
    input.row(2) = gradient->y.segment(iGridStart,n);
    input.row(3) = gradient->z.segment(iGridStart,n);
  } else {
    throw SerenityError("input for requested xc_vars not implemented");
  }
}

template<> void XCFun<Options::SCF_MODES::RESTRICTED>::parseOutput(
    const unsigned int order,
    const xc_vars xcVars,
    const unsigned int iGridStart,
    Eigen::MatrixXd& output,
    FunctionalData<Options::SCF_MODES::RESTRICTED>& funcData){
  unsigned int n = output.cols();
  if (funcData.getType() == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
    (*funcData.epuv).segment(iGridStart,n) = output.row(0);
    (*funcData.potential).segment(iGridStart,n) = output.row(1);
  } else {
    if (xcVars == xc_vars::XC_N) {
      (*funcData.epuv).segment(iGridStart,n) = output.row(0);
      if (order >= 1) (*funcData.dFdRho).segment(iGridStart,n) = output.row(1);
      if (order >= 2) (*funcData.d2FdRho2).segment(iGridStart,n) = output.row(2);
    } else if (xcVars == xc_vars::XC_N_GNN) {
      (*funcData.epuv).segment(iGridStart,n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).segment(iGridStart,n) = output.row(1);
        (*funcData.dFdSigma).segment(iGridStart,n) = output.row(2);
      }
      if (order >= 2) {
        (*funcData.d2FdRho2).segment(iGridStart,n) = output.row(3);
        (*funcData.d2FdRhodSigma).segment(iGridStart,n) = output.row(4);
        (*funcData.d2FdSigma2).segment(iGridStart,n) = output.row(5);
      }
    } else if (xcVars == xc_vars::XC_N_NX_NY_NZ) {
      (*funcData.epuv).segment(iGridStart,n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).segment(iGridStart,n) = output.row(1);
        (*funcData.dFdGradRho).x.segment(iGridStart,n) = output.row(2);
        (*funcData.dFdGradRho).y.segment(iGridStart,n) = output.row(3);
        (*funcData.dFdGradRho).z.segment(iGridStart,n) = output.row(4);
      }
      if (order >= 2) {
        throw SerenityError("GRADIENTS type only implemented up to first order");
      }
    } else {
      throw SerenityError("output for requested xc_vars not implemented");
    }
  }
}

template<>void XCFun<Options::SCF_MODES::UNRESTRICTED>::prepareInput(
    const unsigned int iGridStart,
    const xc_vars xcVars,
    const DensityOnGrid<Options::SCF_MODES::UNRESTRICTED>& density,
    const std::shared_ptr<Gradient<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED> > > gradient,
    const std::shared_ptr<Hessian<DensityOnGrid<Options::SCF_MODES::UNRESTRICTED> > > hessian,
    Eigen::MatrixXd& input){
  unsigned int n = input.cols();
  if (xcVars == xc_vars::XC_A_B) {
    input.row(0) = density.alpha.segment(iGridStart,n);
    input.row(1) = density.beta.segment(iGridStart,n);
  } else if (xcVars == xc_vars::XC_A_B_GAA_GAB_GBB) {
    input.row(0) = density.alpha.segment(iGridStart,n);
    input.row(1) = density.beta.segment(iGridStart,n);
    input.bottomRows(3) = calculateSigma(*gradient,iGridStart,n);
  } else if (xcVars == xc_vars::XC_A_B_2ND_TAYLOR) {
    input.row(0) = density.alpha.segment(iGridStart,n);
    input.row(1) = gradient->x.alpha.segment(iGridStart,n);
    input.row(2) = gradient->y.alpha.segment(iGridStart,n);
    input.row(3) = gradient->z.alpha.segment(iGridStart,n);
    input.row(4) = hessian->xx.alpha.segment(iGridStart,n);
    input.row(5) = hessian->xy.alpha.segment(iGridStart,n);
    input.row(6) = hessian->xz.alpha.segment(iGridStart,n);
    input.row(7) = hessian->yy.alpha.segment(iGridStart,n);
    input.row(8) = hessian->yz.alpha.segment(iGridStart,n);
    input.row(9) = hessian->zz.alpha.segment(iGridStart,n);
    input.row(10) = density.beta.segment(iGridStart,n);
    input.row(11) = gradient->x.beta.segment(iGridStart,n);
    input.row(12) = gradient->y.beta.segment(iGridStart,n);
    input.row(13) = gradient->z.beta.segment(iGridStart,n);
    input.row(14) = hessian->xx.beta.segment(iGridStart,n);
    input.row(15) = hessian->xy.beta.segment(iGridStart,n);
    input.row(16) = hessian->xz.beta.segment(iGridStart,n);
    input.row(17) = hessian->yy.beta.segment(iGridStart,n);
    input.row(18) = hessian->yz.beta.segment(iGridStart,n);
    input.row(19) = hessian->zz.beta.segment(iGridStart,n);
  } else if (xcVars == xc_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ) {
    input.row(0) = density.alpha.segment(iGridStart,n);
    input.row(1) = density.beta.segment(iGridStart,n);
    input.row(2) = gradient->x.alpha.segment(iGridStart,n);
    input.row(3) = gradient->y.alpha.segment(iGridStart,n);
    input.row(4) = gradient->z.alpha.segment(iGridStart,n);
    input.row(5) = gradient->x.beta.segment(iGridStart,n);
    input.row(6) = gradient->y.beta.segment(iGridStart,n);
    input.row(7) = gradient->z.beta.segment(iGridStart,n);
  } else {
    throw SerenityError("input for requested xc_vars not implemented");
  }
}

template<> void XCFun<Options::SCF_MODES::UNRESTRICTED>::parseOutput(
    const unsigned int order,
    const xc_vars xcVars,
    const unsigned int iGridStart,
    Eigen::MatrixXd& output,
    FunctionalData<Options::SCF_MODES::UNRESTRICTED>& funcData){
  unsigned int n = output.cols();
  if (funcData.getType() == FUNCTIONAL_DATA_TYPE::POTENTIAL) {
    (*funcData.epuv).segment(iGridStart,n) = output.row(0);
    (*funcData.potential).alpha.segment(iGridStart,n) = output.row(1);
    (*funcData.potential).beta.segment(iGridStart,n) = output.row(2);

  } else {
    if (xcVars == xc_vars::XC_A_B) {
      (*funcData.epuv).segment(iGridStart,n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).alpha.segment(iGridStart,n) = output.row(1);
        (*funcData.dFdRho).beta.segment(iGridStart,n) = output.row(2);
      }
      if (order >= 2) {
        (*funcData.d2FdRho2).aa.segment(iGridStart,n) = output.row(3);
        (*funcData.d2FdRho2).ab.segment(iGridStart,n) = output.row(4);
        (*funcData.d2FdRho2).bb.segment(iGridStart,n) = output.row(5);
      }

    } else if (xcVars == xc_vars::XC_A_B_GAA_GAB_GBB) {
      (*funcData.epuv).segment(iGridStart,n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).alpha.segment(iGridStart,n) = output.row(1);
        (*funcData.dFdRho).beta.segment(iGridStart,n) = output.row(2);
        (*funcData.dFdSigma).aa.segment(iGridStart,n) = output.row(3);
        (*funcData.dFdSigma).ab.segment(iGridStart,n) = output.row(4);
        (*funcData.dFdSigma).bb.segment(iGridStart,n) = output.row(5);
      }
      if (order >= 2) {
        (*funcData.d2FdRho2).aa.segment(iGridStart,n) = output.row(6);
        (*funcData.d2FdRho2).ab.segment(iGridStart,n) = output.row(7);
        (*funcData.d2FdRho2).bb.segment(iGridStart,n) = output.row(11);
        (*funcData.d2FdRhodSigma).aaa.segment(iGridStart,n) = output.row(8);
        (*funcData.d2FdRhodSigma).aab.segment(iGridStart,n) = output.row(9);
        (*funcData.d2FdRhodSigma).abb.segment(iGridStart,n) = output.row(10);
        (*funcData.d2FdRhodSigma).baa.segment(iGridStart,n) = output.row(12);
        (*funcData.d2FdRhodSigma).bab.segment(iGridStart,n) = output.row(13);
        (*funcData.d2FdRhodSigma).bbb.segment(iGridStart,n) = output.row(14);
        (*funcData.d2FdSigma2).aaaa.segment(iGridStart,n) = output.row(15);
        (*funcData.d2FdSigma2).aaab.segment(iGridStart,n) = output.row(16);
        (*funcData.d2FdSigma2).aabb.segment(iGridStart,n) = output.row(17);
        (*funcData.d2FdSigma2).abab.segment(iGridStart,n) = output.row(18);
        (*funcData.d2FdSigma2).abbb.segment(iGridStart,n) = output.row(19);
        (*funcData.d2FdSigma2).bbbb.segment(iGridStart,n) = output.row(20);
      }

    } else if (xcVars == xc_vars::XC_A_B_AX_AY_AZ_BX_BY_BZ) {
      (*funcData.epuv).segment(iGridStart,n) = output.row(0);
      if (order >= 1) {
        (*funcData.dFdRho).alpha.segment(iGridStart,n) = output.row(1);
        (*funcData.dFdRho).beta.segment(iGridStart,n) = output.row(2);
        (*funcData.dFdGradRho).x.alpha.segment(iGridStart,n) = output.row(3);
        (*funcData.dFdGradRho).y.alpha.segment(iGridStart,n) = output.row(4);
        (*funcData.dFdGradRho).z.alpha.segment(iGridStart,n) = output.row(5);
        (*funcData.dFdGradRho).x.beta.segment(iGridStart,n) = output.row(6);
        (*funcData.dFdGradRho).y.beta.segment(iGridStart,n) = output.row(7);
        (*funcData.dFdGradRho).z.beta.segment(iGridStart,n) = output.row(8);
      }
      if (order >= 2) {
        throw SerenityError("GRADIENTS type only implemented up to first order");
      }

    } else {
      throw SerenityError("output for requested xc_vars not implemented");
    }
  }
}

template<Options::SCF_MODES T> xc_functional XCFun<T>::getFunctional(Functional functional) {
  //Status integer for xcFun
  int iErr;
  //Build functional from set of basic functionals
  auto functionals = functional.getBasicFunctionals();
  auto mixingFactors = functional.getMixingFactors();
  xc_functional func = xc_new_functional();
  for (unsigned int iFunc = 0; iFunc < functionals.size(); ++iFunc) {
    iErr = xc_set(func,getAlias(functionals[iFunc]),mixingFactors[iFunc]);
    if (iErr != 0) {
      std::cout << "\n Functional " << getAlias(functionals[iFunc]) << " unknown to xcFun. Check if alias is set correctly.\n";
      throw SerenityError("Error");
    }
  }
  //Provide additional parameters for RS-Hybrid
  if (functional.isRSHybrid()) {
    iErr = xc_set(func,"cam_alpha", functional.getHfExchangeRatio());
    if (iErr != 0) {
      std::cout << "\n XCFun: Failed to set cam_alpha" << std::endl;
      throw SerenityError("Error");
    }
    iErr = xc_set(func,"cam_beta", functional.getLRExchangeRatio());
    if (iErr != 0) {
      std::cout << "\n XCFun: Failed to set cam_beta" << std::endl;
      throw SerenityError("Error");
    }
    iErr = xc_set(func,"rangesep_mu", functional.getRangeSeparationParameter());
    if (iErr != 0) {
      std::cout << "\n XCFun: Failed to set rangesep_mu" << std::endl;
      throw SerenityError("Error");
    }
  }
  return func;
}

template<Options::SCF_MODES T> double XCFun<T>::calcEnergy(
    std::shared_ptr<GridData<RESTRICTED> > epuv,
    const Eigen::VectorXd& weights) {
  unsigned int nPoints = weights.size();
  unsigned int nBlocks = omp_get_max_threads();
  double energy = 0.0;
#pragma omp parallel for schedule (static,1) reduction(+:energy)
  for (unsigned int iBlock = 0; iBlock<nBlocks;iBlock++){
    unsigned int n = (unsigned int)(nPoints/nBlocks);
    const unsigned int start = iBlock*n;
    if(iBlock==nBlocks-1) n += nPoints%nBlocks;
    energy += (*epuv).segment(start,n).dot(weights.segment(start,n));
  }
  return energy;
}

template class XCFun<Options::SCF_MODES::RESTRICTED>;
template class XCFun<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
