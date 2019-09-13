/**
 * @file   LRSCFTask.cpp
 *
 * @date   Aug 17, 2016
 * @author M. Boeckers, N. Niemeyer, J. Toelle
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
#include "tasks/LRSCFTask.h"

/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"
#include "postHF/LRSCF/Analysis/DipoleIntegrals.h"
#include "postHF/LRSCF/Analysis/ExcitationSpectrum.h"
#include "postHF/LRSCF/Tools/LRSCFSetup.h"

#include "postHF/LRSCF/Kernel/Kernel.h"

#include "postHF/LRSCF/SigmaVectors/CoulombSigmaVector.h"
#include "postHF/LRSCF/SigmaVectors/KSigmaVector.h"
#include "postHF/LRSCF/SigmaVectors/DeltaESigmaVector.h"
#include "postHF/LRSCF/SigmaVectors/KernelSigmaVector.h"
#include "postHF/LRSCF/SigmaVectors/EOSigmaVector.h"

#include "postHF/LRSCF/Tools/EigenvalueSolver.h"

#include "data/ElectronicStructure.h"
#include "potentials/bundles/PotentialBundle.h"

namespace Serenity {

template<Options::SCF_MODES SCFMode>
LRSCFTask<SCFMode>::LRSCFTask(
    const std::vector<std::shared_ptr<SystemController> >& activeSystems,
    const std::vector<std::shared_ptr<SystemController> >& passiveSystems) :
      _act(activeSystems),
      _env(passiveSystems){
}

template<Options::SCF_MODES SCFMode>
void LRSCFTask<SCFMode>::run() {

  printSectionTitle("LRSCF");

  //Set type of calculation
  Options::LRSCF_TYPE type = Options::LRSCF_TYPE::ISOLATED;
  if (_act.size() + _env.size() == 1) {
    type = Options::LRSCF_TYPE::ISOLATED;
  } else if (_act.size() == 1 && _env.size() >=1) {
    type = Options::LRSCF_TYPE::UNCOUPLED;
  } else if (_act.size() > 1) {
    type = Options::LRSCF_TYPE::COUPLED;
  } else {
    assert(false);
  }

  //Creates lrscf controller contains all data needed for the LRSCF calculation
  for (const auto& system : _act){
    _lrscf.push_back(std::make_shared<LRSCFController<SCFMode> >(system,settings));
  }
  //if special coupling pattern is requested
  if(!settings.couplingPattern.empty()){
    const unsigned int numberOfRows = std::sqrt(settings.couplingPattern.size());
    Eigen::MatrixXi couplingPatternMatrix(numberOfRows,numberOfRows);
    couplingPatternMatrix.setZero();   
    for (unsigned int i = 0; i < numberOfRows; i++){
      if(i < _act.size()){
        _referenceLoadingType.push_back(Options::LRSCF_TYPE::COUPLED);
      }else{
        _referenceLoadingType.push_back(settings.loadType);
      }
      for (unsigned int j = 0; j < numberOfRows; j ++){
        couplingPatternMatrix(i,j) = settings.couplingPattern[(i * numberOfRows) + j];
      }
    }
    _couplingPatternMatrix = couplingPatternMatrix;
    printf("  ----------------------------------------------------------------\n");
    printf("                 Special coupling pattern is used!                \n");
    printf("  ----------------------------------------------------------------\n");
    std::cout<<_couplingPatternMatrix<<std::endl;
    //Adds additional lrscfController if one system is used coupled and uncoupled
    if((unsigned int) _couplingPatternMatrix.rows() > _act.size()){
      for(unsigned int iRow = _act.size(); iRow < _couplingPatternMatrix.rows(); iRow ++){
        //the systemcontroller which occurs two times
        unsigned int index = _couplingPatternMatrix(iRow,iRow) - 1;
        _lrscf.push_back(std::make_shared<LRSCFController<SCFMode> >(_act[index],settings));
      }
    }
  //Regular coupled calculation
  }else if(type == Options::LRSCF_TYPE::COUPLED && settings.couplingPattern.empty()){
    Eigen::MatrixXi couplingPatternMatrix(_act.size(),_act.size());
    couplingPatternMatrix.setZero();
    for(unsigned int i = 0; i < _act.size(); i++){
      couplingPatternMatrix(i,i) = i+1;
      _referenceLoadingType.push_back(settings.loadType);
    }
    _couplingPatternMatrix = couplingPatternMatrix;
  }
  
  //This is for off-diagonal Lagrange Mulitplier when performing supersystem TDDFT with local orbitals
  if(settings.localMO){
    if(_lrscf.size() != 1) throw SerenityError("LMO-TDDFT only available for supersystem calculations.");
    
    auto es = _act[0]->getElectronicStructure<SCFMode>();
    auto energyComponentController =  es->getEnergyComponentController();
    auto potentials = _act[0]->getPotentials<SCFMode,Options::ELECTRONIC_STRUCTURE_THEORIES::DFT>(
          Options::GRID_PURPOSES::DEFAULT);
    auto orbitalController = es->getMolecularOrbitals();
    auto coeff = _lrscf[0]->getCoefficients();
    auto eigenvalues = _lrscf[0] -> getEigenvalues();
    auto fock(potentials->getFockMatrix(es->getDensityMatrix(),energyComponentController));
    
    for_spin(coeff,fock){
      auto newF = coeff_spin.transpose() * fock_spin * coeff_spin;
      printf("  ---------------------------------------------------------\n");
      printf("          Local MO Transformed Fock Matrix used!           \n");
      printf("  ---------------------------------------------------------\n");
      fock_spin = newF;
    };
    _lrscf[0]->setFockNonCanon(std::make_shared<MatrixInBasis<SCFMode> >(fock));
  }

  //Setup LRSCFController
  LRSCFSetup<SCFMode>::setupLRSCFController(settings,_couplingPatternMatrix,_act,_env,_lrscf,type);

  //Dimension of LRSCF problem
  unsigned int nDimension = 0;
  for (unsigned int I = 0; I < _lrscf.size(); ++I) {
    auto nOccupied = _lrscf[I]->getNOccupied();
    auto nVirtual = _lrscf[I]->getNVirtual();
    for_spin(nOccupied,nVirtual) {nDimension += nOccupied_spin * nVirtual_spin;};
  }

  /*
   * Adjust number of roots to be determined if the user's input
   * exceeds the actual dimension of the eigenvalue problem
   */
  if(settings.nEigen > nDimension) settings.nEigen = nDimension;

  //Sanity check for gauge-origin input
  if(settings.gaugeOrigin.size() != 3) throw SerenityError("Gauge-origin expects three cartesian coordinates!");

  //Information about the TDDFT calculation
  LRSCFSetup<SCFMode>::printInfo(_lrscf,settings,_env,type);

  //Prepare eigenvectors and eigenvalues
  std::shared_ptr<std::vector<Eigen::MatrixXd> > eigenvectors = nullptr;
  Eigen::VectorXd eigenvalues;

  //Prepare solutions vector in case of a FDEc calculation
  if (type == Options::LRSCF_TYPE::COUPLED){
    eigenvectors = LRSCFSetup<SCFMode>::setupFDEcTransformation(
      settings,_couplingPatternMatrix,_referenceLoadingType,_lrscf,_act,nDimension);
    //Setting the numer of the eigenvalues to be determined in the coupled step
    settings.nEigen = (*eigenvectors)[0].cols();
  }

  //Prepare vector with orbital energy differences
  Eigen::VectorXd diagonal(nDimension);
  unsigned int iStart = 0; 
  for (unsigned int I = 0; I < _lrscf.size(); ++I) {
    auto nOccupied = _lrscf[I]->getNOccupied();
    auto nVirtual = _lrscf[I]->getNVirtual();
    auto orbitalEnergies = _lrscf[I]->getEigenvalues();
    for_spin(nOccupied,nVirtual,orbitalEnergies) {
      for (unsigned int ia = 0; ia < nOccupied_spin*nVirtual_spin; ++ia) {
        unsigned int i = floor(ia/nVirtual_spin);
        unsigned int a = nOccupied_spin + ia - i*nVirtual_spin;
        diagonal(iStart+ia) = orbitalEnergies_spin(a) - orbitalEnergies_spin(i);
      }
      iStart+=nOccupied_spin*nVirtual_spin;
    };
  }

  //Calculate kernel
  bool fxcSigma = false;
  std::shared_ptr<Kernel<SCFMode> > kernel = nullptr;
  //Kernel contributions are needed if DFT is being used in any of the active systems.
  //For environment systems it is assumed that they use DFT.
  for (const auto& sys : _act){
    if(sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT) fxcSigma = true;
  }
  if(_env.size() >= 1) fxcSigma = true;

  if(fxcSigma){
    printf("\n Prepare kernel ..\n\n");
    kernel = std::make_shared<Kernel<SCFMode> >(
      _act,_env,settings.superSystemGrid,
      settings.embedding.naddKinFunc,settings.embedding.naddXCFunc,settings.func);
    printf(" .. done.\n\n");
  }/* Calculate kernel */

  /*
   * Lambda functions
   * 
   * These can be conveniently passed to an iterative eigenvalue solver which
   * can use them to compute response matrix products with guess vectors
   */

  /*
   * Some booleans to determine which sigma vector contributions are necessary
   */
  bool kSigma = false;
  bool eoSigma = false;

  //Exchange sigma to be used
  for(const auto& sys : _act){
    if(sys->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF){
      kSigma = true;
    }else{
      auto func = FunctionalClassResolver::resolveFunctional(sys->getSettings().dft.functional);
      if(func.isHybrid() || func.getLRExchangeRatio() > 0.0) kSigma = true;
    }
  }

  if(settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::LEVELSHIFT){
    eoSigma = true;
  }
  if(settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::HUZINAGA){
    eoSigma = true;
  }
  if(settings.embedding.embeddingMode == Options::KIN_EMBEDDING_MODES::HOFFMANN){
    eoSigma = true;
  }

  double factor = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 4.0 : 2.0;
  double factorEO = (SCFMode == Options::SCF_MODES::RESTRICTED) ? 2.0 : 1.0;

  //Pointer for sigma vectors    
  std::shared_ptr<SigmaVector<SCFMode> > D = nullptr, J = nullptr,
                                         F = nullptr, K = nullptr, EO = nullptr;
  //Convenience typedef
  typedef std::function<std::unique_ptr<std::vector<Eigen::MatrixXd> >(
    std::vector<Eigen::MatrixXd>& guessVectors)> ResponseLambda;

  //Lambda function for OJJ eigenvalue solver
  //Note: (A+B)*b stored in [0,2,4,..], (A-B)*b stored in [1,3,5,..]
  ResponseLambda MKSigma = [&] (std::vector<Eigen::MatrixXd>& guessVectors) {
    
    //Sanity check
    assert(guessVectors.size() % 2 == 0);

    //Build vector with signs between A and B for each set
    std::vector<int> xc(guessVectors.size());
    //Build vector that only contains b for (A+B)*b
    std::vector<Eigen::MatrixXd> guessAPB(guessVectors.size() / 2);
    for(unsigned iSet = 0; iSet < guessVectors.size(); ++iSet){
      //(A+B)
      if(iSet % 2 == 0){
        xc[iSet] = 1;
        guessAPB[iSet/2] = guessVectors[iSet];
      //(A-B)
      }else{
        xc[iSet] = -1;
      }
    }
    
    //Orbital-energy differences
    D = std::make_shared<DeltaESigmaVector<SCFMode> >
      (_lrscf,guessVectors,settings.pseudoDensityThreshold);
    //Coulomb (only needed for first set of guess vectors (A+B))
    J = std::make_shared<CoulombSigmaVector<SCFMode> >
      (_lrscf,guessAPB,settings.pseudoDensityThreshold);
    //Exchange
    if(kSigma) K = std::make_shared<KSigmaVector<SCFMode> >
      (_lrscf,guessVectors,settings.pseudoDensityThreshold,xc);
    //Kernel (only needed for first set of guess vectors (A+B))
    if(fxcSigma) F = std::make_shared<KernelSigmaVector<SCFMode> > 
      (_lrscf,guessAPB,settings.pseudoDensityThreshold,kernel);
    //External orthogonality
    //Note: All sets have the same contributions 
    //      since all EO elements in B are neglected
    //ToDo: correctly implement EO Kernel contribution unrestricted
    if(eoSigma) EO = std::make_shared<EOSigmaVector<SCFMode> > 
      (_lrscf,guessVectors,settings.pseudoDensityThreshold,
      settings.embedding.levelShiftParameter,settings.embedding.embeddingMode);

    //Add everything to sigma vectors
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd> >
      (new std::vector<Eigen::MatrixXd>(guessVectors.size()));

    for(unsigned iSet = 0; iSet < guessVectors.size(); ++iSet){
      //(A+B) and (A-B)
      (*sigma)[iSet] = D->getSigma()[iSet];
      if(kSigma) (*sigma)[iSet] -= K->getSigma()[iSet];
      if(eoSigma) (*sigma)[iSet] += factorEO * EO->getSigma()[iSet];
      //(A+B)
      if(iSet % 2 == 0){
        (*sigma)[iSet] += factor * J->getSigma()[iSet/2];
        if(fxcSigma) (*sigma)[iSet] += factor * F->getSigma()[iSet/2];
      }
    }
    return sigma;
  };/* (A+B)*b and (A-B)*b */

  //sqrt(A-B)*(A+B)*sqrt(A-B)*b
  ResponseLambda HSigma = [&] (std::vector<Eigen::MatrixXd>& guessVectors) {
    //Create sigma vectors and multiply once with sqrt(A-B)
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd> >
      (new std::vector<Eigen::MatrixXd>(guessVectors.size()));
    for(unsigned iSet = 0; iSet < guessVectors.size(); ++iSet){
      (*sigma)[iSet] = diagonal.cwiseSqrt().asDiagonal() * guessVectors[iSet];
    }
    //Orbital-energy differences
    D = std::make_shared<DeltaESigmaVector<SCFMode> > 
      (_lrscf,(*sigma),settings.pseudoDensityThreshold);
    //Coulomb
    J = std::make_shared<CoulombSigmaVector<SCFMode> >
      (_lrscf,(*sigma),settings.pseudoDensityThreshold);
    //Kernel
    if(fxcSigma) F = std::make_shared<KernelSigmaVector<SCFMode> >
      (_lrscf,(*sigma),settings.pseudoDensityThreshold,kernel);

    for (unsigned iSet = 0; iSet < guessVectors.size(); ++iSet) {
      (*sigma)[iSet] = D->getSigma()[iSet];
      (*sigma)[iSet] += factor * J->getSigma()[iSet];
      if(fxcSigma) (*sigma)[iSet] += factor * F->getSigma()[iSet];
      //Multiply again
      (*sigma)[iSet] = diagonal.cwiseSqrt().asDiagonal() * (*sigma)[iSet];
    }

    return sigma;
  };/* sqrt(A-B)*(A+B)*sqrt(A-B)*b */

  //A*b
  ResponseLambda ASigma = [&] (std::vector<Eigen::MatrixXd>& guessVectors) {

    //Orbital-energy differences
    D = std::make_shared<DeltaESigmaVector<SCFMode> >
      (_lrscf,guessVectors,settings.pseudoDensityThreshold);
    //Coulomb
    J = std::make_shared<CoulombSigmaVector<SCFMode> >
      (_lrscf,guessVectors,settings.pseudoDensityThreshold);
    //Exchange
    if(kSigma) K = std::make_shared<KSigmaVector<SCFMode> >
      (_lrscf,guessVectors,settings.pseudoDensityThreshold,std::vector<int>(guessVectors.size(),0));
    //Kernel
    if(fxcSigma) F = std::make_shared<KernelSigmaVector<SCFMode> >
      (_lrscf,guessVectors,settings.pseudoDensityThreshold,kernel);
    //External orthogonality
    if(eoSigma) EO = std::make_shared<EOSigmaVector<SCFMode> >
      (_lrscf,guessVectors, settings.pseudoDensityThreshold,
      settings.embedding.levelShiftParameter,settings.embedding.embeddingMode);

    //Add everything to sigma vectors
    auto sigma = std::unique_ptr<std::vector<Eigen::MatrixXd> >
      (new std::vector<Eigen::MatrixXd>(guessVectors.size()));

    for(unsigned iSet = 0; iSet < guessVectors.size(); ++iSet){
      (*sigma)[iSet] = D->getSigma()[iSet];
      (*sigma)[iSet] += factor / 2.0 * J->getSigma()[iSet];
      if(kSigma) (*sigma)[iSet] -= K->getSigma()[iSet];
      if(fxcSigma) (*sigma)[iSet] += factor / 2.0 * F->getSigma()[iSet];
      if(eoSigma) (*sigma)[iSet] += factorEO * EO->getSigma()[iSet];
    }
    
    return sigma;
  };/* A*b */

  /* End of lambda functions */

  auto& libint = Libint::getInstance();
  libint.keepEngines(libint2::Operator::coulomb,0,2);
  libint.keepEngines(libint2::Operator::coulomb,0,3);
  libint.keepEngines(libint2::Operator::coulomb,0,4);

  /*
    * Setup eigenvalue solver to solve the response eigenvalue problem
    *
    * Note: The EigenvalueSolver object will handle everything by itself,
    * if the TDA or pure TDDFT is used, a balanced Davidson-Liu algorithm is used and
    * if TDHF or hybrid TDDFT is used, a modified OJJ will be used
    * Note: All Hermitian subspace methods use non-orthonormal subspaces for
    *       convergence acceleration
    * 
  */
  //Set mode for eigenvalue solver
  if(settings.tda){
    //TDA: ignore user input and solve symmetrized problem
    settings.responseType = Options::RESPONSE_PROBLEM::TDA;
  }else if(settings.responseType == Options::RESPONSE_PROBLEM::RPA){
    //Keep it if the user wants it that way
  }else{
    if(settings.localMO || eoSigma || kSigma){
      //(A-B) is not diagonal, have to solve symplectic problem
      settings.responseType = Options::RESPONSE_PROBLEM::RPA;
    }else{
      //(A-B) is diagonal, might also solve the symmetrized problem
      settings.responseType = Options::RESPONSE_PROBLEM::TDDFT;
      //Transform eigenvectors if they are about to be used as an initial guess
      if (eigenvectors){
        for(auto& set : (*eigenvectors)){
          set = diagonal.cwiseSqrt().cwiseInverse().asDiagonal() * set;
        }
      }
    }
  }
  
  //Adjust maximum subspace dimension if it was untouched in the input
  if(settings.maxSubspaceDimension == 1e9) settings.maxSubspaceDimension = 25 * settings.nEigen;

  //Set initial subspace size
  //If only very few roots are to be determined, it is safer
  //to include a wider space to ensure convergence to the correct ones
  unsigned initialSubspace = settings.nEigen < 10 ? 3*settings.nEigen : 2*settings.nEigen; 
  if(initialSubspace > nDimension) initialSubspace = nDimension;

  //Note: If this is a coupled calculation and no full FDEc was requested, only perform a single iteration
  auto eigenSolver = std::make_shared<EigenvalueSolver>(settings.saveResponseMatrix,
    nDimension,settings.nEigen,diagonal,settings.convThresh,
    (type == Options::LRSCF_TYPE::COUPLED && !settings.fullFDEc) ? 1 : settings.maxCycles,
    settings.maxSubspaceDimension,initialSubspace,settings.responseType,
    settings.tda ? ASigma : settings.responseType == Options::RESPONSE_PROBLEM::TDDFT ? HSigma : MKSigma,eigenvectors);

  eigenvectors = std::make_shared<std::vector<Eigen::MatrixXd> >(eigenSolver->getEigenvectors());
  eigenvalues = eigenSolver->getEigenvalues();
  
  libint.freeEngines(libint2::Operator::coulomb,0,2);
  libint.freeEngines(libint2::Operator::coulomb,0,3);
  libint.freeEngines(libint2::Operator::coulomb,0,4);

  //Save eigenvectors and eigenvalues
  iStart = 0;
  for(auto lrscf : _lrscf){
    if (type != Options::LRSCF_TYPE::COUPLED) {
      lrscf->setSolution(eigenvectors,std::make_shared<Eigen::VectorXd> (eigenvalues),type);
    } else {
      auto nOccupied = lrscf->getNOccupied();
      auto nVirtual = lrscf->getNVirtual();
      auto vec = std::make_shared<std::vector<Eigen::MatrixXd> >(settings.tda ? 1 : 2);
      unsigned int nDimI = 0;
      for_spin(nOccupied,nVirtual) {nDimI+=nOccupied_spin*nVirtual_spin;};
      (*vec)[0] = (*eigenvectors)[0].block(iStart,0,nDimI,(*eigenvectors)[0].cols());
      if (!settings.tda) (*vec)[1] = (*eigenvectors)[1].block(iStart,0,nDimI,(*eigenvectors)[1].cols());
      iStart += nDimI;
      lrscf->setSolution(vec,std::make_shared<Eigen::VectorXd> (eigenvalues),type);
    }
  }

  //Set gauge-origin
  auto gaugeOrigin = LRSCFSetup<SCFMode>::getGaugeOrigin(settings,_act,_env);

  //Dipole integrals (electric (length/velocity) and magnetic)
  auto dipoles = std::make_shared<DipoleIntegrals<SCFMode> >(_lrscf,gaugeOrigin);

  //Analysis: Dominant Contributions, Excitation Spectrum, CD Spectrum
  if(settings.analysis){
    if(eigenvectors){
      //Prints dominant contributions etc.
      LRSCFAnalysis<SCFMode>::printDominantContributions(_lrscf,(*eigenvectors),eigenvalues,settings.dominantThresh);
      //Prints excitation and CD spectrum (oscillator and rotatory strengths, respectively)
      ExcitationSpectrum<SCFMode>::printSpectrum(dipoles,(*eigenvectors),eigenvalues);
    }
  }
}

template class LRSCFTask<Options::SCF_MODES::RESTRICTED>;
template class LRSCFTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
