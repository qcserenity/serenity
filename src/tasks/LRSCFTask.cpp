/**
 * @file   LRSCFTask.cpp
 *
 * @date   Aug 17, 2016
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
#include "tasks/LRSCFTask.h"
/* Include Serenity Internal Headers */
#include "postHF/LRSCF/Solver/DavidsonSolver.h"
#include "postHF/LRSCF/SigmaVector/DeltaESigmaVector.h"
#include "postHF/LRSCF/Analysis/ExcitationSpectrum.h"
#include "postHF/LRSCF/SigmaVector/GSigmaVector.h"
#include "io/HDF5.h"
#include "postHF/LRSCF/Solver/KDavidsonEigenvalueSolver.h"
#include "postHF/LRSCF/SigmaVector/KSigmaVector.h"
#include "postHF/LRSCF/Kernel.h"
#include "postHF/LRSCF/SigmaVector/KernelSigmaVector.h"
#include "postHF/LRSCF/Analysis/LRSCFAnalysis.h"
#include "settings/Options.h"
#include "data/OrbitalController.h"
#include "system/SystemController.h"



namespace Serenity {
using namespace std;

template<Options::SCF_MODES T> LRSCFTask<T>::LRSCFTask(
    std::shared_ptr<SystemController> activeSystem,
    std::vector<std::shared_ptr<SystemController> > environmentSystems):
      _activeSystem(activeSystem),
      _environmentSystems(environmentSystems){
  //Check input
  assert(settings.multiplicity == Options::MULTIPLICITY::SINGLET or settings.multiplicity == Options::MULTIPLICITY::TRIPLET);
  assert(_activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT or _activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::HF);
}

template<> void LRSCFTask<Options::SCF_MODES::RESTRICTED>::printSettings() {
  printSectionTitle("LRSCF");
  print((string)"Active system                         : " + _activeSystem->getSystemName());
  print((string)"Number of environment systems         : " + _environmentSystems.size());
  print((string)"Number of active electrons            : " + _nElectrons);
  print((string)"Number of active occupied orbitals    : " + _nOccupied);
  print((string)"Number of active virtual orbitals     : " + _nVirtual);
  print((string)"Number of eigenpairs to be determined : " + settings.nEigen);
  if (settings.multiplicity == Options::MULTIPLICITY::SINGLET) {
    print((string)"Multiplicity of excited states        : Singlet");
  } else if (settings.multiplicity == Options::MULTIPLICITY::TRIPLET) {
    print((string)"Multiplicity of excited states        : Triplet");
  }
  if(settings.noNaddKernel and _environmentSystems.size() > 0) {
    print((string) "\n Neglect non-additive contribution to kernel \n");
  }
}

template<> void LRSCFTask<Options::SCF_MODES::UNRESTRICTED>::printSettings() {
  printSectionTitle("LRSCF");
  print((string)"Active system                            : " + _activeSystem->getSystemName());
  print((string)"Number of environment systems            : " + _environmentSystems.size());
  print((string)"Number of active alpha electrons         : " + _nElectrons.alpha);
  print((string)"Number of active beta  electrons         : " + _nElectrons.beta);
  print((string)"Number of active occupied alpha orbitals : " + _nOccupied.alpha);
  print((string)"Number of active occupied beta orbitals  : " + _nOccupied.beta);
  print((string)"Number of active virtual alpha orbitals  : " + _nVirtual.alpha);
  print((string)"Number of active virtual beta orbitals   : " + _nVirtual.beta);
  print((string)"Number of eigenpairs to be determined    : " + settings.nEigen);
  if(settings.noNaddKernel and _environmentSystems.size() > 0) {
    print((string) "\n Neglect non-additive contribution to kernel \n");
  }
}

template<> Eigen::VectorXd LRSCFTask<Options::SCF_MODES::RESTRICTED>::estimateDiagonal() {
  Eigen::VectorXd diagonalElements(_nDimension);
  for (int i = 0; i < _nDimension; ++i) {
    //get index of orbitals given a combined index jb of Matrix M_ia,jb
    int j = floor(i/_nVirtual);
    int b = _nOccupied + i - j * _nVirtual;
    diagonalElements(i) = _orbitalEnergies(b) - _orbitalEnergies(j);
  }
  return diagonalElements;
}

template<>
Eigen::VectorXd LRSCFTask<Options::SCF_MODES::UNRESTRICTED>::estimateDiagonal() {
  Eigen::VectorXd diagonalElements(_nDimension);
  //ALPHA
  for (unsigned int i = 0; i < _nOccupied.alpha * _nVirtual.alpha; ++i) {
    //get index of orbitals given a combined index jb of Matrix M_ia,jb
    int j = floor(i/_nVirtual.alpha);
    int b = _nOccupied.alpha + i - j * _nVirtual.alpha;
    diagonalElements(i) = _orbitalEnergies.alpha(b) - _orbitalEnergies.alpha(j);
  }
  //BETA
  for (unsigned int i = 0; i < _nOccupied.beta * _nVirtual.beta; ++i) {
    //get index of orbitals given a combined index jb of Matrix M_ia,jb
    int j = floor(i/_nVirtual.beta);
    int b = _nOccupied.beta + i - j * _nVirtual.beta;
    diagonalElements(i + _nOccupied.alpha * _nVirtual.alpha) =
        _orbitalEnergies.beta(b) - _orbitalEnergies.beta(j);
  }
  return diagonalElements;
}

template<Options::SCF_MODES T>
void LRSCFTask<T>::run() {
  //get some constants (which are not defined without an electronic structure)
  _nElectrons = _activeSystem->getNElectrons<T>();
  _nOccupied = _activeSystem->getNOccupiedOrbitals<T>();
  _nVirtual = _activeSystem->getNVirtualOrbitals<T>();
  _orbitalEnergies = _activeSystem->getActiveOrbitalController<T>()->getEigenvalues();

  _nDimension = 0;
  for_spin(_nOccupied,_nVirtual) {
    _nDimension += _nOccupied_spin * _nVirtual_spin;
  };

  assert((int)settings.nEigen <= _nDimension);
  //Calculate Kernel if needed
  std::shared_ptr<Kernel<T> > kernel = nullptr;

  if (_activeSystem->getSettings().method == Options::ELECTRONIC_STRUCTURE_THEORIES::DFT || _environmentSystems.size() >= 1) {
    assert(_activeSystem->getSettings().method!=Options::ELECTRONIC_STRUCTURE_THEORIES::HF && "HF in DFT LRSCF not yet implemented");
    kernel = std::make_shared<Kernel<T> >(
        _activeSystem,
        _environmentSystems,
        settings.superSystemGrid,
        settings.noNaddKernel,
        settings.func,
        settings.naddKinFunc,
        settings.naddXCFunc);
  }

  //Print settings
  printSettings();

  //Estimate diagonal elements
  Eigen::VectorXd diagonal = estimateDiagonal();

  //
  //Functions to calculate sigma vectors
  //
  //(A-B)*b
  auto KSigma = [&] (Eigen::MatrixXd& guessVectors) {
    DeltaESigmaVector<T> deltaESigma(_activeSystem,guessVectors);
    KSigmaVector<T> kSigma(_activeSystem,guessVectors,-1.0);
    return Eigen::MatrixXd(deltaESigma.getSigma() + kSigma.getSigma());
  };
  //(A+B)*b
  auto MSigma = [&] (Eigen::MatrixXd& guessVectors) {
    Eigen::MatrixXd sigma;
    if (T == Options::SCF_MODES::UNRESTRICTED) {
      DeltaESigmaVector<T> deltaESigma(_activeSystem,guessVectors);
      GSigmaVector<T> gSigma(_activeSystem,guessVectors);
      KernelSigmaVector<T> kernelSigma(_activeSystem,guessVectors,kernel);
      sigma =  Eigen::MatrixXd(deltaESigma.getSigma() + gSigma.getSigma() + kernelSigma.getSigma());
    } else {
      if (settings.multiplicity == Options::MULTIPLICITY::SINGLET) {
        DeltaESigmaVector<T> deltaESigma(_activeSystem,guessVectors);
        GSigmaVector<T> gSigma(_activeSystem,guessVectors);
        KernelSigmaVector<T> kernelSigma(_activeSystem,guessVectors,kernel);
        sigma =  Eigen::MatrixXd(deltaESigma.getSigma() + gSigma.getSigma() + kernelSigma.getSigma());
      } else if (settings.multiplicity == Options::MULTIPLICITY::TRIPLET) {
        assert(false && "Triplet excitations not yet implemented");
      }
    }
    return sigma;
  };
  //A*b
  auto ASigma = [&] (Eigen::MatrixXd& guessVectors) {
    //ToDo: For pure density functionals, this is fine. For hybrid functionals it is better to calculate A*b directly.
    return 0.5 * (MSigma(guessVectors) + KSigma(guessVectors));
  };
  //dE*b
  auto kssSigma = [&] (Eigen::MatrixXd& guessVectors) {
    DeltaESigmaVector<T> deltaESigma(_activeSystem,guessVectors);
    return Eigen::MatrixXd(deltaESigma.getSigma());
  };


  //Solve eigenvalue problem
  std::shared_ptr<IterativeEigenvalueSolver> solver = nullptr;
  if (settings.rpa) {
    //Eigenvalue solver for RPA problem
    solver = std::make_shared<KDavidsonEigenvalueSolver>(
        _nDimension,
        settings.nEigen,
        diagonal,
        settings.convergenceCriterion,
        settings.maxCycles,
        MSigma,
        KSigma,
        settings.maxSubspaceDimension);
  } else if (settings.tda) {
    //Eigenvalue solver for TDA problem
    solver = std::make_shared<DavidsonSolver>(
        _nDimension,
        settings.nEigen,
        diagonal,
        settings.convergenceCriterion,
        settings.maxCycles,
        ASigma,
        settings.maxSubspaceDimension);
  } else if (settings.kss) {
    //Eigenvalue solver for KSS problem
    solver = std::make_shared<DavidsonSolver>(
        _nDimension,
        settings.nEigen,
        diagonal,
        settings.convergenceCriterion,
        settings.maxCycles,
        kssSigma,
        settings.maxSubspaceDimension);
  } else {
    assert(false && "No method specified. Valid options are rpa, tda or kss");
  }
  solver->solve();
  auto& eigenvectors = solver->getEigenvectors();
  auto& eigenvalues = solver->getEigenvalues();
  auto& residuals = solver->getResiduals();

  //Print excitation energies
  printTableHead(" state  ex. energy      residual norm ");
  for (unsigned int i = 0; i < eigenvalues.rows(); ++i) {
    double maxNorm = 0;
    for (unsigned int j = 0; j < residuals.size(); ++j) {
      maxNorm = std::max(maxNorm,residuals[j].col(i).norm());
    }
    std::cout << "     " << std::setw(3) << i + 1 << " " << std::setw(15) << eigenvalues(i) << " "
         << std::setw(15) << maxNorm << endl;
  }


  //Analysis
  LRSCFAnalysis<T> analysis(_activeSystem,_environmentSystems,eigenvalues,eigenvectors);
  for (unsigned int iState = 0; iState < settings.nEigen; ++iState) {
    printSmallCaption((std::string) "\n State " + (iState + 1));
    analysis.printStateInfo(iState);
    if (settings.mullikenpop) analysis.mullikenPopulationAnalysis(iState);
  }
  if (settings.writeAOVectors) analysis.printAOExcitationVectors();

  ExcitationSpectrum<T>  spec(_activeSystem,eigenvectors,eigenvalues);
  spec.printSpectrum();


  //Wirte eigenpairs to HDF5 file for restart/analysis
  std::string mode = (T==RESTRICTED) ? "res":"unres";
  std::string fName = _activeSystem->getSettings().path+_activeSystem->getSettings().name+".lrscf." + mode + ".h5";
  HDF5::H5File file(fName.c_str(), H5F_ACC_TRUNC);
  HDF5::save_scalar_attribute(file,"ID",_activeSystem->getSettings().identifier);
  HDF5::save(file,"X",eigenvectors[0]);
  if(eigenvectors.size() == 2 && settings.rpa) {
    HDF5::save(file,"Y",eigenvectors[1]);
  } else if (eigenvectors.size() == 1 && settings.tda) {
    Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(eigenvectors[0].rows(),eigenvectors[0].cols());
    HDF5::save(file,"Y",zero);
  } else {
    assert(false);
  }
  HDF5::save(file,"EIGENVALUES",eigenvalues);
  file.close();
}

template class LRSCFTask<Options::SCF_MODES::RESTRICTED> ;
template class LRSCFTask<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
