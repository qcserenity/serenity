/**
 * @file NTOCalculator.cpp
 *
 * @date Jul 13, 2017
 * @author L. Hellmann and J. Gießbach
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
#include "postHF/LRSCF/Analysis/NTOCalculator.h"
/* Include Serenity Internal Headers */
#include "data/matrices/CoefficientMatrix.h"
#include "io/HDF5.h"
#include "data/OrbitalController.h"
/* Include Std and External Headers */
#include <stdlib.h>


namespace Serenity {

template<Options::SCF_MODES T>
NTOCalculator<T>::NTOCalculator(
    std::shared_ptr<SystemController> activeSystem,
    const double plottingThreshold) :
    _activeSystem(activeSystem),
    _plottingThreshold(plottingThreshold),
    _hasBeenCalculated(false),
    _state(999999){
  //read Informations printed by the LRSCF
  printBigCaption("NTO Calculation");
  Eigen::MatrixXd x;
  Eigen::MatrixXd y;
  //Read transition densities
  if (readFromHDF5(x,y,_eigenvalues)) {
    std::cout << "Plot RPA NTOs..." << std::endl;
    _XPY = x+y;
  } else {
    assert(false && "ERROR: Could not find excitation vectors");
  }
}

template<Options::SCF_MODES T>
bool NTOCalculator<T>::readFromHDF5(
    Eigen::MatrixXd& X,
    Eigen::MatrixXd& Y,
    Eigen::VectorXd& eigenvalues) {
  auto id = _activeSystem->getSettings().identifier;
  std::string loadPath;
  std::string mode = (T==RESTRICTED) ? ".res":".unres";
  try {
    loadPath = _activeSystem->getSettings().load
        + _activeSystem->getSettings().name;
    std::string fName = loadPath + ".lrscf"+mode+".h5";
    std::cout << "\n FN " << fName << std::endl;
    HDF5::Filepath name(fName);
    HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    HDF5::attribute_exists(file, "ID");
    HDF5::check_attribute(file, "ID", id);
    std::string dataset = "X";
    HDF5::dataset_exists(file, dataset.c_str());
    HDF5::load(file, dataset.c_str(), X);
    dataset = "Y";
    HDF5::dataset_exists(file, dataset.c_str());
    HDF5::load(file, dataset.c_str(), Y);
    dataset = "EIGENVALUES";
    HDF5::dataset_exists(file, dataset.c_str());
    HDF5::load(file, dataset.c_str(), eigenvalues);
    file.close();
    return true;
  } catch (...) {
    try {
      loadPath = _activeSystem->getSettings().path
          + _activeSystem->getSettings().name;
      std::string fName = loadPath + ".lrscf"+mode+".h5";
      std::cout << "\n FN " << fName << std::endl;
      HDF5::Filepath name(fName);
      HDF5::H5File file(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      HDF5::attribute_exists(file, "ID");
      HDF5::check_attribute(file, "ID", id);
      std::string dataset = "X";
      HDF5::dataset_exists(file, dataset.c_str());
      HDF5::load(file, dataset.c_str(), X);
      dataset = "Y";
      HDF5::dataset_exists(file, dataset.c_str());
      HDF5::load(file, dataset.c_str(), Y);
      dataset = "EIGENVALUES";
      HDF5::dataset_exists(file, dataset.c_str());
      HDF5::load(file, dataset.c_str(), eigenvalues);
      file.close();
      return true;
    } catch (...) {
      return false;
    }
  }
}

template<>
void NTOCalculator<Options::SCF_MODES::RESTRICTED>::calcNTOs(int iState) {
  //Get number of occupied and virtual orbitals
  unsigned int nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::RESTRICTED>();
  unsigned int nVirt = _activeSystem->getNVirtualOrbitals<Options::SCF_MODES::RESTRICTED>();
  //Build transition density g
  Eigen::MatrixXd g(nOcc,nVirt);
  g.setZero();
  for (unsigned int i = 0, ia = 0; i < nOcc; i++) {
    for (unsigned int a = 0; a < nVirt; a++, ++ia) {
      g(i, a) = _XPY(ia, iState);
    }
  }
  //Calculate g g^T and g^T g
  Eigen::MatrixXd ggT = g * g.transpose();
  Eigen::MatrixXd gTg = g.transpose() * g;
  //Determine eigenpairs
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  //occupied problem
  es.compute(ggT);
  Eigen::MatrixXd u = es.eigenvectors();
  _occEigenvalues = es.eigenvalues();
  //virtual problem
  es.compute(gTg);
  Eigen::MatrixXd v = es.eigenvectors();
  _virtEigenvalues = es.eigenvalues();
  //Calculate NTOs
  CoefficientMatrix<Options::SCF_MODES::RESTRICTED> coeff = _activeSystem->getActiveOrbitalController<Options::SCF_MODES::RESTRICTED>()->getCoefficients();
  //occupied NTOs
  _occNTOs.resize(coeff.rows(), nOcc);
  _occNTOs.setZero();
  for (unsigned int k = 0; k < nOcc; ++k) {
    for (unsigned int i = 0; i < nOcc; ++i) {
      _occNTOs.col(k) += coeff.col(i) * u(i, k);
    }
  }
  //virtual NTOs
  unsigned int nDim = nOcc;
  if (nVirt < nOcc) nDim = nVirt;
  _virtNTOs.resize(coeff.rows(), nDim);
  _virtNTOs.setZero();
  for (unsigned int l = 0; l < nDim; ++l) {
   for (unsigned int a = 0; a < nVirt; ++a) {
     _virtNTOs.col(l) += coeff.col(a + nOcc) * v(a, v.cols() - l - 1);
   }
  }
  //print additional data to file
  printNTOInfo("RESTRICTED", iState, _occEigenvalues, _virtEigenvalues, u, v);
  _hasBeenCalculated = true;
  _state = iState;
}


template<>
void NTOCalculator<Options::SCF_MODES::UNRESTRICTED>::calcNTOs(int iState) {
  //Get number of occupied and virtual orbitals
  auto nOcc = _activeSystem->getNOccupiedOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  auto nVirt = _activeSystem->getNVirtualOrbitals<Options::SCF_MODES::UNRESTRICTED>();
  //Build transition density g
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> g;
  for_spin(g,nOcc,nVirt) {
    g_spin.resize(nOcc_spin,nVirt_spin);
    g_spin.setZero();
  };
  for (unsigned int i = 0, ia = 0; i < nOcc.alpha; i++) {
    for (unsigned int a = 0; a < nVirt.alpha; a++, ++ia) {
      g.alpha(i, a) = _XPY(ia, iState);
    }
  }
  for (unsigned int i = 0, ia = nOcc.alpha * nVirt.alpha; i < nOcc.beta; i++) {
    for (unsigned int a = 0; a < nVirt.beta; a++, ++ia) {
      g.beta(i, a) = _XPY(ia, iState);
    }
  }
  //Calculate g g^T and g^T g
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> ggT;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> gTg;
  for_spin(g,ggT,gTg) {
    ggT_spin = g_spin * g_spin.transpose();
    gTg_spin = g_spin.transpose() * g_spin;
  };
  //Determine eigenpairs
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> u;
  SpinPolarizedData<Options::SCF_MODES::UNRESTRICTED, Eigen::MatrixXd> v;
  //occupied problem
  for_spin(ggT,u,_occEigenvalues) {
    es.compute(ggT_spin);
    u_spin = es.eigenvectors();
    _occEigenvalues_spin = es.eigenvalues();
  };
  //virtual problem
  for_spin(gTg,v,_virtEigenvalues) {
    es.compute(gTg_spin);
    v_spin = es.eigenvectors();
    _virtEigenvalues_spin = es.eigenvalues();
  };
  //Calculate NTOs
  CoefficientMatrix<Options::SCF_MODES::UNRESTRICTED> coeff = _activeSystem->getActiveOrbitalController<Options::SCF_MODES::UNRESTRICTED>()->getCoefficients();
  //occupied NTOs
  for_spin(_occNTOs,coeff,nOcc,u) {
    _occNTOs_spin.resize(coeff_spin.rows(),nOcc_spin);
    _occNTOs_spin.setZero();
    for (unsigned int k = 0; k < nOcc_spin; ++k) {
      for (unsigned int i = 0; i < nOcc_spin; ++i) {
        _occNTOs_spin.col(k) += coeff_spin.col(i) * u_spin(i, k);
      }
    }
  };
  //virtual NTOs
  for_spin(_virtNTOs,coeff,nOcc,nVirt,v) {
    unsigned int nDim = nOcc_spin;
    if (nVirt_spin < nOcc_spin) nDim = nVirt_spin;
    _virtNTOs_spin.resize(coeff_spin.rows(),nDim);
    _virtNTOs_spin.setZero();
    for (unsigned int l = 0; l < nDim; ++l) {
      for (unsigned int a = 0; a < nVirt_spin; ++a) {
        _virtNTOs_spin.col(l) += coeff_spin.col(a + nOcc_spin) * v_spin(a, v_spin.cols() - l - 1);
      }
    }
  };
  printNTOInfo("ALPHA", iState, _occEigenvalues.alpha, _virtEigenvalues.alpha, u.alpha, v.alpha);
  printNTOInfo("BETA", iState, _occEigenvalues.beta, _virtEigenvalues.beta, u.beta, v.beta);
  _hasBeenCalculated = true;
  _state = iState;
}

template<Options::SCF_MODES T>
void NTOCalculator<T>::printNTOInfo(
    std::string spin,
    const unsigned int iState,
    Eigen::VectorXd& oEigenvalues,
    Eigen::VectorXd& vEigenvalues,
    Eigen::MatrixXd& u,
    Eigen::MatrixXd& v){
  printSmallCaption((std::string)"\n NTOs for state " + (iState + 1) + " " + spin);
  //Make directory
  std::string dirName = _activeSystem->getSettings().path+"/NTOS" + "/NTO" + std::to_string(iState+1) + "/";
  std::string command="mkdir -p "+ dirName;
  auto stat = system(command.c_str());
  (void) stat;
  print((std::string)"Plotting NTOs to " + dirName + "... \n");
  //create or update readme file
  FILE *readme;
  readme = fopen ((dirName+"README").c_str(),"a");
  //Print occupied eigenvalues to readme
  fprintf (readme, "\n %s : \n",spin.c_str());
  fprintf (readme, "%s \n","OCCUPIED EIGENVALUES");
  print((std::string)"Occupied eigenvalues : ");
  for (unsigned int i = 0; i < oEigenvalues.rows(); ++i) {
    if (oEigenvalues(i) > _plottingThreshold) {
      printf("   NTO %3i %3.2f \n ", i + 1, oEigenvalues(i));
      fprintf (readme, "NTO %d %f \n",i+1,oEigenvalues(i));
    }
  }
  //Print occupied dominant contributions to readme
  fprintf (readme, "\n %s ","DOMINANT CONTRIBUTIONS");
  print((std::string)"\n Occupied dominant contributions : ");
  for (unsigned int i = 0; i < u.cols(); ++i) {
    if (oEigenvalues(i) > _plottingThreshold) {
      printf("   NTO %3.i :", i+1);
      fprintf (readme, "\n NTO %d : ",i+1);
      for (unsigned int j = 0; j < u.rows(); ++j) {
        if (fabs(u(j,i)) > _plottingThreshold) {
          printf("%3.i (%5.3f) ", j+1, u(j,i));
          fprintf (readme, "%d (%5.3f)  ",j+1,u(j,i));
        }
      }
      printf("\n");
    }
  }
  //Print virtual eigenvalues to readme
  fprintf (readme, "\n \n %s \n","VIRTUAL EIGENVALUES");
  print((std::string)"\n Virtual eigenvalues : ");
  for (unsigned int a = 0; a < vEigenvalues.rows(); ++a) {
    if (vEigenvalues((int)(vEigenvalues.rows() - a - 1)) > _plottingThreshold) {
      printf("   NTO %3i %3.2f \n ", (int)(oEigenvalues.rows() + a+1) , vEigenvalues((int)(vEigenvalues.rows() - a - 1)));
      fprintf (readme, "NTO %d %f \n",(int)(oEigenvalues.rows() + a+1),vEigenvalues((int)(vEigenvalues.rows() - a - 1)));
    }
  }
  //Print virtual dominant contributions to readme
  fprintf (readme, "\n %s ","DOMINANT CONTRIBUTIONS");
  print((std::string)"\n Virtual dominant contributions : ");
  for (unsigned int a = 0; a < vEigenvalues.rows(); ++a) {
    if (vEigenvalues((int)(vEigenvalues.rows() - a - 1)) > _plottingThreshold) {
      printf("   NTO %3.i :", (int)(oEigenvalues.rows() + a+1));
      fprintf (readme, "\n NTO %d : ",(int)(oEigenvalues.rows() + a+1));
      for (unsigned int j = 0; j < v.rows(); ++j) {
        if (fabs(v(j,v.cols() - a - 1)) > _plottingThreshold) {
          printf("%3.i (%5.3f) ",(int)(u.cols() + j+1),v(j,(int)(v.cols() - a - 1)));
          fprintf (readme, "%d (%5.3f)  ",(int)(u.cols() + j+1),v(j,(int)(v.cols() - a - 1)));
        }
      }
      printf("\n");
    }
  }
  fclose(readme);
}

template<Options::SCF_MODES T>
unsigned int NTOCalculator<T>::getNumberOfStates(){
  return _eigenvalues.rows();
}

template<Options::SCF_MODES T>
std::string NTOCalculator<T>::getDir(unsigned int i){
  std::string dirName = _activeSystem->getSettings().path+"/NTOS" + "/NTO" + std::to_string(i+1) + "/";
  return dirName;
}

template class NTOCalculator<Options::SCF_MODES::RESTRICTED>;
template class NTOCalculator<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
