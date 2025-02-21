/**
 * @file NumericalHessianCalc.cpp
 *
 * @date Feb 23, 2017
 * @author Kevin Klahr
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
#include "geometry/gradients/NumericalHessianCalc.h"
/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "geometry/Atom.h"
#include "geometry/AtomType.h"
#include "geometry/gradients/NumericalGeomGradCalc.h"
#include "integrals/wrappers/Libint.h"
#include "io/HDF5.h"
#include "math/Matrix.h"
#include "math/linearAlgebra/Orthogonalization.h"
#include "parameters/Constants.h"
#include "settings/EmbeddingSettings.h"
#include "settings/Settings.h"
#include "system/SystemController.h"
#include "tasks/FreezeAndThawTask.h"
#include "tasks/GradientTask.h"
#include "tasks/ScfTask.h"
/* Include Std and External Headers */
#include <unistd.h>
#include <iomanip>

namespace Serenity {

template<Options::SCF_MODES SCFMode>
NumericalHessianCalc<SCFMode>::NumericalHessianCalc(double stepsizeGrad, double stepsizeHess, bool printToFile)
  : _deltaGrad(stepsizeGrad), _deltaHess(stepsizeHess), _printToFile(printToFile) {
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd NumericalHessianCalc<SCFMode>::calcHessian(std::shared_ptr<SystemController> systemController) {
  auto gradtype = Options::GRADIENT_TYPES::ANALYTICAL;

  if (_deltaGrad > 0.0) {
    std::cout << "Numerical Hessian will be calculated!" << std::endl;
    gradtype = Options::GRADIENT_TYPES::NUMERICAL;
  }
  else {
    std::cout << "Semiumerical Hessian will be calculated with analytical first derivatives!" << std::endl;
  }

  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 2);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 3);
  libint.keepEngines(LIBINT_OPERATOR::coulomb, 1, 4);

  auto geometry = systemController->getGeometry();

  double delta = _deltaHess;
  unsigned int nAtoms = geometry->getNAtoms();

  Matrix<double> hessian(3 * nAtoms, 3 * nAtoms);
  hessian.setZero();

  for (unsigned int i = 0; i != nAtoms; ++i) {
    geometry->getAtoms()[i]->addToX(delta);
    auto scf = ScfTask<SCFMode>(systemController);
    scf.run();
    {
      GradientTask<SCFMode> task({systemController});
      task.settings.gradType = gradtype;
      task.settings.numGradStepSize = _deltaGrad;
      task.settings.print = false;
      task.run();
    }
    auto grad = systemController->getGeometry()->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i, 3 * j + 0) += grad(j, 0); // XX
      hessian(3 * i, 3 * j + 1) += grad(j, 1); // XY
      hessian(3 * i, 3 * j + 2) += grad(j, 2); // XZ
    }
    geometry->getAtoms()[i]->addToX(-2.0 * delta);
    scf.run();
    {
      GradientTask<SCFMode> task({systemController});
      task.settings.gradType = gradtype;
      task.settings.numGradStepSize = _deltaGrad;
      task.settings.print = false;
      task.run();
    }
    grad = systemController->getGeometry()->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i, 3 * j + 0) -= grad(j, 0);
      hessian(3 * i, 3 * j + 1) -= grad(j, 1);
      hessian(3 * i, 3 * j + 2) -= grad(j, 2);
      hessian(3 * i, 3 * j + 0) /= 2.0 * delta;
      hessian(3 * i, 3 * j + 1) /= 2.0 * delta;
      hessian(3 * i, 3 * j + 2) /= 2.0 * delta;
    }
    geometry->getAtoms()[i]->addToX(delta);
  }

  for (unsigned int i = 0; i != nAtoms; ++i) {
    geometry->getAtoms()[i]->addToY(delta);
    auto scf = ScfTask<SCFMode>(systemController);
    scf.run();
    {
      GradientTask<SCFMode> task({systemController});
      task.settings.gradType = gradtype;
      task.settings.numGradStepSize = _deltaGrad;
      task.settings.print = false;
      task.run();
    }
    auto grad = systemController->getGeometry()->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 1, 3 * j + 0) += grad(j, 0); // YX
      hessian(3 * i + 1, 3 * j + 1) += grad(j, 1); // YY
      hessian(3 * i + 1, 3 * j + 2) += grad(j, 2); // YZ
    }
    geometry->getAtoms()[i]->addToY(-2.0 * delta);
    scf.run();
    {
      GradientTask<SCFMode> task({systemController});
      task.settings.gradType = gradtype;
      task.settings.numGradStepSize = _deltaGrad;
      task.settings.print = false;
      task.run();
    }
    grad = systemController->getGeometry()->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 1, 3 * j + 0) -= grad(j, 0);
      hessian(3 * i + 1, 3 * j + 1) -= grad(j, 1);
      hessian(3 * i + 1, 3 * j + 2) -= grad(j, 2);
      hessian(3 * i + 1, 3 * j + 0) /= 2.0 * delta;
      hessian(3 * i + 1, 3 * j + 1) /= 2.0 * delta;
      hessian(3 * i + 1, 3 * j + 2) /= 2.0 * delta;
    }
    geometry->getAtoms()[i]->addToY(delta);
  }

  for (unsigned int i = 0; i != nAtoms; ++i) {
    geometry->getAtoms()[i]->addToZ(delta);
    auto scf = ScfTask<SCFMode>(systemController);
    scf.run();
    {
      GradientTask<SCFMode> task({systemController});
      task.settings.gradType = gradtype;
      task.settings.numGradStepSize = _deltaGrad;
      task.settings.print = false;
      task.run();
    }
    auto grad = systemController->getGeometry()->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 2, 3 * j + 0) += grad(j, 0); // ZX
      hessian(3 * i + 2, 3 * j + 1) += grad(j, 1); // ZY
      hessian(3 * i + 2, 3 * j + 2) += grad(j, 2); // ZZ
    }
    geometry->getAtoms()[i]->addToZ(-2.0 * delta);
    scf.run();
    {
      GradientTask<SCFMode> task({systemController});
      task.settings.gradType = gradtype;
      task.settings.numGradStepSize = _deltaGrad;
      task.settings.print = false;
      task.run();
    }
    grad = systemController->getGeometry()->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 2, 3 * j + 0) -= grad(j, 0);
      hessian(3 * i + 2, 3 * j + 1) -= grad(j, 1);
      hessian(3 * i + 2, 3 * j + 2) -= grad(j, 2);
      hessian(3 * i + 2, 3 * j + 0) /= 2.0 * delta;
      hessian(3 * i + 2, 3 * j + 1) /= 2.0 * delta;
      hessian(3 * i + 2, 3 * j + 2) /= 2.0 * delta;
    }
    geometry->getAtoms()[i]->addToZ(delta);
  }
  // reset electronic structure.
  {
    GradientTask<SCFMode> task({systemController});
    task.settings.gradType = gradtype;
    task.settings.numGradStepSize = _deltaGrad;
    task.settings.print = false;
    task.run();
  }
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 0, 4);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 2);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 3);
  libint.freeEngines(LIBINT_OPERATOR::coulomb, 1, 4);

  auto& settings = systemController->getSettings();

  // print hessian
  std::cout << "\n Cartesian Hessian:" << std::endl;
  for (unsigned int hessBlock = 0; hessBlock < nAtoms; hessBlock++) {
    std::cout << "            " << 3 * hessBlock + 1 << "                   " << 3 * hessBlock + 2
              << "                   " << 3 * hessBlock + 3 << std::endl;
    for (unsigned int i = 0; i != 3 * nAtoms; ++i) {
      std::cout << i + 1 << ": " << hessian(i, 3 * hessBlock + 0) << " " << hessian(i, 3 * hessBlock + 1) << " "
                << hessian(i, 3 * hessBlock + 2) << std::endl;
    }
    std::cout << "\n";
  }
  frequencyCalculation(hessian, geometry, {settings});
  return std::move(hessian);
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd NumericalHessianCalc<SCFMode>::calcFaTHessian(std::vector<std::shared_ptr<SystemController>> activeSystems,
                                                              std::vector<std::shared_ptr<SystemController>> passiveSystems,
                                                              EmbeddingSettings embedding, int FatmaxCycles,
                                                              double FaTenergyConvThresh, double FaTgridCutOff) {
  double delta = _deltaHess;
  unsigned int nAtoms;
  std::shared_ptr<Geometry> systemGeometry;

  auto supersystemGeometry = std::make_shared<Geometry>();
  for (auto sys : activeSystems) {
    *supersystemGeometry += *sys->getGeometry();
  }
  systemGeometry = supersystemGeometry;

  nAtoms = systemGeometry->getNAtoms();

  Matrix<double> hessian(3 * nAtoms, 3 * nAtoms);
  hessian.setZero();
  for (unsigned int i = 0; i != nAtoms; ++i) {
    systemGeometry->getAtoms()[i]->addToX(delta);

    // Active system cycles
    {
      GradientTask<SCFMode> task(activeSystems, passiveSystems);
      task.settings.embedding = embedding;
      task.settings.FaTmaxCycles = FatmaxCycles;
      task.settings.FaTenergyConvThresh = FaTenergyConvThresh;
      task.settings.FDEgridCutOff = FaTgridCutOff;
      task.settings.print = false;
      task.run();
    }

    auto grad = systemGeometry->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i, 3 * j + 0) += grad(j, 0); // XX
      hessian(3 * i, 3 * j + 1) += grad(j, 1); // XY
      hessian(3 * i, 3 * j + 2) += grad(j, 2); // XZ
    }
    systemGeometry->getAtoms()[i]->addToX(-2.0 * delta);

    // Active system cycles
    {
      GradientTask<SCFMode> task(activeSystems, passiveSystems);
      task.settings.embedding = embedding;
      task.settings.FaTmaxCycles = FatmaxCycles;
      task.settings.FaTenergyConvThresh = FaTenergyConvThresh;
      task.settings.FDEgridCutOff = FaTgridCutOff;
      task.settings.print = false;
      task.run();
    }

    grad = systemGeometry->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i, 3 * j + 0) -= grad(j, 0);
      hessian(3 * i, 3 * j + 1) -= grad(j, 1);
      hessian(3 * i, 3 * j + 2) -= grad(j, 2);
      hessian(3 * i, 3 * j + 0) /= 2.0 * delta;
      hessian(3 * i, 3 * j + 1) /= 2.0 * delta;
      hessian(3 * i, 3 * j + 2) /= 2.0 * delta;
    }
    systemGeometry->getAtoms()[i]->addToX(delta);
  }

  for (unsigned int i = 0; i != nAtoms; ++i) {
    systemGeometry->getAtoms()[i]->addToY(delta);

    // Active system cycles
    {
      GradientTask<SCFMode> task(activeSystems, passiveSystems);
      task.settings.embedding = embedding;
      task.settings.FaTmaxCycles = FatmaxCycles;
      task.settings.FaTenergyConvThresh = FaTenergyConvThresh;
      task.settings.FDEgridCutOff = FaTgridCutOff;
      task.settings.print = false;
      task.run();
    }

    auto grad = systemGeometry->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 1, 3 * j + 0) += grad(j, 0); // YX
      hessian(3 * i + 1, 3 * j + 1) += grad(j, 1); // YY
      hessian(3 * i + 1, 3 * j + 2) += grad(j, 2); // YZ
    }
    systemGeometry->getAtoms()[i]->addToY(-2.0 * delta);

    // Active system cycles
    {
      GradientTask<SCFMode> task(activeSystems, passiveSystems);
      task.settings.embedding = embedding;
      task.settings.FaTmaxCycles = FatmaxCycles;
      task.settings.FaTenergyConvThresh = FaTenergyConvThresh;
      task.settings.FDEgridCutOff = FaTgridCutOff;
      task.settings.print = false;
      task.run();
    }

    grad = systemGeometry->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 1, 3 * j + 0) -= grad(j, 0);
      hessian(3 * i + 1, 3 * j + 1) -= grad(j, 1);
      hessian(3 * i + 1, 3 * j + 2) -= grad(j, 2);
      hessian(3 * i + 1, 3 * j + 0) /= 2.0 * delta;
      hessian(3 * i + 1, 3 * j + 1) /= 2.0 * delta;
      hessian(3 * i + 1, 3 * j + 2) /= 2.0 * delta;
    }
    systemGeometry->getAtoms()[i]->addToY(delta);
  }

  for (unsigned int i = 0; i != nAtoms; ++i) {
    systemGeometry->getAtoms()[i]->addToZ(delta);

    // Active system cycles
    {
      GradientTask<SCFMode> task(activeSystems, passiveSystems);
      task.settings.embedding = embedding;
      task.settings.FaTmaxCycles = FatmaxCycles;
      task.settings.FaTenergyConvThresh = FaTenergyConvThresh;
      task.settings.FDEgridCutOff = FaTgridCutOff;
      task.settings.print = false;
      task.run();
    }

    auto grad = systemGeometry->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 2, 3 * j + 0) += grad(j, 0); // ZX
      hessian(3 * i + 2, 3 * j + 1) += grad(j, 1); // ZY
      hessian(3 * i + 2, 3 * j + 2) += grad(j, 2); // ZZ
    }
    systemGeometry->getAtoms()[i]->addToZ(-2.0 * delta);

    // Active system cycles
    {
      GradientTask<SCFMode> task(activeSystems, passiveSystems);
      task.settings.embedding = embedding;
      task.settings.FaTmaxCycles = FatmaxCycles;
      task.settings.FaTenergyConvThresh = FaTenergyConvThresh;
      task.settings.FDEgridCutOff = FaTgridCutOff;
      task.settings.print = false;
      task.run();
    }

    grad = systemGeometry->getGradients();
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 2, 3 * j + 0) -= grad(j, 0);
      hessian(3 * i + 2, 3 * j + 1) -= grad(j, 1);
      hessian(3 * i + 2, 3 * j + 2) -= grad(j, 2);
      hessian(3 * i + 2, 3 * j + 0) /= 2.0 * delta;
      hessian(3 * i + 2, 3 * j + 1) /= 2.0 * delta;
      hessian(3 * i + 2, 3 * j + 2) /= 2.0 * delta;
    }
    systemGeometry->getAtoms()[i]->addToZ(delta);
  }

  std::vector<Settings> settings;
  for (unsigned int i = 0; i < activeSystems.size(); i++) {
    settings.push_back(activeSystems[i]->getSettings());
  }

  // print hessian
  std::cout << "\n Cartesian Hessian:" << std::endl;
  for (unsigned int hessBlock = 0; hessBlock < nAtoms; hessBlock++) {
    std::cout << "            " << 3 * hessBlock + 1 << "                   " << 3 * hessBlock + 2
              << "                   " << 3 * hessBlock + 3 << std::endl;
    for (unsigned int i = 0; i != 3 * nAtoms; ++i) {
      std::cout << i + 1 << ": " << hessian(i, 3 * hessBlock + 0) << " " << hessian(i, 3 * hessBlock + 1) << " "
                << hessian(i, 3 * hessBlock + 2) << std::endl;
    }
    std::cout << "\n";
  }

  frequencyCalculation(hessian, systemGeometry, settings);
  return std::move(hessian);
}

template<Options::SCF_MODES SCFMode>
void NumericalHessianCalc<SCFMode>::frequencyCalculation(Eigen::MatrixXd hessian, std::shared_ptr<Geometry> geometry,
                                                         const std::vector<Settings> settings) {
  Eigen::MatrixXd rawHessian = hessian;
  unsigned int nAtoms = geometry->getNAtoms();
  // Mass weighting
  auto atoms = geometry->getAtoms();
  for (unsigned int i = 0; i != nAtoms; ++i) {
    for (unsigned int j = 0; j != nAtoms; ++j) {
      hessian(3 * i + 0, 3 * j + 0) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 0, 3 * j + 1) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 0, 3 * j + 2) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 1, 3 * j + 0) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 1, 3 * j + 1) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 1, 3 * j + 2) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 2, 3 * j + 0) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 2, 3 * j + 1) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
      hessian(3 * i + 2, 3 * j + 2) *=
          1.0 / (sqrt(atoms[i]->getAtomType()->getMass()) * sqrt(atoms[j]->getAtomType()->getMass()));
    }
  }

  // Move center of mass of the system to cartesian origin

  std::vector<double> mCoords(3, 0.0);
  double molMass = 0.0;
  for (unsigned int i = 0; i != nAtoms; ++i) {
    mCoords[0] += atoms[i]->getX() * atoms[i]->getAtomType()->getMass();
    mCoords[1] += atoms[i]->getY() * atoms[i]->getAtomType()->getMass();
    mCoords[2] += atoms[i]->getZ() * atoms[i]->getAtomType()->getMass();
    molMass += atoms[i]->getAtomType()->getMass();
  }
  mCoords[0] = mCoords[0] / molMass;
  mCoords[1] = mCoords[1] / molMass;
  mCoords[2] = mCoords[2] / molMass;
  for (unsigned int i = 0; i != nAtoms; ++i) {
    atoms[i]->setX(atoms[i]->getX() - mCoords[0]);
    atoms[i]->setY(atoms[i]->getY() - mCoords[1]);
    atoms[i]->setZ(atoms[i]->getZ() - mCoords[2]);
  }

  Eigen::MatrixXd transModes = geometry->getTransModes();

  Eigen::MatrixXd rotModes = geometry->getRotModes();

  // eliminate trans components from projection operator

  Eigen::MatrixXd c2p(3 * nAtoms, 3 * nAtoms);
  c2p = Eigen::MatrixXd::Identity(3 * nAtoms, 3 * nAtoms);

  for (unsigned int i = 0; i < 3; i++) {
    Eigen::Map<Eigen::VectorXd> colVec(transModes.col(i).data(), transModes.rows());
    c2p -= colVec * colVec.transpose();
  }

  // eliminate rot components from projection operator

  for (unsigned int i = 0; i < 3; i++) {
    Eigen::Map<Eigen::VectorXd> colVec(rotModes.col(i).data(), rotModes.rows());
    c2p -= colVec * colVec.transpose();
  }

  Eigen::MatrixXd c2p_ortho(3 * nAtoms, 3 * nAtoms);

  // Gram Schmidt procedure

  Orthogonalization::modifiedGramSchmidtHessianProjection(c2p, c2p_ortho);

  Eigen::MatrixXd c;
  unsigned int dim1 = 3 * nAtoms;
  unsigned int dim2 = 3 * nAtoms;

  bool linearCheck = geometry->isLinear();
  c = c2p_ortho.block(0, 0, dim1, dim2 - (linearCheck ? 5 : 6));

  Eigen::MatrixXd hessian_p = c.transpose() * hessian * c;

  // Symmetrize Hessian

  for (unsigned int j = 0; j < (dim2 - (linearCheck ? 5 : 6)); j++) {
    for (unsigned int i = 0; i < j; i++) {
      hessian_p(i, j) = 0.5 * (hessian_p(i, j) + hessian_p(j, i));
      hessian_p(j, i) = hessian_p(i, j);
    }
  }

  // Diagonalizing

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenvalueSolver(3 * nAtoms - 6);
  eigenvalueSolver.compute(hessian_p);
  Eigen::VectorXd exactEigenvalues = eigenvalueSolver.eigenvalues();
  auto eigenVectors = eigenvalueSolver.eigenvectors();

  // Convert frequencies units to wavenumbers

  Eigen::VectorXd frequencies(exactEigenvalues.size());
  frequencies.setZero();
  double atu2sec = (1E-06 * (BOHR * BOHR)) / HARTREE_TO_KJ_PER_MOL;

  for (unsigned int i = 0; i < exactEigenvalues.size(); i++) {
    if (exactEigenvalues[i] < 0.0) {
      frequencies[i] = exactEigenvalues[i] / atu2sec;
      frequencies[i] = -((1.0 / (2.0 * PI)) * (sqrt(-frequencies[i]) / (100 * SPEEDOFLIGHT)));
    }
    else {
      frequencies[i] = exactEigenvalues[i] / atu2sec;
      frequencies[i] = (1.0 / (2.0 * PI)) * (sqrt(frequencies[i]) / (100 * SPEEDOFLIGHT));
    }
  }

  // Get normal modes
  Eigen::MatrixXd normalModes(c.rows(), exactEigenvalues.size());
  for (unsigned int i = 0; i < exactEigenvalues.size(); i++) {
    Eigen::Map<Eigen::VectorXd> eVec(eigenVectors.col(i).data(), eigenVectors.rows());
    normalModes.col(i) = c * eVec;
  }

  // Print modes
  std::cout << "\nFrequencies and mass-weighted normal coordinates:" << std::endl;
  std::cout << std::fixed << std::setprecision(5);
  for (unsigned int freq = 0; freq < frequencies.size(); freq++) {
    std::cout << "\n\nRoot #: " << freq + 1 << std::endl;
    std::cout << "Eigval: " << frequencies(freq) << "\n\n";
    for (unsigned int normCoord = 0; normCoord < nAtoms; normCoord++) {
      std::cout << "x  " << geometry->getAtoms()[normCoord]->getAtomType()->getElementSymbol() << "   " << normCoord + 1
                << "   " << normalModes(3 * normCoord + 0, freq) << std::endl;
      std::cout << "y  " << geometry->getAtoms()[normCoord]->getAtomType()->getElementSymbol() << "   " << normCoord + 1
                << "   " << normalModes(3 * normCoord + 1, freq) << std::endl;
      std::cout << "z  " << geometry->getAtoms()[normCoord]->getAtomType()->getElementSymbol() << "   " << normCoord + 1
                << "   " << normalModes(3 * normCoord + 2, freq) << std::endl;
    }
  }
  std::cout << std::scientific << std::setprecision(12);
  std::vector<const char*> systemNames;

  // write Eigenvectors/eigenvalues to HDF5
  std::string name = settings[0].path + settings[0].name + ".hess.h5";
  HDF5::H5File file(name.c_str(), H5F_ACC_TRUNC);
  HDF5::save(file, "eigenvectors", eigenVectors);
  HDF5::save(file, "eigenvalues", exactEigenvalues);
  HDF5::save(file, "frequencies", frequencies);
  HDF5::save(file, "normalModes", normalModes);
  HDF5::save_scalar_attribute(file, "ID", settings[0].identifier);

  // make sure there is a link to the supersystem hessian data in all the subsystem settings folders
  systemNames.push_back(settings[0].name.c_str());
  for (unsigned int i = 1; i < settings.size(); i++) {
    std::string linkName = settings[i].path + settings[i].name + ".hess.h5";
    systemNames.push_back(settings[i].name.c_str());
    assert(0 == link((name).c_str(), linkName.c_str()) && "Error: Failed to create hardlink for hessian output.");
  }

  HDF5::save_std_vector(file, "names", systemNames);

  // write Hessian to HDF5

  HDF5::save(file, "hessian", hessian);
  HDF5::save(file, "rawHessian", rawHessian);
  file.close();

  // write molden input

  if (_printToFile) {
    std::ofstream output;
    auto moldenFilePath = settings[0].path + "molden.input";
    output.open(moldenFilePath);
    output << "[Molden Format]\n";
    output << "[Title]\n";
    output << "\n";
    output << "[Atoms] AU\n";
    unsigned int atomIndex = 0;
    for (auto atom : geometry->getAtoms()) {
      atomIndex++;
      output << atom->getAtomType()->getElementSymbol() << "   " << atomIndex << "   "
             << atom->getAtomType()->getNuclearCharge() << " " << std::to_string(atom->getX()) << " "
             << std::to_string(atom->getY()) << " " << std::to_string(atom->getZ()) << "\n";
    }
    output << "[FREQ]\n";
    for (unsigned int freq = 0; freq < frequencies.size(); freq++) {
      output << frequencies[freq] << "\n";
    }
    output << "[FR-COORD]\n";
    for (auto atom : geometry->getAtoms()) {
      output << atom->getAtomType()->getElementSymbol() << " " << std::to_string(atom->getX()) << " "
             << std::to_string(atom->getY()) << " " << std::to_string(atom->getZ()) << "\n";
    }
    output << "[FR-NORM-COORD]\n";
    for (unsigned int i = 0; i < frequencies.size(); i++) {
      output << "vibration   " << i + 1 << "\n";
      atomIndex = 0;
      for (auto atom : geometry->getAtoms()) {
        output << normalModes(3 * atomIndex + 0, i) << " " << normalModes(3 * atomIndex + 1, i) << " "
               << normalModes(3 * atomIndex + 2, i) << "\n";
        atomIndex++;
      }
    }
    output.close();
  }
}

template class NumericalHessianCalc<Options::SCF_MODES::RESTRICTED>;
template class NumericalHessianCalc<Options::SCF_MODES::UNRESTRICTED>;
} /* namespace Serenity */
