/**
 * @file ElectronicStructure_python.cpp
 *
 * @date Mar 20, 2017
 * @author Jan Unsleber
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

/* Include Serenity Internal Headers */
#include "data/ElectronicStructure.h"
#include "data/OrbitalController.h"
#include "energies/EnergyComponentController.h"
#include "integrals/OneElectronIntegralController.h"
/* Include Std and External Headers */
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace Serenity;

Eigen::MatrixXd overlap(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>> es) {
  return es->getOneElectronIntegralController()->getOverlapIntegrals();
}
/*
 * Density Matrix
 */
template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd denstotal(std::shared_ptr<ElectronicStructure<SCFMode>> es) {
  return es->getDensityMatrix().total();
}
Eigen::MatrixXd densalpha(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getDensityMatrix().alpha;
}
Eigen::MatrixXd densbeta(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getDensityMatrix().beta;
}
/*
 * Orbital Coefficients
 */
Eigen::MatrixXd coeff(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>> es) {
  return es->getMolecularOrbitals()->getCoefficients();
}
Eigen::MatrixXd coeffalpha(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getMolecularOrbitals()->getCoefficients().alpha;
}
Eigen::MatrixXd coeffbeta(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getMolecularOrbitals()->getCoefficients().beta;
}
/*
 * Orbital Energies
 */
Eigen::VectorXd orbEn(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>> es) {
  return es->getMolecularOrbitals()->getEigenvalues();
}
Eigen::VectorXd orbEnalpha(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getMolecularOrbitals()->getEigenvalues().alpha;
}
Eigen::VectorXd orbEnbeta(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getMolecularOrbitals()->getEigenvalues().beta;
}
/*
 * Fock Matrix
 */
Eigen::MatrixXd fock(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>> es) {
  return es->getFockMatrix();
}
Eigen::MatrixXd fockalpha(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getFockMatrix().alpha;
}
Eigen::MatrixXd fockbeta(std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>> es) {
  return es->getFockMatrix().beta;
}
/*
 * Class
 */
void export_ElectronicStructure(py::module& spy) {
  py::class_<ElectronicStructure<Options::SCF_MODES::RESTRICTED>, std::shared_ptr<ElectronicStructure<Options::SCF_MODES::RESTRICTED>>>(
      spy, "ElectronicStructure_R", "@brief A restricted electronic structure.")
      .def("getDensityMatrix", &ElectronicStructure<Options::SCF_MODES::RESTRICTED>::getDensityMatrixController,
           "@brief Returns the current density matrix(controller).\n"
           "\n"
           "@returns The density matrix.")
      .def("totalDens", &denstotal<RESTRICTED>)
      .def("overlap", &overlap)
      .def("coeff", &coeff)
      .def("orbEn", &orbEn)
      .def("fock", &fock);
  py::class_<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>, std::shared_ptr<ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>>>(
      spy, "ElectronicStructure_U", "@brief An unrestricted electronic structure.")
      .def("getDensityMatrix", &ElectronicStructure<Options::SCF_MODES::UNRESTRICTED>::getDensityMatrixController,
           "@brief Returns the current density matrix(controller).\n"
           "\n"
           "@returns The density matrix.")
      .def("totalDens", &denstotal<UNRESTRICTED>)
      .def("alphaDens", &densalpha)
      .def("betaDens", &densbeta)
      .def("overlap", &overlap)
      .def("alphaCoeff", &coeffalpha)
      .def("betaCoeff", &coeffbeta)
      .def("alphaOrbEn", &orbEnalpha)
      .def("betaOrbEn", &orbEnbeta)
      .def("alphaFock", &fockalpha)
      .def("betaFock", &fockbeta);
}
