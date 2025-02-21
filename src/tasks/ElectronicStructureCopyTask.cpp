/**
 * @file ElectronicStructureCopyTask.cpp
 *
 * @author Moritz Bensberg
 * @date Feb 13, 2020
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
#include "tasks/ElectronicStructureCopyTask.h"
/* Include Serenity Internal Headers */
#include "basis/AtomCenteredBasisControllerFactory.h" //Target basis construction before orbital projection.
#include "basis/Basis.h"                              //Shell-wise rotations.
#include "basis/BasisController.h"                    //Access to the basis.
#include "basis/SphericalHarmonicsRotations.h"        //Rotation of spherical harmonics.
#include "data/ElectronicStructure.h"                 //Constructor.
#include "data/OrbitalController.h"                   //Orbital rotations.
#include "geometry/Atom.h"                            //Cross checks for atoms and frame construction.
#include "geometry/Geometry.h"                        //Cross checks for atoms and frame construction.
#include "integrals/OneElectronIntegralController.h"  //A--A overlap matrix for projection.
#include "integrals/wrappers/Libint.h"                //A--B overlap matrix for projection.
#include "io/FormattedOutput.h"                       //Captions.
#include "io/FormattedOutputStream.h"                 //Filtered output streams.
#include "math/linearAlgebra/MatrixFunctions.h"       //Pseudo inverse square root for orthogonalisation.
#include "misc/SerenityError.h"                       //Error messages.
#include "settings/Settings.h"                        //Settings definition.
#include "system/SystemController.h"                  //Access to system properties.
/* Include Std and External Headers */
#include <math.h> //acos
#include <limits> //std::numeric_limits<double>::infinity()

namespace Serenity {

template<Options::SCF_MODES SCFMode>
ElectronicStructureCopyTask<SCFMode>::ElectronicStructureCopyTask(std::shared_ptr<SystemController> sourceSystem,
                                                                  std::vector<std::shared_ptr<SystemController>> targetSystems)
  : _sourceSystem(sourceSystem), _targetSystems(targetSystems) {
}

template<Options::SCF_MODES SCFMode>
std::vector<Eigen::Vector3d> ElectronicStructureCopyTask<SCFMode>::getInternalFrame(std::shared_ptr<Geometry> geometry) {
  auto atoms = geometry->getAtoms();
  unsigned int nAtoms = atoms.size();
  // Get coordinates of first point.
  unsigned int atomIndex0 = settings.atomFrameIndices.size() < 3 ? 0 : settings.atomFrameIndices[0];
  unsigned int atomIndex1 = settings.atomFrameIndices.size() < 3 ? 1 : settings.atomFrameIndices[1];
  unsigned int atomIndex2 = settings.atomFrameIndices.size() < 3 ? 2 : settings.atomFrameIndices[2];

  Eigen::Vector3d vecA = atoms[atomIndex0]->coords();
  Eigen::Vector3d e_x;
  Eigen::Vector3d e_y;
  Eigen::Vector3d e_z;
  if (nAtoms > 1) {
    // More than one atom --> x-axis will be the connection between the atoms.
    Eigen::Vector3d vecB = atoms[atomIndex1]->coords();
    Eigen::Vector3d vecAB = vecB - vecA;
    double normAB = vecAB.norm();
    if (normAB < 1e-4)
      throw SerenityError("ERROR: The first two atoms have identical coordinates! Check your geometry!");
    // Calculate x-unit vector from AB axis.
    e_x = 1.0 / normAB * vecAB;
    double normBCortho = 0;
    Eigen::Vector3d vecBCortho;
    // Try to construct the y-axis as the connection to the next atom that is not linear to the first two atoms.
    bool linear = true;
    std::vector<unsigned int> atomIndices;
    if (nAtoms > atomIndex2)
      atomIndices.insert(atomIndices.begin(), atomIndex2);
    for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom)
      if (iAtom != atomIndex2 and iAtom != atomIndex1 and iAtom != atomIndex0)
        atomIndices.push_back(iAtom);
    for (const auto iAtom : atomIndices) {
      Eigen::Vector3d vecC = atoms[iAtom]->coords();
      Eigen::Vector3d vecBC = vecC - vecB;
      Eigen::Vector3d vecBCNormed = vecBC * 1.0 / vecBC.norm();
      double linearDependentCheck = std::fabs(vecBCNormed.transpose() * e_x);
      if (std::fabs(linearDependentCheck - 1.0) > 1e-5) {
        // BC vector is linear independent of e_x. Orthogonalize to e_x, normalize and break the loop.
        vecBCortho = vecBC - (e_x.transpose() * vecBC) * e_x;
        normBCortho = vecBCortho.norm();
        linear = false;
        break;
      }
    } // for iAtom
    // If there is no third atom which has a connection to the second atom that is linear independent of e_x, the
    // molecule has to be linear. The e_y-axis will be chosen to be the to e_x orthogonalized Cartesian x- or y-axis.
    if (linear) {
      Eigen::Vector3d u;
      u << 1, 0, 0;
      vecBCortho = u - (e_x.transpose() * u) * e_x;
      normBCortho = vecBCortho.norm();
      if (normBCortho < 1e-4) {
        u << 0, 1, 0;
        vecBCortho = u - (e_x.transpose() * u) * e_x;
      }
    } // if linear
    // The e_z axis is constructed from the cross product in order to ensure that a right handed coordinate system
    // is constructed.
    e_y = 1.0 / normBCortho * vecBCortho;
    e_z = e_x.cross(e_y);
    e_z *= 1.0 / e_z.norm();
  }
  else {
    // If there is only one atom. The frame will be the Cartesian frame.
    e_x << 1, 0, 0;
    e_y << 0, 1, 0;
    e_z << 0, 0, 1;
  } // else nAtoms > 1
  std::vector<Eigen::Vector3d> frame(3);
  frame[0] = e_x;
  frame[1] = e_y;
  frame[2] = e_z;
  return frame;
}

template<Options::SCF_MODES SCFMode>
Eigen::Vector3d ElectronicStructureCopyTask<SCFMode>::getEulerAngles(const std::vector<Eigen::Vector3d>& sourceFrame,
                                                                     const std::vector<Eigen::Vector3d>& targetFrame) {
  const Eigen::Vector3d& e_x = sourceFrame[0];
  const Eigen::Vector3d& e_y = sourceFrame[1];
  const Eigen::Vector3d& e_z = sourceFrame[2];
  const Eigen::Vector3d& u_x = targetFrame[0];
  const Eigen::Vector3d& u_y = targetFrame[1];
  const Eigen::Vector3d& u_z = targetFrame[2];
  double eZ_uY = e_z.transpose() * u_y;
  double eZ_uX = e_z.transpose() * u_x;
  bool orthoYz = std::fabs(eZ_uY) < 1e-6;
  bool orthoXz = std::fabs(eZ_uX) < 1e-6;
  Eigen::Vector3d xyXYIntersection;
  if (orthoYz && orthoXz) {
    xyXYIntersection = e_y;
  }
  else {
    if (orthoXz) {
      xyXYIntersection = u_x;
    }
    else {
      xyXYIntersection = u_y - (eZ_uY / eZ_uX) * u_x;
    }
    xyXYIntersection *= 1.0 / xyXYIntersection.norm();
  }
  double sign = (e_x.transpose() * xyXYIntersection > 0) ? -1.0 : 1.0;
  double alpha = sign * acos(e_y.transpose() * xyXYIntersection); // Sign?
  sign = (e_z.cross(u_z).transpose() * xyXYIntersection > 0) ? 1.0 : -1.0;
  double beta = sign * acos(e_z.transpose() * u_z);
  sign = (u_x.transpose() * xyXYIntersection < 0) ? -1.0 : 1.0;
  double gamma = sign * acos(xyXYIntersection.transpose() * u_y);
  Eigen::Vector3d angles;
  angles << alpha, beta, gamma;
  return angles;
}

template<Options::SCF_MODES SCFMode>
void ElectronicStructureCopyTask<SCFMode>::checkSystems() {
  auto sourceBasisSet = _sourceSystem->getBasisController();
  auto sourceGeometry = _sourceSystem->getGeometry();
  auto sourceAtoms = sourceGeometry->getAtoms();
  if (!sourceBasisSet->isPureSpherical())
    throw SerenityError("ERROR: Copying and rotating systems is only supported for purely spherical basis sets.");
  for (auto& targetSystem : _targetSystems) {
    auto targetBasisSet = _sourceSystem->getBasisController();
    if (!targetBasisSet->isPureSpherical())
      throw SerenityError("ERROR: Copying and rotating systems is only supported for purely spherical basis sets.");
    auto targetGeometry = targetSystem->getGeometry();
    auto targetAtoms = targetGeometry->getAtoms();
    if (targetAtoms.size() != sourceAtoms.size())
      throw SerenityError((std::string) "ERROR: The number of atoms is not identical in source and target system " +
                          targetSystem->getSystemName());
    for (unsigned int iAtom = 0; iAtom < targetAtoms.size(); ++iAtom) {
      auto atomTypeSource = sourceAtoms[iAtom]->getAtomType();
      auto atomTypeTarget = targetAtoms[iAtom]->getAtomType();
      if (atomTypeSource != atomTypeTarget)
        throw SerenityError(
            (std::string) "ERROR: Different atom-types detected for the source system and target system " +
            targetSystem->getSystemName() + " atom index " + iAtom);
    } // for iAtom
  }   // for targetSystem
  // Check if the first two atom indices are identical.
  // For the third atom index linear dependency checks are used any way, so an identical index
  // should not result into difficulties.
  if (settings.atomFrameIndices.size() > 1) {
    if (settings.atomFrameIndices[0] == settings.atomFrameIndices[1])
      throw SerenityError(
          "The atom indices that are used for the internal frame construction are not allowed to be identical!");
  }
}
template<Options::SCF_MODES SCFMode>
CoefficientMatrix<SCFMode>
ElectronicStructureCopyTask<SCFMode>::rotateMatrixIntoCartesianFrame(const std::vector<Eigen::Vector3d>& sourceFrame,
                                                                     const CoefficientMatrix<SCFMode>& coefficients) {
  CoefficientMatrix<SCFMode> newCoefficients(coefficients.getBasisController());
  std::vector<Eigen::Vector3d> cartesianFrame(3);
  cartesianFrame[0] << 1, 0, 0;
  cartesianFrame[1] << 0, 1, 0;
  cartesianFrame[2] << 0, 0, 1;
  const Eigen::Vector3d eulerAngles = getEulerAngles(cartesianFrame, sourceFrame);
  auto basis = coefficients.getBasisController()->getBasis();
  unsigned int nOrbitals = newCoefficients.cols();
  unsigned int startRow = 0;
  for (const auto& shell : basis) {
    unsigned int angularMomentum = shell->getAngularMomentum();
    unsigned int nContracted = shell->getNContracted();
    Eigen::MatrixXd rotationMatrix = SphericalHarmonicsRotations::getTransformationMatrix(
        angularMomentum, -eulerAngles[2], -eulerAngles[1], -eulerAngles[0]);
    for_spin(newCoefficients, coefficients) {
      newCoefficients_spin.block(startRow, 0, nContracted, nOrbitals) =
          rotationMatrix * coefficients_spin.block(startRow, 0, nContracted, nOrbitals);
    };
    startRow += nContracted;
  } // for shell
  return newCoefficients;
}

template<Options::SCF_MODES SCFMode>
CoefficientMatrix<SCFMode>
ElectronicStructureCopyTask<SCFMode>::rotateMatrixIntoTargetFrame(const std::vector<Eigen::Vector3d>& targetFrame,
                                                                  const CoefficientMatrix<SCFMode>& coefficients,
                                                                  std::shared_ptr<BasisController> targetBasis) {
  CoefficientMatrix<SCFMode> newCoefficients(targetBasis);
  std::vector<Eigen::Vector3d> cartesianFrame(3);
  cartesianFrame[0] << 1, 0, 0;
  cartesianFrame[1] << 0, 1, 0;
  cartesianFrame[2] << 0, 0, 1;
  const Eigen::Vector3d eulerAngles = getEulerAngles(cartesianFrame, targetFrame);
  auto basis = coefficients.getBasisController()->getBasis();
  unsigned int nOrbitals = newCoefficients.cols();
  unsigned int startRow = 0;
  for (const auto& shell : basis) {
    unsigned int angularMomentum = shell->getAngularMomentum();
    unsigned int nContracted = shell->getNContracted();
    Eigen::MatrixXd rotationMatrix = SphericalHarmonicsRotations::getTransformationMatrix(angularMomentum, eulerAngles[0],
                                                                                          eulerAngles[1], eulerAngles[2]);
    for_spin(newCoefficients, coefficients) {
      newCoefficients_spin.block(startRow, 0, nContracted, nOrbitals) =
          rotationMatrix * coefficients_spin.block(startRow, 0, nContracted, nOrbitals);
    };
    startRow += nContracted;
  } // for shell
  return newCoefficients;
}
template<Options::SCF_MODES SCFMode>
void ElectronicStructureCopyTask<SCFMode>::loewdingOrthogonalization(CoefficientMatrix<SCFMode>& coefficients,
                                                                     const MatrixInBasis<RESTRICTED>& S) {
  for_spin(coefficients) {
    Eigen::MatrixXd S_MO = coefficients_spin.transpose() * S * coefficients_spin;
    coefficients_spin = coefficients_spin * pseudoInversSqrt_Sym(S_MO);
  };
}

template<Options::SCF_MODES SCFMode>
void ElectronicStructureCopyTask<SCFMode>::run() {
  for (auto& sys : _targetSystems) {
    if (sys->getSCFMode() != SCFMode)
      throw SerenityError("ERROR: The target system '" + sys->getSystemName() + "' needs the same SCFMode as the source system.");
  }
  if (GLOBAL_PRINT_LEVEL > Options::GLOBAL_PRINT_LEVELS::MINIMUM)
    printSubSectionTitle("Copying Electronic Structures");
  checkSystems();
  OutputControl::nOut << "  Source system:\n"
                      << "    " << _sourceSystem->getSystemName() << std::endl;
  OutputControl::nOut << "  Copying target systems:" << std::endl;
  const CoefficientMatrix<SCFMode>& sourceCoefficients =
      _sourceSystem->getActiveOrbitalController<SCFMode>()->getCoefficients();
  const auto& sourceEigenvalues = _sourceSystem->getActiveOrbitalController<SCFMode>()->getEigenvalues();
  const auto& sourceCoreOrbitals = _sourceSystem->getActiveOrbitalController<SCFMode>()->getOrbitalFlags();
  // Construct internal frame for the source system.
  const std::vector<Eigen::Vector3d> sourceFrame = getInternalFrame(_sourceSystem->getGeometry());
  // Rotate the Electronic structure of the source system into the Cartesian frame by rotating the Cartesian axes.
  const CoefficientMatrix<SCFMode> sourceCoeffInCartesian = rotateMatrixIntoCartesianFrame(sourceFrame, sourceCoefficients);
  auto sourceBasisController = _sourceSystem->getBasisController();
  const int sourceSpin = _sourceSystem->getSpin();
  const int sourceCharge = _sourceSystem->getCharge();
  std::string sourceBasisLabel = _sourceSystem->getSettings().basis.label;
  for (auto& targetSystem : _targetSystems) {
    OutputControl::nOut << "    " << targetSystem->getSystemName() << "     ...";
    OutputControl::nOut.flush();
    if (settings.copyCharges) {
      targetSystem->setSpin(sourceSpin);
      targetSystem->setCharge(sourceCharge);
    }
    std::shared_ptr<BasisController> targetBasisController = targetSystem->getBasisController();
    const unsigned int nBasisFunctions = targetBasisController->getNBasisFunctions();
    auto rotatedCoefficientsPtr =
        std::unique_ptr<CoefficientMatrix<SCFMode>>(new CoefficientMatrix<SCFMode>(targetBasisController));
    auto eigenvaluesPtr = std::unique_ptr<SpinPolarizedData<SCFMode, Eigen::VectorXd>>(
        new SpinPolarizedData<SCFMode, Eigen::VectorXd>(Eigen::VectorXd::Zero(nBasisFunctions)));
    auto coreOrbitalsPtr =
        std::make_unique<SpinPolarizedData<SCFMode, Eigen::VectorXi>>(Eigen::VectorXi::Zero(nBasisFunctions));
    // Construct the internal frame of the target system.
    const std::vector<Eigen::Vector3d> targetFrame = getInternalFrame(targetSystem->getGeometry());
    // Rotate the source coefficients, expressed in the Cartesian frame to the target frame by rotating the Cartesian
    // axes.
    const CoefficientMatrix<SCFMode> rotatedCoefficients =
        rotateMatrixIntoTargetFrame(targetFrame, sourceCoeffInCartesian, targetBasisController);

    const Settings& targetSystemSettings = targetSystem->getSettings();
    if (targetSystemSettings.basis.label == sourceBasisLabel) {
      *rotatedCoefficientsPtr = rotateMatrixIntoTargetFrame(targetFrame, sourceCoeffInCartesian, targetBasisController);
      *eigenvaluesPtr = sourceEigenvalues;
      *coreOrbitalsPtr = sourceCoreOrbitals;
    }
    else {
      // Project the orbitals.
      unsigned int nSourceBasisFunctions = sourceBasisController->getNBasisFunctions();
      unsigned int nParsed = std::min(nBasisFunctions, nSourceBasisFunctions);
      const CoefficientMatrix<SCFMode> rotatedCoefficients =
          rotateMatrixIntoTargetFrame(targetFrame, sourceCoeffInCartesian, sourceBasisController);
      CoefficientMatrix<SCFMode>& parseCoefficients = *rotatedCoefficientsPtr;
      SpinPolarizedData<SCFMode, Eigen::VectorXd>& parseEigenvalues = *eigenvaluesPtr;
      SpinPolarizedData<SCFMode, Eigen::VectorXi>& parseCoreOrbitals = *coreOrbitalsPtr;
      auto translatedSourceBasis = AtomCenteredBasisControllerFactory::produce(
          targetSystem->getGeometry(), targetSystemSettings.basis.basisLibPath,
          targetSystemSettings.basis.makeSphericalBasis, false, targetSystemSettings.basis.firstECP, sourceBasisLabel);
      auto& libint = Libint::getInstance();
      Eigen::MatrixXd targetOverlap = targetSystem->getOneElectronIntegralController()->getOverlapIntegrals(); // libint.compute1eInts(
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(targetOverlap, Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::MatrixXd overlapSourceTarget =
          libint.compute1eInts(LIBINT_OPERATOR::overlap, translatedSourceBasis, targetBasisController);
      Eigen::MatrixXd projection = svd.solve(overlapSourceTarget);
      for_spin(rotatedCoefficients, parseCoefficients, parseEigenvalues, sourceEigenvalues, parseCoreOrbitals,
               sourceCoreOrbitals) {
        parseEigenvalues_spin = Eigen::VectorXd::Constant(nBasisFunctions, std::numeric_limits<double>::infinity());
        parseCoefficients_spin.leftCols(nParsed) = projection * rotatedCoefficients_spin.leftCols(nParsed);
        parseEigenvalues_spin.segment(0, nParsed) = sourceEigenvalues_spin.segment(0, nParsed);
        parseCoreOrbitals_spin.segment(0, nParsed) = sourceCoreOrbitals_spin.segment(0, nParsed);
      };
    }
    if (settings.orthogonalize)
      loewdingOrthogonalization(*rotatedCoefficientsPtr,
                                targetSystem->getOneElectronIntegralController()->getOverlapIntegrals());
    // Construct the new electronic structure of the target system and write it to disk.
    auto newOrbitalController = std::make_shared<OrbitalController<SCFMode>>(
        std::move(rotatedCoefficientsPtr), targetBasisController, std::move(eigenvaluesPtr), std::move(coreOrbitalsPtr));
    auto targetElectronicStructure =
        std::make_shared<ElectronicStructure<SCFMode>>(newOrbitalController, targetSystem->getOneElectronIntegralController(),
                                                       targetSystem->template getNOccupiedOrbitals<SCFMode>());
    targetSystem->template setElectronicStructure<SCFMode>(targetElectronicStructure);
    targetElectronicStructure->toHDF5(targetSystem->getHDF5BaseName(), targetSystem->getSystemIdentifier());
    OutputControl::nOut << " done" << std::endl;
  } // for targetSystem
  OutputControl::nOut << std::endl;
}

template class ElectronicStructureCopyTask<Options::SCF_MODES::RESTRICTED>;
template class ElectronicStructureCopyTask<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
