/**
 * @file PCMPotential.cpp
 *
 * @author Moritz Bensberg
 * @date Apr 28, 2020
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
#include "potentials/PCMPotential.h"
/* Include Serenity Internal Headers */
#include "data/OrbitalController.h"                           //getEigenvalues().
#include "data/grid/ElectrostaticPotentialOnGridController.h" //getDensityMatrixController()
#include "data/matrices/DensityMatrixController.h"            //getDensityMatrix()
#include "geometry/Geometry.h"                                //getNAtoms()
#include "geometry/MolecularSurfaceController.h"              //Point to atom mapping.
#include "integrals/wrappers/Libint.h"                        //compute1eInts()
#include "misc/SerenityError.h"                               //throw SerenityError
#include "misc/Timing.h"                                      //Timings.
#include "settings/PCMSettings.h"                             //Member .use.
#include "solvation/ContinuumModel.h"                         //getPCMCharges() and getActiveEnergy().

namespace Serenity {

template<Options::SCF_MODES SCFMode>
PCMPotential<SCFMode>::PCMPotential(const PCMSettings& pcmSettings, std::shared_ptr<BasisController> basisController,
                                    std::shared_ptr<const Geometry> geometry,
                                    std::shared_ptr<MolecularSurfaceController> molecularSurface,
                                    std::shared_ptr<MolecularSurfaceController> vdWmolecularSurface,
                                    std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>> activePotential,
                                    std::vector<std::shared_ptr<ElectrostaticPotentialOnGridController<SCFMode>>> environmentPotentials)
  : Potential<SCFMode>(basisController),
    _calculateG_cav(pcmSettings.cavityFormation && vdWmolecularSurface),
    _geometry(geometry),
    _vdWmolecularSurface(vdWmolecularSurface),
    _activePotential(activePotential),
    _environmentPotentials(environmentPotentials) {
  if (pcmSettings.use) {
    if (!molecularSurface || !activePotential)
      throw SerenityError("ERROR: PCM potential construction for a non-existing system.");
    _continuumModel =
        std::make_shared<ContinuumModel<SCFMode>>(pcmSettings, molecularSurface, activePotential, environmentPotentials);
    activePotential->getDensityMatrixController()->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
    for (auto envPotential : environmentPotentials)
      envPotential->getDensityMatrixController()->addSensitiveObject(ObjectSensitiveClass<DensityMatrix<SCFMode>>::_self);
  }
  else {
    _continuumModel = nullptr;
  }
  basisController->addSensitiveObject(ObjectSensitiveClass<Basis>::_self);
  auto& libint = Libint::getInstance();
  libint.keepEngines(LIBINT_OPERATOR::nuclear, 0, 2);
}

template<Options::SCF_MODES SCFMode>
PCMPotential<SCFMode>::~PCMPotential() {
  auto& libint = Libint::getInstance();
  libint.freeEngines(LIBINT_OPERATOR::nuclear, 0, 2);
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd PCMPotential<SCFMode>::integrateToFockMatrix(const GridPotential<RESTRICTED>& pcmCharges) {
  Timings::takeTime(" Tech. - PCM Charge Integration");
  Eigen::MatrixXd result = _activePotential->integrateFockMatrix(pcmCharges);
  Timings::timeTaken(" Tech. - PCM Charge Integration");
  return result;
}

template<Options::SCF_MODES SCFMode>
double PCMPotential<SCFMode>::getEnergy(const DensityMatrix<SCFMode>& P) {
  (void)P;
  double energy = 0.0;
  if (!_continuumModel)
    return energy;
  if (_calculateG_cav)
    energy += _vdWmolecularSurface->getCavityEnergy();
  energy += _continuumModel->getActivePCMEnergy();
  return energy;
}

template<Options::SCF_MODES SCFMode>
double PCMPotential<SCFMode>::getTotalEnergy() {
  double energy = 0.0;
  if (!_continuumModel)
    return 0.0;
  if (_calculateG_cav)
    energy += _vdWmolecularSurface->getCavityEnergy();
  energy += _continuumModel->getTotalPCMEnergy();
  return energy;
}

template<Options::SCF_MODES SCFMode>
FockMatrix<SCFMode>& PCMPotential<SCFMode>::getMatrix() {
  if (!_potential) {
    _potential.reset(new FockMatrix<SCFMode>(this->_basis));
    auto& pot = *_potential;
    if (_continuumModel) {
      Eigen::MatrixXd pcmFockContribution = integrateToFockMatrix(_continuumModel->getPCMCharges());
      for_spin(pot) {
        pot_spin.setZero();
        pot_spin += pcmFockContribution;
      };
    }
    else {
      for_spin(pot) {
        pot_spin.setZero();
      };
    }
  }
  return *_potential;
}

template<Options::SCF_MODES SCFMode>
Eigen::MatrixXd PCMPotential<SCFMode>::getGeomGradients() {
  const unsigned int nAtoms = _geometry->getNAtoms();
  Eigen::MatrixXd gradient(nAtoms, 3);
  gradient.setZero();
  if (!_continuumModel)
    return gradient;

  if (_continuumModel->getPCMSettings().solverType == Options::PCM_SOLVER_TYPES::CPCM && _environmentPotentials.size() == 0) {
    /*
     * This implementation currently works only for single systems, since
     * atom mapping (surface geometry != active geometry) would be necessary.
     * Furthermore, some terms would be missing which are associated with
     * the derivative of the 1/r operator for grid points of environment atoms.
     */
    auto molSurfaceController = _continuumModel->getMolecularSurfaceController();
    const std::vector<unsigned int>& pointWiseIndices = molSurfaceController->getPointToSphereMapping();
    const std::vector<std::pair<unsigned int, unsigned int>>& atomWiseIndices =
        molSurfaceController->getSphereToPointMapping();
    const auto& surfaceCoords = molSurfaceController->getGridPoints();
    const auto& atoms = _geometry->getAtoms();
    const unsigned int nPoints = molSurfaceController->getNGridPoints();
    const auto& pcmCharges = _continuumModel->getPCMCharges();
    // Terms from nabla V: Derivative of the nuclear potential.
    // No costly integrals are needed!
    for (unsigned int iAtom = 0; iAtom < nAtoms; ++iAtom) {
      const double z_I = atoms[iAtom]->getEffectiveCharge();
      Eigen::Vector3d r_I;
      r_I << atoms[iAtom]->getX(), atoms[iAtom]->getY(), atoms[iAtom]->getZ();
      for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint) {
        Eigen::Vector3d cont = r_I - surfaceCoords.col(iPoint);
        const double cont_norm = cont.norm();
        cont /= cont_norm * cont_norm * cont_norm;
        cont *= pcmCharges[iPoint] * z_I;
        // Contribution for delta(a,I)
        gradient.row(iAtom) -= cont.transpose();
        // Contribution for delta(s,a)
        const unsigned int aAtom = pointWiseIndices[iPoint];
        gradient.row(aAtom) += cont.transpose();
      } // for iPoint
    }   // for iAtom#
    // Terms from nabla V: Derivative of the potential from the electrons
    //(derivative of basis functions and operators are needed)!
    gradient += _activePotential->calculateGradient(pointWiseIndices, pcmCharges);
    // Derivative of the apparent surface charge interaction matrix.
    // Diagonal is 0.
    const double scaling = _continuumModel->getCPCMScaling();
    for (unsigned int aAtom = 0; aAtom < nAtoms; ++aAtom) {
      for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint) {
        for (unsigned int jPoint = atomWiseIndices[aAtom].first; jPoint < atomWiseIndices[aAtom].second; ++jPoint) {
          if (iPoint == jPoint)
            continue;
          const Eigen::Vector3d r_i = surfaceCoords.col(iPoint);
          const Eigen::Vector3d r_j = surfaceCoords.col(jPoint);
          const Eigen::Vector3d ri_rj = r_i - r_j;
          const double diffNorm = ri_rj.norm();
          const double prefactor = pcmCharges[iPoint] * pcmCharges[jPoint] / (scaling * diffNorm * diffNorm * diffNorm);
          gradient.row(aAtom) += prefactor * ri_rj.transpose();
        } // for jPoint
      }   // for iPoint
    }     // for aAtom
  }
  else {
    throw SerenityError("Analytical, geometrical gradients for IEFPCM and FDE-embedded CPCM are not implemented yet.");
  }
  return gradient;
}

template class PCMPotential<Options::SCF_MODES::RESTRICTED>;
template class PCMPotential<Options::SCF_MODES::UNRESTRICTED>;

} /* namespace Serenity */
