/**
 * @file   LRSCFTask.h
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
#ifndef LRSCFTASK_H_
#define LRSCFTASK_H_
/* Include Serenity Internal Headers */
#include "settings/Reflection.h"
#include "data/SpinPolarizedData.h"
#include "tasks/Task.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <vector>


namespace Serenity {
/* Forward declarations */
class SystemController;

using namespace::Serenity::Reflection;
struct LRSCFTaskSettings {
  LRSCFTaskSettings():
    nEigen(1),
    maxCycles(100),
    maxSubspaceDimension(1.0e6),
    convergenceCriterion(1.0e-5),
    rpa(false),
    tda(false),
    kss(false),
    noNaddKernel(false),
    superSystemGrid(false),
    mullikenpop(false),
    writeAOVectors(false),
    multiplicity(Options::MULTIPLICITY::SINGLET),
    func(Options::FUNCTIONALS::LDA),
    naddKinFunc(Options::FUNCTIONALS::TF),
    naddXCFunc(Options::FUNCTIONALS::LDA){}
  REFLECTABLE(
    (unsigned int) nEigen,
    (unsigned int) maxCycles,
    (unsigned int) maxSubspaceDimension,
    (double) convergenceCriterion,
    (bool) rpa,
    (bool) tda,
    (bool) kss,
    (bool) spa,
    (bool) noNaddKernel,
    (bool) superSystemGrid,
    (bool) mullikenpop,
    (bool) writeAOVectors,
    (Options::MULTIPLICITY) multiplicity,
    (Options::FUNCTIONALS) func,
    (Options::FUNCTIONALS) naddKinFunc,
    (Options::FUNCTIONALS) naddXCFunc
  )
};

/**
 * @class  LRSCFTask LRSCF.h
 * @brief  Perform a LRSCF calculation
 */
template<Options::SCF_MODES T>
class LRSCFTask: public Task {
public:
  /**
   * @brief Constructor.
   * @param activeSystem The main (active) system.
   * @param environmentSystems Environment systems.
   */
  LRSCFTask(
      std::shared_ptr<SystemController> activeSystem,
      std::vector<std::shared_ptr<SystemController> > environmentSystems);
  /**
   * @brief Default destructor.
   */
  virtual ~LRSCFTask() = default;
  /**
   * @see Task.h
   */
  void run();

  LRSCFTaskSettings settings;
private:
  std::shared_ptr<SystemController> _activeSystem;
  std::vector<std::shared_ptr<SystemController> > _environmentSystems;

  //The SCF type, i.e. HF or DFT
  Options::ELECTRONIC_STRUCTURE_THEORIES _scfType;

  //The number of electrons
  SpinPolarizedData<T,unsigned int> _nElectrons;

  //The number of occupied orbitals
  SpinPolarizedData<T,unsigned int> _nOccupied;

  //The number of virtual orbitals
  SpinPolarizedData<T,unsigned int> _nVirtual;

  //The orbital energies
  SpinPolarizedData<T,Eigen::VectorXd > _orbitalEnergies;

  //The dimension of the eigenvalue problem
  int _nDimension;

  //True if hybrid functional is used
  bool _hybrid;


  void printSettings();
  Eigen::VectorXd estimateDiagonal();

};

} /* namespace Serenity */

#endif  /* LRSCFTASK_H_*/
