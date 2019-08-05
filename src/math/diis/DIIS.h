/**
 * @file   DIIS.h
 *
 * @date   Nov 18, 2013
 * @author Thomas Dresselhaus
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
#ifndef DIIS_H_
#define DIIS_H_
/* Include Serenity Internal Headers */
#include "math/optimizer/Optimizer.h"
/* Experimental FaT convergence acceleration */
#include "misc/VectorOnDiskStorageController.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <vector>


namespace Serenity {
/* Forward calculations */
template<class T>class Matrix;
/**
 * @class DIIS DIIS.h
 * @brief An implementation of the DIIS by Pulay.
 *
 * According to:
 * Journal of Computational Chemistry 3 (4): 556–560; doi:10.1002/jcc.540030413
 */
class DIIS {
public:
  /**
   * @param maxStore                       How many error vectors and target matrices are stored at max.
   *                                       If reached (at least) the oldest ones are deleted in the
   *                                       following cycle.
   * @param conditionNumberThreshold       If the B-Matrix is ill-conditioned numerical instabilities can
   *                                       lead to huge problems. Thus old entries in the B-Matrix are
   *                                       thrown out if it is ill-conditioned, i.e. the condition number
   *                                       is above this threshold. Easily said I think this means that
   *                                       far too different (i.e. too old) target matrices are not
   *                                       considered while constructing new ones.
   */
  DIIS(
    unsigned int maxStore,
    double conditionNumberThreshold);

  /**
   * @param maxStore                       How many error vectors and target matrices are stored at max.
   *                                       If reached (at least) the oldest ones are deleted in the
   *                                       following cycle.
   * @param conditionNumberThreshold       If the B-Matrix is ill-conditioned numerical instabilities can
   *                                       lead to huge problems. Thus old entries in the B-Matrix are
   *                                       thrown out if it is ill-conditioned, i.e. the condition number
   *                                       is above this threshold. Easily said I think this means that
   *                                       far too different (i.e. too old) target matrices are not
   *                                       considered while constructing new ones.
   * @param fockMatrixOnly                 A bool specifying that only Fock matrices or other symmetric
   *                                       matrices will be handed over, here the DIIS will further bias
   *                                       towards the best available matrix and therefore further
   *                                       accelerate convergence.
   */
  DIIS(
    unsigned int maxStore,
    double conditionNumberThreshold,
    bool fockMatrixOnly,
	bool diskMode = false);
  virtual ~DIIS() = default;

  /**
   * @brief Perform a DIIS step
   *
   * @param value       energy for the parameters corresponding to the gradients
   * @param parameters  the target vector of the current cycle. It will be mutated to be
   *                    hopefully closer to the optimum as the incoming vector.
   * @param gradients   the error vector of the current optimization cycle. It is some measure
   *                    for the gradient.
   */
  void optimize(double value,
      Eigen::Ref<Eigen::VectorXd> parameters,
      const Eigen::Ref<const Eigen::VectorXd>& gradients);

  /**
   * @brief Perform a DIIS step
   *
   * @param value       energy for the parameters corresponding to the gradients
   * @param parameters  the target vector of the current cycle. It will be mutated to be
   *                    hopefully closer to the optimum as the incoming vector.
   * @param gradients   the error vector of the current optimization cycle. It is some measure
   *                    for the gradient.
   */
  void optimize(double value,
                Eigen::MatrixXd& parameters,
                const Eigen::MatrixXd& gradients){
    this->optimize(value,
        Eigen::Map<Eigen::VectorXd>(parameters.data(), parameters.cols()*parameters.rows()),
        Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.cols()*gradients.rows()));
  }

  /**
   * @brief Perform a DIIS step
   *
   * @param value       energy for the parameters corresponding to the gradients
   * @param parameters  the target vector of the current cycle. It will be mutated to be
   *                    hopefully closer to the optimum as the incoming vector.
   * @param gradients   the error vector of the current optimization cycle. It is some measure
   *                    for the gradient.
   */
  void optimize(double value,
                std::vector<double>& parameters,
                const std::vector<double>& gradients){
    this->optimize(value,Eigen::Map<Eigen::VectorXd>(parameters.data(), parameters.size()),
                         Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.size()));
  }
  /**
   * @brief Perform a DIIS step
   *
   * @param energy           energy for the parameters corresponding to the gradients
   * @param targetVector     the target vector of the current cycle. It will be mutated to be
   *                         hopefully closer to the optimum as the incoming vector.
   * @param newErrorVector   the error vector of the current optimization cycle. It is some measure
   *                         for the gradient.
   */
  void optimize(double energy,
		  VectorOnDiskStorageController& targetVector,
		  VectorOnDiskStorageController& newErrorVector);

  /**
   * (Re-)initializes this object. All stored data is erased.
   */
  void reinit();

private:
  /*
   * How many error vectors and target matrices are stored at max
   */
  unsigned int _maxStore;

  /*
   * A list of old error vectors indicating how far the algorithm is
   * from convergence.
   */
  std::vector<std::unique_ptr< Eigen::VectorXd > > _errorVectors;
  std::vector<std::unique_ptr< VectorOnDiskStorageController> > _errorDiskVectors;

  /*
   * A list of old target vectors. Optimized vectors are constructed
   * with respect to these. The aim of the algorithm is to find the
   * ideal targetVector which reduces the errorVector to zero.
   */
  std::vector<std::unique_ptr< Eigen::VectorXd> > _targetVectors;
  std::vector<std::unique_ptr< VectorOnDiskStorageController> > _targetDiskVectors;

  /*
   * A list of the calculated energies for the corresponding _errorVectors
   * and _targetVectors
   */
  Eigen::VectorXd _energies;

  /*
   * How many error vectors and target matrices are currently stored
   */
  unsigned int _nStored;

  /*
   * Holds the scalar products of the error vectors. It is used in an
   * eigenvalue equation to determine the optimum amounts of how much
   * certain target matrices should get into the new one.
   */
  std::unique_ptr<Matrix<double> > _B;

  /*
   * If the B-Matrix is ill-conditioned numerical instabilities can lead
   * to huge problems. Thus old entries in the B-Matrix are thrown out
   * if it is ill-conditioned, i.e. the condition number is above this
   * threshold. Easily said I think this means that far too different
   * (i.e. too old) target vectors are not considered while constructing
   * new ones.
   */
  const double _conditionNumberThreshold;

  const bool _fockMatrixOnly;

  /*
   * A flag for the use of VectorOnDiskStorageController instead of
   * conventional Eigen vectors.
   */
  bool _diskMode;

  /*
   * Erase the oldest stuff and shift up the other entries in the vectors.
   */
  void shiftVectors();

  /*
   * Initialize a new B matrix including the normalization condition in
   * the first row.
   */
  std::unique_ptr<Matrix<double> > initNewB();

  void cleanUp();

  unsigned int _cycle=0;
};

} /* namespace Serenity */
#endif /* DIIS_H_ */
