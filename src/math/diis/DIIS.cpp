/**
 * @file   DIIS.cpp
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
/* Include Class Header*/
#include "math/diis/DIIS.h"
/* Include Serenity Internal Headers */
#include "math/Matrix.h"
/* Include Std and External Headers */
#include <algorithm>
#include <cmath>
#include <memory>


namespace Serenity {
using namespace std;

DIIS::DIIS(
    const unsigned int maxStore,
    const double conditionNumberThreshold) :
          _maxStore(maxStore),
          _errorVectors(maxStore + 1),
          _targetVectors(maxStore + 1),
          _energies(maxStore + 1),
          _nStored(0),
          _conditionNumberThreshold(conditionNumberThreshold),
          _fockMatrixOnly(false){
}

DIIS::DIIS(
    const unsigned int maxStore,
    const double conditionNumberThreshold,
    const bool fockMatrixOnly,
    const bool diskMode) :
          _maxStore(maxStore),
          _errorVectors(maxStore + 1),
          _targetVectors(maxStore + 1),
          _energies(maxStore + 1),
          _nStored(0),
          _conditionNumberThreshold(conditionNumberThreshold),
          _fockMatrixOnly(fockMatrixOnly),
          _diskMode(diskMode){
  if(diskMode) {
    _errorDiskVectors.resize(maxStore+1);
    _targetDiskVectors.resize(maxStore+1);
  }
}


void DIIS::optimize(double energy,
    Eigen::Ref<Eigen::VectorXd> targetVector,
    const Eigen::Ref<const Eigen::VectorXd>& newErrorVector){
  if (_nStored > 0) {
    assert (_errorVectors[_nStored-1]->size()==int(newErrorVector.size()));
  }


  /*
   * Store the initial target vector with all the others and
   * do the same for the energy and the error vector.
   * New data is stored because it will be mutated throughout the cycles.
   */
  _targetVectors[_nStored].reset(new Eigen::VectorXd(targetVector));
  _energies[_nStored] = energy;
  _errorVectors[_nStored].reset(new Eigen::VectorXd(newErrorVector));
  unsigned int minEnergy = 0;
  for (unsigned int i=1; i<_nStored; ++i) if (_energies[i] < _energies[minEnergy]) minEnergy=i;
  ++_nStored;

  /*
   * Take care of the storage, delete vectors that are beyond the limit.
   */
  //bool shifted = false;
  if (_nStored > _maxStore) {
    //shifted = true;
    //shiftVectors();
    this->reinit();
    _targetVectors[_nStored].reset(new Eigen::VectorXd(targetVector));
    _energies[_nStored] = energy;
    _errorVectors[_nStored].reset(new Eigen::VectorXd(newErrorVector));
    ++_nStored;
  }

  /*
   * This is an improvement described in the manual of orca.
   * Manual for version 3.0, p.293
   *
   * In essence:
   * All but the 'best' matrix are disfavored by multiplying the diagonal element
   * by a factor of 1.05 this speeds up convergence, because the given 'bad 'vectors are
   * populated less.
   *
   * This trick is not general as it only works for a diis applied to a square matrices.
   * Therefore it is assumed that any vector with an inter sqrt of elements is a square matrix.
   */
  double testSqMat = sqrt(newErrorVector.size());
  if (_fockMatrixOnly) {
    assert(floor(testSqMat)==testSqMat);
    unsigned int nCol = (unsigned int)testSqMat;
    for (unsigned int i=0; i<nCol;++i){
      (*_errorVectors[_nStored-1])[i+i*nCol] *= 1.05;
    }
  }

  /*
   * Initialize the new B matrix
   * (  0  -1  -1  .. )
   * ( -1   ?   ?  .. )
   * ( -1   ?   ?  .. )
   * ( ..  ..  ..  .. )
   */
  auto newB = initNewB();

  /*
   * Compute error vector scalar products
   * (incl. second part of the trick described in ORCAs manual)
   */
  for (unsigned int i = 0; i < _nStored; ++i) {
    auto tmpI =  (*_errorVectors[i]);
    if (minEnergy==i && _fockMatrixOnly){
      unsigned int nCol = (unsigned int)testSqMat;
      for (unsigned int k=0; k<nCol;++k){
        tmpI[k+k*nCol] /= (1.05);
      }
    }
    for (unsigned int j = 0; j <=i; ++j) {
      auto tmpJ =  (*_errorVectors[j]);
      if (minEnergy==j && _fockMatrixOnly){
        unsigned int nCol = (unsigned int)testSqMat;
        for (unsigned int k=0; k<nCol;++k){
          tmpJ[k+k*nCol] /= (1.05);
        }
      }
      const double prod =  tmpI.cwiseProduct(tmpJ).sum();
      (*newB)(i+1,j+1) = prod;
      (*newB)(j+1,i+1) = prod;
    }
  }
  _B = std::move(newB);

  /*
   * Check whether B is ill-conditioned (iteratively, until its fine).
   */
  while (true) {
    /*
     * 1. Find max and min diagonal element (except first -> normalization).
     */
    double maxB = (*_B)(1, 1);
    double minB = (*_B)(1, 1);
    for (unsigned int i = 1; i < _nStored; ++i) {
      if ((*_B)(i + 1, i + 1) > maxB) {
        maxB = (*_B)(i + 1, i + 1);
      }
      if ((*_B)(i + 1, i + 1) < minB) {
        minB = (*_B)(i + 1, i + 1);
      }
    }
    /*
     * 2. Determine lower bound for the condition number.
     * 3. Check whether its above some threshold -> ill-conditioned
     */
    if (maxB / minB > _conditionNumberThreshold) {
      /*
       * 4. Shift up the vectors...
       */
      shiftVectors();
      /*
       * ...and delete first (i.e. second, first is just normalization)
       * row in B.
       */
      newB = initNewB();
      for (unsigned int i = 1; i < _nStored + 1; ++i) {
        for (unsigned int j = 1; j <= i; ++j) {
          (*newB)(i, j) = (*_B)(i + 1, j + 1);
        }
      }
      _B = std::move(newB);
    } else {
      /*
       * Matrix is well-conditioned.
       */
      break;
    }
  }

  /*
   * B is set up. Solve the linear equation system
   */
  if (_nStored > 1) {
    Eigen::VectorXd coefficients(_nStored + 1);
    Eigen::VectorXd rhs(_nStored + 1);
    rhs.fill(-1.0);
    coefficients = (*_B).jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
    /*
     * Construct a new optimized vector
     */
    targetVector.fill(0.0);
    for (unsigned int i = 0; i < _nStored; ++i) {
      for (unsigned int j = 0; j < targetVector.size(); ++j) {
        targetVector[j] += ((*_targetVectors[i])[j] * coefficients[i + 1]);
      }
    }
  }
  /*
   * Done.
   */
}

void DIIS::optimize(double energy,
    VectorOnDiskStorageController& targetVector,
    VectorOnDiskStorageController& newErrorVector) {
  _cycle++;
  if (_nStored > 0) {
    assert (_errorDiskVectors[_nStored-1]->size()== newErrorVector.size());
  }
  assert(_diskMode == true);

  /*
   * Store the initial target vector with all the others and
   * do the same for the energy and the error vector.
   * New data is stored because it will be mutated throughout the cycles.
   */
  auto nameTag = "DIISStored_"+ std::to_string(_cycle);
  _targetDiskVectors[_nStored].reset(new VectorOnDiskStorageController(targetVector, nameTag+targetVector.getHDF5FileName()));
  _energies[_nStored] = energy;
  _errorDiskVectors[_nStored].reset(new VectorOnDiskStorageController(newErrorVector, nameTag+newErrorVector.getHDF5FileName()));
  //	unsigned int minEnergy = 0;
  //	for (unsigned int i=1; i<_nStored; ++i) if (_energies[i] < _energies[minEnergy]) minEnergy=i;
  ++_nStored;

  /*
   * Take care of the storage, delete vectors that are beyond the limit.
   */
  //bool shifted = false;
  if (_nStored > _maxStore) {
    //shifted = true;
    //shiftVectors();
    this->reinit();
    _targetDiskVectors[_nStored].reset(new VectorOnDiskStorageController(targetVector, nameTag+targetVector.getHDF5FileName()));
    _energies[_nStored] = energy;
    _errorDiskVectors[_nStored].reset(new VectorOnDiskStorageController(newErrorVector, nameTag+newErrorVector.getHDF5FileName()));
    ++_nStored;
  }

  /*
   * Initialize the new B matrix
   * (  0  -1  -1  .. )
   * ( -1   ?   ?  .. )
   * ( -1   ?   ?  .. )
   * ( ..  ..  ..  .. )
   */
  auto newB = initNewB();

  /*
   * Compute error vector scalar products
   */
  for (unsigned int i = 0; i < _nStored; ++i) {
    auto& tmpI =  (*_errorDiskVectors[i]);
    for (unsigned int j = 0; j <=i; ++j) {
      auto& tmpJ =  (*_errorDiskVectors[j]);
      const double prod =  tmpI*tmpJ;
      (*newB)(i+1,j+1) = prod;
      (*newB)(j+1,i+1) = prod;
    }
  }
  _B = std::move(newB);

  /*
   * Check whether B is ill-conditioned (iteratively, until its fine).
   */
  while (true) {
    /*
     * 1. Find max and min diagonal element (except first -> normalization).
     */
    double maxB = (*_B)(1, 1);
    double minB = (*_B)(1, 1);
    for (unsigned int i = 1; i < _nStored; ++i) {
      if ((*_B)(i + 1, i + 1) > maxB) {
        maxB = (*_B)(i + 1, i + 1);
      }
      if ((*_B)(i + 1, i + 1) < minB) {
        minB = (*_B)(i + 1, i + 1);
      }
    }
    /*
     * 2. Determine lower bound for the condition number.
     * 3. Check whether its above some threshold -> ill-conditioned
     */
    if (maxB / minB > _conditionNumberThreshold) {
      /*
       * 4. Shift up the vectors...
       */
      shiftVectors();
      /*
       * ...and delete first (i.e. second, first is just normalization)
       * row in B.
       */
      newB = initNewB();
      for (unsigned int i = 1; i < _nStored + 1; ++i) {
        for (unsigned int j = 1; j <= i; ++j) {
          (*newB)(i, j) = (*_B)(i + 1, j + 1);
        }
      }
      _B = std::move(newB);
    } else {
      /*
       * Matrix is well-conditioned.
       */
      break;
    }
  }

  /*
   * B is set up. Solve the linear equation system
   */
  if (_nStored > 1) {
    Eigen::VectorXd coefficients(_nStored + 1);
    Eigen::VectorXd rhs(_nStored + 1);
    rhs.fill(-1.0);
    coefficients = (*_B).jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);
    /*
     * Construct a new optimized vector
     */
    const auto labelList = targetVector.getLabelList();
    for(const auto& label : labelList) {
      auto optSegment = std::make_shared<Eigen::VectorXd>(*targetVector.getVectorSegment(label));
      optSegment->setZero();
      for (unsigned int i = 0; i < _nStored; ++i) {
        *optSegment += (coefficients[i + 1] * (*_targetDiskVectors[i]->getVectorSegment(label)));
      }
      targetVector.storeVectorSegment(optSegment,label);
    }
  }
  /*
   * Done.
   */
}


void DIIS::reinit() {
  /*
   * Just to make sure: in case of reinitialization kill leftovers.
   */
  cleanUp();
  _B.reset(new Matrix<double>(0,0));
  _targetVectors.resize(_maxStore + 1);
  _targetDiskVectors.resize(_maxStore + 1);
  _errorVectors.resize(_maxStore + 1);
  _errorDiskVectors.resize(_maxStore + 1);
  _energies.resize(_maxStore + 1);
};

void DIIS::shiftVectors() {
  --_nStored;
  /*
   * Shift up the other data in the list
   */
  for (unsigned int i = 0; i < _nStored; i++) {
    if (_diskMode) {
      _errorDiskVectors[i] = move(_errorDiskVectors[i + 1]);
      _targetDiskVectors[i] = move(_targetDiskVectors[i + 1]);
    } else {
      _errorVectors[i] = move(_errorVectors[i + 1]);
      _targetVectors[i] = move(_targetVectors[i + 1]);
    }
    _energies[i] = _energies[i + 1];
  }
  _energies[_nStored] = 0.0;
}

unique_ptr<Matrix<double> > DIIS::initNewB() {
  auto newB = unique_ptr<Matrix<double> >(new Matrix<double>(_nStored + 1,_nStored + 1));
  newB->fill(0.0);
  for (unsigned int i = 1; i < _nStored + 1; i++) {
    (*newB)(0, i) = -1.0;
    (*newB)(i, 0) = -1.0;
  }
  return newB;
}

void DIIS::cleanUp() {
  _B.reset(nullptr);
  if (_nStored > 0) {
    _targetVectors.resize(0);
    _targetDiskVectors.resize(0);
    _errorVectors.resize(0);
    _errorDiskVectors.resize(0);
    _energies.resize(0);
    _nStored = 0;
  }
}

} /* namespace Serenity */
