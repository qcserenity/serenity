/**
 * @file KernelSigmaVector.h
 *
 * @date Oct 09, 2017
 * @author Michael Boeckers
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

#ifndef POSTHF_RESPONSE_KERNELSIGMAVECTOR_H_
#define POSTHF_RESPONSE_KERNELSIGMAVECTOR_H_

/* Include Serenity Internal Headers */
#include "data/grid/BasisFunctionOnGridController.h"
#include "postHF/LRSCF/Kernel.h"
#include "postHF/LRSCF/SigmaVector/SigmaVector.h"

namespace Serenity {

template<Options::SCF_MODES T> class KernelSigmaVector: public SigmaVector<T> {
public:
  KernelSigmaVector(
      std::shared_ptr<SystemController> system,
      Eigen::MatrixXd& guessVector,
      std::shared_ptr<Kernel<T> > kernel);
  virtual ~KernelSigmaVector() = default;

protected:
  Eigen::MatrixXd calculateBlock(
      Eigen::MatrixXd& guess,
      std::vector<SpinPolarizedData<T,Eigen::MatrixXd> >& dens);

private:
  std::shared_ptr<Kernel<T> > _kernel;
  std::shared_ptr<BasisFunctionOnGridController> _basisFunctionOnGridController;

  std::shared_ptr<d2F_dRho2<T> > _d2FdRho2;
  std::shared_ptr<dF_dSigma<T> > _dFdSigma;
  std::shared_ptr<d2F_dSigma2<T> > _d2FdSigma2;
  std::shared_ptr<d2F_dRhodSigma<T> > _d2FdRhodSigma;
  std::shared_ptr<Gradient<DensityOnGrid<T> > > _densityGradient;

  std::shared_ptr<dF_dSigma<T> > _totalDFdSigma;
  std::shared_ptr<d2F_dSigma2<T> > _totalD2FdSigma2;
  std::shared_ptr<d2F_dRhodSigma<T> > _totalD2FdRhodSigma;
  std::shared_ptr<Gradient<DensityOnGrid<T> > > _totalDensityGradient;

  bool _gga;
  bool _nAddGGA;

  void contractGridBlock(
      unsigned int iBlock,
      const std::shared_ptr<BasisFunctionOnGridController::BasisFunctionBlockOnGridData> blockData,
      std::vector<SpinPolarizedData<T,Eigen::MatrixXd> >& dens,
      std::vector<SpinPolarizedData<T,Eigen::VectorXd> >& scalarPart,
      std::vector<Gradient<SpinPolarizedData<T,Eigen::VectorXd> > >* gradientPart = nullptr);
};

} /* namespace Serenity */

#endif /* POSTHF_RESPONSE_KERNELSIGMAVECTOR_H_ */
