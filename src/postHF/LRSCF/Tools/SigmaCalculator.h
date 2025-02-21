/**
 * @file   SigmaCalculator.h
 * @author Niklas Niemeyer
 * @date   Oct 13, 2021
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
#ifndef SIGMACALCULATOR_H
#define SIGMACALCULATOR_H
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <algorithm>
#include <memory>
#include <vector>

namespace Serenity {
/**
 * @class SigmaCalculator
 * @brief Used to perform linear matrix--vector multiplications for TDDFT.
 */
using SigmaCalculator =
    std::function<std::unique_ptr<std::vector<Eigen::MatrixXd>>(std::vector<Eigen::MatrixXd>& guessVectors)>;

/**
 * @class NonlinearSigmaCalculator
 * @brief Used to perform nonlinear matrix--vector multiplications for ADC(2)/CC2.
 */
using NonlinearSigmaCalculator =
    std::function<std::unique_ptr<Eigen::MatrixXd>(Eigen::MatrixXd& guessVectors, Eigen::VectorXd eigenvalues)>;

} /* namespace Serenity */
#endif /* SIGMACALCULATOR_H */
