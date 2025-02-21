/**
 * @file   CDIntegralController.cpp
 *
 * @date   Jun 28, 2018
 * @author Lars Hellmann
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
#include "integrals/CDIntegralController.h"
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
#include "basis/CombinedShellPair.h"
#include "basis/CustomBasisController.h"
#include "basis/Shell.h"
#include "dft/Functional.h"
#include "geometry/Geometry.h"
#include "integrals/CDStorageController.h"
#include "integrals/Normalization.h"
#include "integrals/RI_J_IntegralController.h"
#include "integrals/decomposer/CholeskyDecomposer.h"
#include "integrals/decomposer/TwoElecFourCenterIntDecomposer.h"
#include "integrals/decomposer/TwoElecTwoCenterIntDecomposer.h"
#include "integrals/looper/TwoElecFourCenterIntLooper.h"
#include "integrals/wrappers/Libint.h"
#include "io/FormattedOutputStream.h"
#include "math/IntegerMaths.h"
#include "math/RegularRankFourTensor.h"
#include "misc/WarningTracker.h"
/* Include Std and External Headers */
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Eigen/Eigen>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace Serenity {

CDIntegralController::CDIntegralController(const Settings& settings)
  : _settings(settings), _decompositionThreshold(_settings.basis.cdThreshold), _diskMode(std::make_shared<bool>(false)) {
}

bool CDIntegralController::getACDVectors(std::shared_ptr<BasisController> basisController,
                                         std::shared_ptr<BasisController> auxBasisController) {
  //  return false;
  auto aoController = this->getStorageController("ACDAO");

  auto n = basisController->getNBasisFunctions();
  auto m = auxBasisController->getNBasisFunctions();

  double memDemand = n * n;
  memDemand *= m;
  memDemand *= 8;
  memDemand *= 1e-9;

  std::cout << "\nApproximated memory demand to store the full AC(C)D Vectors (GB): " << memDemand << std::endl;

  auto memManager = MemoryManager::getInstance();
  double availableMem = memManager->getAvailableSystemMemory();
  availableMem *= 1e-9;

  if ((availableMem - memDemand) < 2) {
    std::cout << "\nStarting integral direct Cholesky routines." << std::endl;
    return false;
  }

  if (!aoController->getUpToDate())
    generateACDVectors(basisController, auxBasisController);
  if (aoController->getNVectors())
    return true;
  return false;
}

void CDIntegralController::generateACDVectors(std::shared_ptr<BasisController> basisController,
                                              std::shared_ptr<BasisController> auxBasisController) {
  bool normAux = !(auxBasisController->isAtomicCholesky());
  // Get the Storage controller to handle the generated vectors
  auto aoController = this->getStorageController("ACDAO");

  aoController->addSensitiveBasis(basisController);
  aoController->addSensitiveBasis(auxBasisController);

  OutputControl::nOut << "\tGenerating atomic(-compact) Cholesky Vectors\n" << std::endl;
  Timings::takeTime("Chol. -   generate aCD Vectors");

  double prescreeningThreshold = _settings.basis.integralThreshold;

  auto n = basisController->getNBasisFunctions();
  auto m = auxBasisController->getNBasisFunctions();
  auto nShells = basisController->getReducedNBasisFunctions();
  auto nAuxShells = auxBasisController->getReducedNBasisFunctions();

  // Calculate the diagonal of the ERI-matrix
  Eigen::VectorXd diagonal(n * n);
  diagonal.setZero();
  TwoElecFourCenterIntLooper diagLooper(LIBINT_OPERATOR::coulomb, 0, basisController, prescreeningThreshold);

  auto const calcDiagonal = [&](const unsigned int& i, const unsigned int& j, const Eigen::VectorXd& integral,
                                const unsigned int threadId) {
    (void)threadId; // no warnings, please
    unsigned int ij = i * n + j;
    unsigned int ji = j * n + i;
    diagonal(ij) = integral(0);
    diagonal(ji) = integral(0);
  };

  diagLooper.loopDiagonal(calcDiagonal);
  aoController->storeDiag(std::make_shared<Eigen::VectorXd>(diagonal));

  // calculate three center integrals and evaluate aCD Vectors for every element
  // Utilize existing RI framework to evaluate three index integrals
  auto& basis = basisController->getBasis();
  auto& auxbasis = auxBasisController->getBasis();
  const auto& riPrescreeningFactors = auxBasisController->getRIPrescreeningFactors();

  auto riints = std::make_shared<RI_J_IntegralController>(basisController, auxBasisController);

  Eigen::MatrixXd minv = riints->getInverseMSqrt();

  riints = nullptr;

  auto& libint = Libint::getInstance();
  libint.initialize(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.initialize(LIBINT_OPERATOR::coulomb, 0, 4);

  aoController->allocateVectors(m, n * n);

  for (unsigned int I = 0; I < nShells; I++) {
    const auto& basI = *basis[I];
    const unsigned int nI = basis[I]->getNContracted();
    const auto firstI = basisController->extendedIndex(I);
    Eigen::MatrixXd threeCInts(nI * n, m);
    threeCInts.setZero();

#ifdef _OPENMP
    std::vector<Eigen::MatrixXd> integrals(omp_get_max_threads());
#else
    std::vector<Eigen::MatrixXd> integrals(1);
#endif
#pragma omp parallel for schedule(dynamic)
    for (unsigned int J = 0; J < nShells; J++) {
#ifdef _OPENMP
      const unsigned int threadId = omp_get_thread_num();
#else
      const unsigned int threadId = 0;
#endif

      const auto& basJ = *basis[J];
      const unsigned int nJ = basis[J]->getNContracted();
      const auto firstJ = basisController->extendedIndex(J);

      bool significant = libint.compute(LIBINT_OPERATOR::coulomb, 0, basI, basJ, basI, basJ, integrals[threadId]);
      if (!significant)
        continue;

      double factor = sqrt(integrals[threadId].maxCoeff());
      for (unsigned int K = 0; K < nAuxShells; K++) {
        auto& k = (*riPrescreeningFactors)[K];

        assert(k.bf1 == K);

        if (factor * k.factor < prescreeningThreshold)
          continue;

        const auto& auxbasK = *auxbasis[K];
        const unsigned int nK = auxbasis[K]->getNContracted();
        const auto firstK = auxBasisController->extendedIndex(K);

        bool significant = libint.compute(LIBINT_OPERATOR::coulomb, 0, auxbasK, basI, basJ, integrals[threadId], normAux);

        if (significant) {
          unsigned int koffset = n * nI - nK;
          unsigned int ioffset = n - nJ;
          auto dataptr = threeCInts.data() + firstK * n * nI + firstJ;
          auto intptr = integrals[threadId].data();
          for (unsigned int k = 0; k < nK; k++) {
            for (unsigned int i = 0; i < nI; i++) {
              for (unsigned int j = 0; j < nJ; j++, intptr++) {
                (*dataptr) = (*intptr);
                dataptr++;
              }
              dataptr += ioffset;
            }
          }
          dataptr += koffset;
        }
      }
    }

    Eigen::MatrixXd acdCholVecBlocks = threeCInts * minv;

    // Store segment to disk to ensure efficient reordering for larger systems
    for (unsigned int l = 0; l < m; l++) {
      std::shared_ptr<std::vector<double>> tmpVecPtr = std::make_shared<std::vector<double>>(
          acdCholVecBlocks.data() + l * n * nI, acdCholVecBlocks.data() + (l + 1) * n * nI);

      aoController->storeSegment(l, firstI * n, tmpVecPtr);
    }
  }

  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 3);
  libint.finalize(LIBINT_OPERATOR::coulomb, 0, 4);

  aoController->setUpToDate();
  Timings::timeTaken("Chol. -   generate aCD Vectors");
}

void CDIntegralController::generateACDBasis(std::shared_ptr<Geometry> geom, std::string op_label, LIBINT_OPERATOR op) {
  for (unsigned int i = 0; i < geom->getNAtoms(); i++) {
    if ((*geom)[i]->getPrimaryBasisLabel() != (*geom)[0]->getPrimaryBasisLabel()) {
      throw SerenityError("Primary basis labels differ on atoms");
    }
  }
  this->generateACDBasis(geom, (*geom)[0]->getPrimaryBasisLabel(), op_label, op);
}

void CDIntegralController::generateACDBasis(std::shared_ptr<Geometry> geom, std::string label, std::string op_label,
                                            LIBINT_OPERATOR op) {
  OutputControl::nOut << "\tGenerating atomic Cholesky Basis" << std::endl;

  if (checkBasisFile(_settings.path + "/ACD-" + label + op_label, getAtomTypes(geom, label)))
    return;

  auto functional = _settings.customFunc.basicFunctionals.size() ? Functional(_settings.customFunc)
                                                                 : resolveFunctional(_settings.dft.functional);
  const double mu = functional.getRangeSeparationParameter();

  auto& libint = Libint::getInstance();
  libint.initialize(op, 0, 4);
  libint.initialize(op, 0, 2);

  Timings::takeTime("Chol. -     generate aCD Basis");

  auto atomTypes = getAtomTypes(geom, label);

  // Generate Basis file
  std::ofstream file;
  file.open(_settings.path + "/ACD-" + label + op_label);
  // store cout flags to restore later

  std::ios_base::fmtflags f(std::cout.flags());
  // Write file header and relevant settings to file
  file << "#Automatically generated atomic Cholesky Basis\n";
  file << "#Basis generated for Cholesky decomposition threshold = " << _settings.basis.cdThreshold << "\n";
  file << "#Basis generated for range-separation factor = " << mu << "\n";
  file << "\n\n\n$basis\n";
  // Loop over all found combinations
  for (auto const& atom : atomTypes) {
    // Perform aCD for found combination
    std::string atomLabel = atom.first;
    atomLabel += op_label;

    std::vector<std::shared_ptr<Shell>> shells = (*geom)[atom.second]->getBasisFunctions();

    // Generate Custom Basis Controller
    std::vector<std::shared_ptr<const Shell>> constShells;
    for (auto shell : shells) {
      const Shell tmpShell = *(shell);
      constShells.push_back(std::make_shared<const Shell>(tmpShell));
    }

    auto customBasisController = std::make_shared<CustomBasisController>(constShells, atomLabel);
    // Generate the Cholesky Basis
    std::shared_ptr<TwoElecFourCenterIntDecomposer> eriDecomposer = std::make_shared<TwoElecFourCenterIntDecomposer>(
        _settings, customBasisController, this->shared_from_this(), atomLabel);
    if (op == LIBINT_OPERATOR::erf_coulomb) {
      eriDecomposer = std::make_shared<TwoElecFourCenterIntDecomposer>(_settings, customBasisController,
                                                                       this->shared_from_this(), atomLabel, op, mu);
    }
    std::vector<unsigned int> cBasis = eriDecomposer->getCholeskyBasis();
    // sort Cholesky Basis
    std::sort(cBasis.begin(), cBasis.end());

    unsigned int nShells = (*geom)[atom.second]->getNBasisFunctions();

    // define Eigen Vector<unsigned int>
    typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> VectorUi;
    VectorUi firstShell(nShells + 1);

    unsigned int nbfs = customBasisController->getNBasisFunctions();
    firstShell.setZero();
    for (unsigned int I = 1; I < nShells + 1; ++I) {
      firstShell[I] = firstShell[I - 1] + shells[I - 1]->getNContracted();
    }

    // convert the combined indices (ij) that appear in the Cholesky basis into unique pairs of single indices (i and j)
    std::vector<std::pair<unsigned int, unsigned int>> selectedPairs;
    for (auto c : cBasis) {
      unsigned int i = std::floor((c) / nbfs);
      unsigned int j = c - (i * nbfs);
      if (i * nbfs + j != c)
        throw SerenityError("i and j not calculated correctly!");
      unsigned int counter = 0;
      while (firstShell[counter + 1] <= i)
        counter++;
      i = counter;
      counter = 0;
      while (firstShell[counter + 1] <= j)
        counter++;
      j = counter;

      auto tmpPair = std::make_pair(i, j);
      auto it = std::find_if(selectedPairs.begin(), selectedPairs.end(),
                             [&tmpPair](const std::pair<unsigned int, unsigned int>& element) {
                               return ((element.first == tmpPair.first && element.second == tmpPair.second) ||
                                       (element.first == tmpPair.second && element.second == tmpPair.first));
                             });
      if (it == selectedPairs.end()) {
        selectedPairs.push_back(tmpPair);
      }
    }

    // This vector keeps track of the shell splitting.
    // In the following blocks for every shell with the index k the following is done:
    //          1. If that shell contains more than 20 primitive functions it is split in blocks of max 20 primitive
    //          functions
    //               Each block forms a new shell and for each block the index k is pushed to the back of this vector.
    //          2. If this is for a spherical basis shells for lower angular momentum are added for a complete shell
    //          structure.
    //               For each shell generated this way the index k is pushed to the back of this vector.
    // This information allows for a recombination of these shells (i.e for the ACCD algorithm)
    std::vector<unsigned int> shellSplit;

    if (_settings.basis.makeSphericalBasis) {
      //==================================== //
      //         BLOCK FOR SPHERICAL BASIS SETS             //
      //====================================//
      std::vector<std::shared_ptr<const Shell>> expandedCombinedShells;
      for (unsigned int k = 0; k < selectedPairs.size(); k++) {
        std::shared_ptr<CombinedShellPair> tmpShellPair = std::make_shared<CombinedShellPair>(
            shells[selectedPairs[k].first], shells[selectedPairs[k].second], _settings.basis.makeSphericalBasis);

        std::vector<std::shared_ptr<CombinedShellPair>> completeSpherical =
            generateCompleteCombinedSphericalShell(tmpShellPair);
        for (unsigned int l = 0; l < completeSpherical.size(); l++) {
          std::vector<std::shared_ptr<const CombinedShellPair>> completeSplitSpherical =
              splitPrimitives(completeSpherical[l]);
          for (unsigned int m = 0; m < completeSplitSpherical.size(); m++) {
            expandedCombinedShells.push_back(completeSplitSpherical[m]);
            shellSplit.push_back(k);
          }
        }
      }

      std::shared_ptr<CustomBasisController> customExpandedBasisController(
          new CustomBasisController(expandedCombinedShells, atomLabel + "_EXPANDED"));

      // Generate the Cholesky Basis
      std::shared_ptr<TwoElecTwoCenterIntDecomposer> eriDecomposer1 = std::make_shared<TwoElecTwoCenterIntDecomposer>(
          _settings, customExpandedBasisController, this->shared_from_this(), atomLabel + "_EXPANDED");
      if (op == LIBINT_OPERATOR::erf_coulomb) {
        eriDecomposer1 = std::make_shared<TwoElecTwoCenterIntDecomposer>(
            _settings, customExpandedBasisController, this->shared_from_this(), atomLabel + "_EXPANDED", op, mu);
      }
      eriDecomposer1->setThreshold(_settings.basis.secondCD);
      std::vector<unsigned int> cBasisExpanded = eriDecomposer1->getCholeskyBasis();

      // sort Cholesky Basis
      std::sort(cBasisExpanded.begin(), cBasisExpanded.end());

      auto extendedBasis = customExpandedBasisController->getBasis();
      unsigned int nShells = extendedBasis.size();

      // define Eigen Vector<unsigned int>
      typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> VectorUi;
      VectorUi firstShell(nShells + 1);

      firstShell.setZero();
      for (unsigned int I = 1; I < nShells + 1; ++I) {
        firstShell[I] = firstShell[I - 1] + extendedBasis[I - 1]->getNContracted();
      }

      // convert the combined indices (ij) that appear in the Cholesky basis into unique pairs of single indices (i and j)
      std::vector<unsigned int> selectedShells;
      for (auto c : cBasisExpanded) {
        unsigned int counter = 0;
        while (firstShell[counter + 1] <= c)
          counter++;
        auto it = std::find(selectedShells.begin(), selectedShells.end(), counter);
        if (it == selectedShells.end()) {
          selectedShells.push_back(counter);
        }
      }

      // select elements present in selectedShells from the shellSplit
      std::vector<unsigned int> selectedSplit;
      for (auto i : selectedShells)
        selectedSplit.push_back(shellSplit[i]);

      // Construct the aCD Basis for each combination
      // First write the shell splitting
      file << "#The following shells might have been split during construction and can be \n"
           << "#recombined if wanted (only for matching angular momentum)(spherical):\n"
           << "#Split: ";
      for (auto i : selectedSplit)
        file << i << " ";
      file << "\n";
      // Now write the actual basis
      file << "*\n"
           << atom.first << "   "
           << "ACD-" << label << op_label << "\n*\n";
      // Loop over selectedPairs
      for (unsigned int k = 0; k < selectedShells.size(); k++) {
        // Limitation due to the limits of libint. However, contributions are so small they can be neglected
        if (extendedBasis[selectedShells[k]]->getAngularMomentum() > AM_MAX) {
          continue;
        }

        auto expo = extendedBasis[selectedShells[k]]->getExponents();
        auto cont = extendedBasis[selectedShells[k]]->getContractions();
        auto nPrim = extendedBasis[selectedShells[k]]->getNPrimitives();
        char angularMomChar = getAngularMomentumChar(extendedBasis[selectedShells[k]]->getAngularMomentum());

        file << "\t" << int(nPrim) << "  " << angularMomChar << "\n";
        for (unsigned int j = 0; j < nPrim; j++) {
          file << std::fixed << std::setprecision(12) << "\t  " << expo[j] << "\t\t\t" << std::fixed
               << std::setprecision(14) << cont[j] << "\n";
        }

      } // end of loop over selected pairs
    }
    else if (!_settings.basis.makeSphericalBasis) {
      //==================================== //
      //         BLOCK FOR CARTESIAN BASIS SETS             //
      //====================================//

      // Hardcoded this to 20 because changes to this number in Libint.h caused the LRSCFTaskCC2Test.CD_ADC2 to fail,
      // which was traced back to different shell splitting behaviour of the aCD basis (Niklas Göllmann, Oct. 2024)

      // unsigned int maxNPrimCart = Libint::getNPrimMax() - (Libint::getNPrimMax() % 3);
      unsigned int maxNPrimCart = 20 - (20 % 3);

      // Write shell split
      for (unsigned int k = 0; k < selectedPairs.size(); k++) {
        // create product-shell of the two shells in the selected shellpair

        unsigned int am = shells[selectedPairs[k].first]->getAngularMomentum();
        am += shells[selectedPairs[k].second]->getAngularMomentum();
        while (am > AM_MAX)
          am -= 2;

        CombinedShellPair tmpShellPair(shells[selectedPairs[k].first], shells[selectedPairs[k].second], am,
                                       _settings.basis.makeSphericalBasis);
        // Limitation due to the limits of libint. However, contributions are so small they can be neglected
        if (tmpShellPair.getAngularMomentum() > AM_MAX) {
          continue;
        }
        int nPrim = tmpShellPair.getNPrimitives();
        while (nPrim > 0) {
          shellSplit.push_back(k);
          nPrim -= maxNPrimCart;
        }
      }
      // First write the shell splitting
      file << "#The following shells might have been split during construction and can be \n"
           << "#recombined if wanted (only for matching angular momentum) (cartesian):\n"
           << "#Split: ";
      for (auto i : shellSplit)
        file << i << " ";
      file << "\n";

      // Construct the aCD Basis for each combination
      file << "*\n"
           << atom.first << "   "
           << "ACD-" << label << op_label << "\n*\n";
      // Loop over selectedPairs
      for (unsigned int k = 0; k < selectedPairs.size(); k++) {
        // creat product-shell of the two shells in the selected shallpair
        unsigned int am = shells[selectedPairs[k].first]->getAngularMomentum();
        am += shells[selectedPairs[k].second]->getAngularMomentum();
        while (am > AM_MAX)
          am -= 2;

        CombinedShellPair tmpShellPair(shells[selectedPairs[k].first], shells[selectedPairs[k].second], am,
                                       _settings.basis.makeSphericalBasis);
        // Limitation due to the limits of libint. However, contributions are so small they can be neglected
        if (tmpShellPair.getAngularMomentum() > AM_MAX) {
          continue;
        }
        unsigned int nPrim = tmpShellPair.getNPrimitives();
        if (nPrim > maxNPrimCart) { // split the shell if nPrim is larger than libint->N_PRIM_MAX
          libint2::svector<double> contr = tmpShellPair.getNormContractions();
          libint2::svector<double> expo = tmpShellPair.getExponents();
          unsigned int i = 0;
          unsigned int step = maxNPrimCart;
          while (i < nPrim) {
            unsigned int am = tmpShellPair.getAngularMomentum();

            if (i + step >= nPrim)
              step = nPrim - i;
            std::vector<double> split_contr(contr.begin() + i, contr.begin() + i + step);
            std::vector<double> split_expo(expo.begin() + i, expo.begin() + i + step);
            std::shared_ptr<const CombinedShellPair> splitShellPair =
                std::make_shared<const CombinedShellPair>(shells[selectedPairs[k].first], shells[selectedPairs[k].second],
                                                          split_expo, split_contr, am, _settings.basis.makeSphericalBasis);
            auto expo = splitShellPair->getExponents();
            auto cont = splitShellPair->getContractions();
            auto nPrim = splitShellPair->getNPrimitives();
            char angularMomChar = getAngularMomentumChar(splitShellPair->getAngularMomentum());

            file << "\t" << int(nPrim) << "  " << angularMomChar << "\n";
            for (unsigned int j = 0; j < nPrim; j++) {
              file << std::fixed << std::setprecision(12) << "\t  " << expo[j] << "\t\t\t" << std::fixed
                   << std::setprecision(14) << cont[j] << "\n";
            }

            i += step;
          }
        }
        else {
          auto expo = tmpShellPair.getExponents();
          auto cont = tmpShellPair.getContractions();
          auto nPrim = tmpShellPair.getNPrimitives();
          char angularMomChar = getAngularMomentumChar(tmpShellPair.getAngularMomentum());

          file << "\t" << int(nPrim) << "  " << angularMomChar << "\n";
          for (unsigned int j = 0; j < nPrim; j++) {
            file << std::fixed << std::setprecision(12) << "\t  " << expo[j] << "\t\t\t" << std::fixed
                 << std::setprecision(14) << cont[j] << "\n";
          }
        }

      } // end of loop over selected pairs
    }
    else {
      throw SerenityError("custom ACD basis is not pure spherical or cartesian");
    }
  } // end loop over atom types
  file << "*\n$end";
  // restore state of cout
  std::cout.flags(f);
  // close file
  file.close();

  libint.finalize(op, 0, 4);
  libint.finalize(op, 0, 2);

  Timings::timeTaken("Chol. -     generate aCD Basis");
}

void CDIntegralController::generateACCDBasis(std::shared_ptr<Geometry> geom, std::string op_label, LIBINT_OPERATOR op) {
  for (unsigned int i = 0; i < geom->getNAtoms(); i++) {
    if ((*geom)[i]->getPrimaryBasisLabel() != (*geom)[0]->getPrimaryBasisLabel()) {
      throw SerenityError("Primary basis labels differ on atoms");
    }
  }
  this->generateACCDBasis(geom, (*geom)[0]->getPrimaryBasisLabel(), op_label, op);
}

void CDIntegralController::generateACCDBasis(std::shared_ptr<Geometry> geom, std::string label, std::string op_label,
                                             LIBINT_OPERATOR op) {
  OutputControl::nOut << "\tGenerating atomic-compact Cholesky Basis" << std::endl;
  if (_settings.basis.cdThreshold > 1e-4)
    WarningTracker::printWarning(
        "Threshold chosen for the accd procedure is larger than 1E-4. This might lead to unreliable results.", true);

  if (checkBasisFile(_settings.path + "/ACCD-" + label + op_label, getAtomTypes(geom, label)))
    return;

  Timings::takeTime("Chol. -    generate acCD Basis");
  auto atomTypes = getAtomTypes(geom, label);

  auto functional = _settings.customFunc.basicFunctionals.size() ? Functional(_settings.customFunc)
                                                                 : resolveFunctional(_settings.dft.functional);
  const double mu = functional.getRangeSeparationParameter();
  auto& libint = Libint::getInstance();
  if (op == LIBINT_OPERATOR::erf_coulomb) {
    libint.initialize(op, 0, 4, std::vector<std::shared_ptr<Atom>>(0), mu);
    libint.initialize(op, 0, 2, std::vector<std::shared_ptr<Atom>>(0), mu);
  }
  else {
    libint.initialize(op, 0, 4);
    libint.initialize(op, 0, 2);
  }
  Shell::do_enforce_unit_normalization(false);
  //  Shell::do_enforce_unit_normalization(true);

  // 2. Generate Basis file
  // Create basis file
  std::ofstream file;
  file.open(_settings.path + "/ACCD-" + label + op_label);
  // store cout flags to restore later
  std::ios_base::fmtflags f(std::cout.flags());
  // Write file header and relevant settings to file
  file << "#Automatically generated atomic Cholesky Basis\n";
  file << "#Basis generated for Cholesky decomposition threshold = " << _settings.basis.cdThreshold << "\n";
  file << "#Basis generated for range-separation factor = " << mu << "\n";
  file << "\n\n\n$basis\n";

  // Loop over all atoms combinations
  for (auto atom : atomTypes) {
    file << "*\n"
         << atom.first << "   "
         << "ACCD-" << label << op_label << "\n*\n";
    // loop over all basis types for the atom
    // setup basis to collect all relevant primitive shells for one atom
    std::vector<std::shared_ptr<const Shell>> newAtomBas;

    std::vector<std::shared_ptr<Shell>> shells = (*geom)[atom.second]->getBasisFunctions("ACD-" + label + op_label);

    {
      // recombine the shells as defined in the ACD basis set file

      // get the shell split for the current atom
      std::vector<unsigned int> shellSplit = getShellSplit(_settings.path + "/ACD-" + label + op_label, atom.first);

      // recombine the shells with the same splitting index and angular momentum
      std::vector<std::shared_ptr<Shell>> tmpShells;
      libint2::svector<double> contr = shells[0]->getContractions();
      libint2::svector<double> expo = shells[0]->getExponents();
      unsigned int am = shells[0]->getAngularMomentum();
      unsigned int shellNo = shellSplit[0];
      unsigned int i = 1;
      while (i < shellSplit.size()) {
        while (i < shellSplit.size() and am == shells[i]->getAngularMomentum() and shellNo == shellSplit[i]) {
          libint2::svector<double> tmpcontr = shells[i]->getContractions();
          libint2::svector<double> tmpexpo = shells[i]->getExponents();
          for (unsigned int j = 0; j < tmpexpo.size(); j++) {
            contr.push_back(tmpcontr[j]);
            expo.push_back(tmpexpo[j]);
          }
          i++;
        }
        std::array<double, 3> coords = {shells[i - 1]->getX(), shells[i - 1]->getY(), shells[i - 1]->getZ()};
        std::shared_ptr<Shell> newShell =
            std::make_shared<Shell>(expo, contr, am, shells[i - 1]->isSpherical(), coords, shells[i - 1]->getElement());
        tmpShells.push_back(newShell);
        if (i < shellSplit.size()) {
          am = shells[i]->getAngularMomentum();
          shellNo = shellSplit[i];
        }
        contr.clear();
        expo.clear();
        //      i++;
      }
      // replace the split shells with the recombined ones
      shells = tmpShells;
    }

    for (unsigned int i = 0; i < shells.size(); i++) {
      /*
       * For each shell setup a custom temporary basiscontroller that contains
       * each basis function from the original shell in a new shell on its own.
       */
      std::vector<std::shared_ptr<const Shell>> newBas;

      const auto& shell = shells[i];

      auto contr = shell->getContractions();
      auto expo = shell->getExponents();

      libint2::svector<double> newExpo;
      newExpo.push_back(0.0);
      libint2::svector<double> newContr;
      newContr.push_back(0.0);

      unsigned int angularMomentum = shell->getAngularMomentum();
      bool isSpherical = shell->isSpherical();
      std::array<double, 3> coords;
      coords[0] = shell->getX();
      coords[1] = shell->getY();
      coords[2] = shell->getZ();
      std::string element = shell->getElement();

      for (unsigned int j = 0; j < contr.size(); j++) {
        newExpo[0] = expo[j];
        newContr[0] = contr[j];
        newBas.push_back(std::make_shared<const Shell>(newExpo, newContr, angularMomentum, isSpherical, coords, element));
      }

      std::shared_ptr<CustomBasisController> basController(new CustomBasisController(newBas, "customShellWiseBasis"));

      auto nShells = basController->getReducedNBasisFunctions();

      auto nbfs = basController->getNBasisFunctions();

      const auto& basis = basController->getBasis();

      /*
       * calculate the two-center coulomb integrals. (Remember the basis is the ACD-basis, meaning that every index
       * corresponds to two original indices. Thus, the two-center integral calculated here actually correspond to
       * four-center integrals.)
       */

      if (op == LIBINT_OPERATOR::erf_coulomb) {
        libint.initialize(op, 0, 2, std::vector<std::shared_ptr<Atom>>(0), mu);
      }
      else {
        libint.initialize(op, 0, 2);
      }

      Eigen::MatrixXd twoC(nShells, nShells);
      twoC.setZero();

      Eigen::MatrixXd twoCFull(nbfs, nbfs);
      twoCFull.setZero();

      std::vector<Eigen::MatrixXd> ints(1);

      for (unsigned int i = 0; i < nShells; i++) {
        const auto& basI = *basis[i];
        const unsigned int nI = basis[i]->getNContracted();
        const auto& firstI = basController->extendedIndex(i);

        for (unsigned int j = 0; j < nShells; j++) {
          const auto& basJ = *basis[j];
          const unsigned int nJ = basis[j]->getNContracted();
          const auto& firstJ = basController->extendedIndex(j);

          bool significant = libint.compute(op, 0, basI, basJ, ints[0], false);
          if (!significant)
            continue;

          for (unsigned int ii = 0; ii < nI; ii++) {
            for (unsigned int jj = 0; jj < nJ; jj++) {
              twoC(i, j) += std::abs(ints[0].col(0)[jj + nJ * ii]);
              twoCFull(firstI + ii, firstJ + jj) = ints[0].col(0)[jj + nJ * ii];
            }
          }
        }
      }

      libint.finalize(op, 0, 2);

      std::vector<unsigned int> cbasShells;

      {
        auto matCopy = twoCFull;
        auto vec = twoCFull;
        vec.setZero();
        Eigen::VectorXd diag = matCopy.diagonal();

        double thresh = _settings.basis.cdThreshold;

        unsigned int index = 0;
        unsigned int counter = 0;
        std::vector<unsigned int> cbas;
        while (true) {
          if (diag.maxCoeff(&index) <= thresh)
            break;

          cbas.push_back(index);
          double factor = 1.0 / std::sqrt(diag[index]);
          vec.col(counter) = factor * matCopy.col(index);
          matCopy -= vec.col(counter) * vec.col(counter).transpose();
          matCopy.col(index).setZero();
          matCopy.row(index).setZero();

          diag = matCopy.diagonal();

          counter++;
        }

        for (auto i : cbas) {
          cbasShells.push_back(basController->reducedIndex(i));
        }
        std::sort(cbasShells.begin(), cbasShells.end());
        cbasShells.erase(std::unique(cbasShells.begin(), cbasShells.end()), cbasShells.end());
      }

      // Decompose reduced integral matrix for comparison
      // This decomposition is done in a direct implementation, since the matrix twoC will always fit in memory
      auto matCopy = twoC;
      auto vec = twoC;
      vec.setZero();
      Eigen::VectorXd diag = matCopy.diagonal();

      double thresh = _settings.basis.cdThreshold;

      unsigned int index = 0;
      unsigned int counter = 0;
      std::vector<unsigned int> cbas;
      while (true) {
        if (diag.maxCoeff(&index) <= thresh)
          break;

        cbas.push_back(index);
        double factor = 1.0 / std::sqrt(diag[index]);
        vec.col(counter) = factor * matCopy.col(index);
        matCopy -= vec.col(counter) * vec.col(counter).transpose();
        matCopy.col(index).setZero();
        matCopy.row(index).setZero();

        diag = matCopy.diagonal();

        counter++;
      }

      cbas = cbasShells;
      // sort Cholesky Basis
      std::sort(cbas.begin(), cbas.end());

      if (cbas.size() == 0) {
        continue;
      }

      // create custom basis controller of accd basis to determine fitted contraction coefficients
      std::vector<std::shared_ptr<const Shell>> newACCDBas;
      for (unsigned int k = 0; k < cbas.size(); k++) {
        newACCDBas.push_back(basis[cbas[k]]);
      }
      std::shared_ptr<CustomBasisController> accdNonFittedBasCont(new CustomBasisController(newACCDBas, "newACCDBas"));

      Eigen::MatrixXd mixedTwoC;
      {
        auto nShells = basController->getReducedNBasisFunctions();
        const auto& basis = basController->getBasis();

        auto nShellsACCD = accdNonFittedBasCont->getReducedNBasisFunctions();
        const auto& basisACCD = accdNonFittedBasCont->getBasis();

        if (op == LIBINT_OPERATOR::erf_coulomb) {
          libint.initialize(op, 0, 2, std::vector<std::shared_ptr<Atom>>(0), mu);
        }
        else {
          libint.initialize(op, 0, 2);
        }

        Eigen::MatrixXd twoC(nShellsACCD, nShells);
        twoC.setZero();

        std::vector<Eigen::MatrixXd> ints(1);

        for (unsigned int i = 0; i < nShellsACCD; i++) {
          const auto& basI = *basisACCD[i];
          const unsigned int nI = basisACCD[i]->getNContracted();

          for (unsigned int j = 0; j < nShells; j++) {
            const auto& basJ = *basis[j];
            const unsigned int nJ = basis[j]->getNContracted();

            bool significant = libint.compute(op, 0, basI, basJ, ints[0], false);
            if (!significant)
              continue;

            for (unsigned int ii = 0; ii < nI; ii++) {
              for (unsigned int jj = 0; jj < nJ; jj++) {
                twoC(i, j) += std::abs(ints[0].col(0)[jj + nJ * ii]);
              }
            }
          }
        }

        libint.finalize(op, 0, 2);

        mixedTwoC = twoC.transpose();
      }

      // Calculate the new scalled coefficients after primitive basis functions have been removed
      Eigen::MatrixXd coeffs = mixedTwoC.fullPivHouseholderQr().solve(twoC);
      Eigen::VectorXd redCoeffs = coeffs.rowwise().sum();

      libint2::svector<double> contractions;
      libint2::svector<double> exponents;

      for (unsigned int k = 0; k < cbas.size(); k++) {
        contractions.push_back(basis[cbas[k]]->getContractions()[0] * redCoeffs[k]);
        exponents.push_back(basis[cbas[k]]->getExponents()[0]);
      }

      auto angularMom = basis[0]->getAngularMomentum();
      auto angMomChar = getAngularMomentumChar(angularMom);
      file << "\t" << int(contractions.size()) << "  " << angMomChar << "\n";
      for (unsigned int j = 0; j < contractions.size(); j++) {
        file << std::fixed << std::setprecision(12) << "\t  " << exponents[j] << "\t\t\t" << std::fixed
             << std::setprecision(14) << contractions[j] << "\n";
      }
    } // loop over shells
  }
  file << "*\n$end";
  // restore state of cout
  std::cout.flags(f);
  // close output file
  file.close();
  libint.finalize(op, 0, 4);
  libint.finalize(op, 0, 2);
  Shell::do_enforce_unit_normalization(true);
  Timings::timeTaken("Chol. -    generate acCD Basis");
}

char CDIntegralController::getAngularMomentumChar(unsigned int angularMom) {
  char angularMomChar = 'z';
  switch (angularMom) {
    case 0:
      angularMomChar = 's';
      break;
    case 1:
      angularMomChar = 'p';
      break;
    case 2:
      angularMomChar = 'd';
      break;
    case 3:
      angularMomChar = 'f';
      break;
    case 4:
      angularMomChar = 'g';
      break;
    case 5:
      angularMomChar = 'h';
      break;
    case 6:
      angularMomChar = 'i';
      break;
    case 7:
      angularMomChar = 'k';
      break;
    case 8:
      angularMomChar = 'm';
      break;
    case 9:
      angularMomChar = 'n';
      break;
    case 10:
      angularMomChar = 'o';
      break;
  }
  if (angularMomChar == 'z')
    throw SerenityError("A angular momentum larger than the maximum angular momentum supported by the currently active "
                        "Libint-settings was detected during aCD-basis construction!");
  return angularMomChar;
}

std::vector<std::pair<std::string, unsigned int>> CDIntegralController::getAtomTypes(std::shared_ptr<Geometry> geom,
                                                                                     std::string label) {
  std::vector<std::pair<std::string, unsigned int>> atomTypes;
  for (unsigned int i = 0; i < geom->getNAtoms(); i++) {
    auto name = (*geom)[i]->getAtomType()->getName();

    // Check for ghostatoms and convert to normal atom
    if (name.substr(name.size() - 1, name.size()) == ":") {
      name = name.substr(0, name.size() - 1);
    }

    if ((*geom)[i]->basisFunctionsExist(label)) {
      std::transform(name.begin(), name.end(), name.begin(), ::tolower);
      // Check if atom is not already present
      auto tmpPair = std::make_pair(name, i);
      auto it = std::find_if(atomTypes.begin(), atomTypes.end(), [&tmpPair](const std::pair<std::string, unsigned int>& atom) {
        return ((atom.first == tmpPair.first));
      });
      if (it == atomTypes.end()) {
        atomTypes.push_back(tmpPair);
      }
    }
    else {
      WarningTracker::printWarning("\nWarning: No basis (" + label + ") for " + name + "!", true);
    }
  }
  return atomTypes;
}

bool CDIntegralController::checkBasisFile(std::string path, std::vector<std::pair<std::string, unsigned int>> atomTypes) {
  auto functional = _settings.customFunc.basicFunctionals.size() ? Functional(_settings.customFunc)
                                                                 : resolveFunctional(_settings.dft.functional);
  const double mu = functional.getRangeSeparationParameter();
  // check if there is already a correct basis in the systemfolder
  std::ifstream fileCheck(path);
  if (fileCheck.good()) {
    // check for correct settings
    try {
      std::string line;
      std::getline(fileCheck, line);
      std::getline(fileCheck, line);
      unsigned int c = line.size();
      unsigned int i;
      for (i = c; i > 0; i--) {
        if (line[i] == ' ')
          break;
      }
      double thresh = std::stod(line.substr(i + 1, c));

      std::getline(fileCheck, line);
      double muRead;
      {
        unsigned int c = line.size();
        unsigned int i;
        for (i = c; i > 0; i--) {
          if (line[i] == ' ')
            break;
        }
        muRead = std::stod(line.substr(i + 1, c));
      }

      if (thresh == _settings.basis.cdThreshold and muRead == mu) {
        // check for all atoms
        std::vector<std::string> atomsInBasis;
        while (std::getline(fileCheck, line)) {
          if (line[0] == '*') {
            std::getline(fileCheck, line);
            if (isalpha(line[0])) {
              std::string element = line.substr(0, 2);
              std::string::iterator end_pos = std::remove(element.begin(), element.end(), ' ');
              element.erase(end_pos, element.end());
              atomsInBasis.push_back(element);
            }
            std::getline(fileCheck, line);
          }
        }

        for (auto atom : atomTypes) {
          std::vector<std::string>::iterator it;
          it = find(atomsInBasis.begin(), atomsInBasis.end(), atom.first);
          if (it == atomsInBasis.end()) {
            return false;
          }
        }

        // if good return true
        return true;
      }
    }
    catch (...) {
      return false;
    }
  }
  return false;
}

std::vector<unsigned int> CDIntegralController::getShellSplit(std::string path, std::string element) {
  // The array to hold the read shell split
  std::vector<unsigned int> shellSplit;
  std::ifstream fileCheck(path);
  if (fileCheck.good()) {
    try {
      std::string line;
      std::getline(fileCheck, line);

      // Read all lines
      while (std::getline(fileCheck, line)) {
        // If the split is found read it into the return vector
        if (line.substr(0, 7) == "#Split:") {
          shellSplit.clear();

          std::stringstream linestream(line);
          std::istream_iterator<std::string> begin(linestream);
          std::istream_iterator<std::string> end;
          std::vector<std::string> vstrings(begin, end);
          for (unsigned int i = 1; i < vstrings.size(); i++) {
            shellSplit.push_back(std::stoul(vstrings[i]));
          }
        }

        // If the element requested is found return the current shellSplit
        //(This works since the split is always written before the actual basis in the file)
        if (line[0] == '*') {
          std::getline(fileCheck, line);
          if (isalpha(line[0])) {
            std::string readElement = line.substr(0, 2);
            std::string::iterator end_pos = std::remove(readElement.begin(), readElement.end(), ' ');
            readElement.erase(end_pos, readElement.end());
            if (element == readElement)
              return shellSplit;
          }
          std::getline(fileCheck, line);
        }
      }
    }
    catch (...) {
      throw SerenityError("Could not fetch shell split from ACD basis file");
    }
  }
  // If you get here something failed while reading. We return an empty string as the basis still can be used without
  // recombining.
  shellSplit.clear();
  return shellSplit;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> CDIntegralController::generatePseudoCoefficients(Eigen::MatrixXd p) {
  // Generate pseudo Coefficients from the DensityMatrix
  // Since for Densitydifferences this can lead to negative eigenvalues
  // we need to introduce a sign vector than can later be used to reconstruct
  // the correct density matrix from only the left-hand side.
  Eigen::EigenSolver<Eigen::MatrixXd> es(p);

  Eigen::VectorXd eigenVals = es.eigenvalues().real();

  Eigen::MatrixXd eigenVec = es.eigenvectors().real();
  Eigen::VectorXd esSigns(eigenVals.size());
  double numThreshold = 1e-14;
  for (unsigned int i = 0; i < eigenVals.size(); i++) {
    if (eigenVals[i] > numThreshold)
      esSigns[i] = 1;
    else if (eigenVals[i] < -1 * numThreshold) {
      eigenVals[i] *= -1;
      esSigns[i] = -1;
    }
    else {
      esSigns[i] = 0;
      eigenVals[i] = 0.0;
    }
  }

  { // Truncate the vectors and eigenvalues
    unsigned int count = (esSigns.array() != 0.0).count();
    Eigen::VectorXd eigenValsTrunc(count);
    Eigen::MatrixXd eigenVecTrunc(eigenVec.rows(), count);
    Eigen::VectorXd esSignsTrunc(count);
    eigenValsTrunc.setZero();
    eigenVecTrunc.setZero();

    unsigned int j = 0;
    for (unsigned int i = 0; i < eigenVals.size(); i++) {
      if (esSigns[i] != 0.0) {
        eigenValsTrunc[j] = eigenVals[i];
        eigenVecTrunc.col(j) = eigenVec.col(i);
        esSignsTrunc[j] = esSigns[i];
        j++;
      }
    }
    eigenVals = eigenValsTrunc;
    eigenVec = eigenVecTrunc;
    esSigns = esSignsTrunc;
  }

  Eigen::VectorXd eigenValsSqrt = eigenVals.cwiseSqrt();
  Eigen::MatrixXd occCoeff = (eigenVec * eigenValsSqrt.asDiagonal());
  return std::make_pair(esSigns, occCoeff);
}

std::vector<std::shared_ptr<CombinedShellPair>>
CDIntegralController::generateCompleteCombinedSphericalShell(std::shared_ptr<CombinedShellPair> tmpShellPair) {
  bool expandExponents = true;
  if (_settings.basis.extendSphericalACDShells == Options::EXTEND_ACD::SIMPLE)
    expandExponents = false;
  std::vector<std::shared_ptr<CombinedShellPair>> completeShell;

  auto shellI = tmpShellPair->getBaseShellA();
  auto shellJ = tmpShellPair->getBaseShellB();

  completeShell.push_back(tmpShellPair);

  if (_settings.basis.extendSphericalACDShells == Options::EXTEND_ACD::NONE)
    return completeShell;

  if (tmpShellPair->isCartesian())
    return completeShell;

  unsigned int am = tmpShellPair->getAngularMomentum();

  libint2::svector<double> lastLContr = tmpShellPair->getNormContractions();
  libint2::svector<double> lastLExpo = tmpShellPair->getExponents();
  std::vector<double> lastContr(lastLContr.begin(), lastLContr.end());
  std::vector<double> lastExpo(lastLExpo.begin(), lastLExpo.end());

  double offset = _settings.basis.cdOffset;
  double op = 1.0 + offset;
  double om = 1 / (1.0 + offset);

  int lA = tmpShellPair->getBaseShellA()->getAngularMomentum();
  int lB = tmpShellPair->getBaseShellB()->getAngularMomentum();
  unsigned int minL = std::abs(lA - lB);

  while (am > minL + 1) {
    am -= 2;

    std::vector<double> addExpo;
    std::vector<double> addContr;

    for (unsigned int i = 0; i < lastContr.size(); i++) {
      // generate additional contractions and exponents
      if (expandExponents) {
        addExpo.push_back(lastExpo[i]);
        addExpo.push_back(op * lastExpo[i]);
        addExpo.push_back(om * lastExpo[i]);

        double frac = (om * exp(-1 * om) - exp(-1)) / (exp(-1) - op * exp(-1 * op));
        double d = lastContr[i] / (lastExpo[i] * (1 - om + (1 - op) * frac));

        double c = frac * d;
        double b = -1 * (c + d);

        addContr.push_back(b);
        addContr.push_back(c);
        addContr.push_back(d);
      }
      else {
        addExpo.push_back(lastExpo[i]);
        addContr.push_back(lastContr[i]);
      }
    }
    // Check for identical exponents and summarize
    std::vector<std::pair<double, double>> pairVector;
    for (unsigned int i = 0; i < addExpo.size(); i++) {
      pairVector.push_back(std::make_pair(addExpo[i], addContr[i]));
    }
    std::sort(pairVector.begin(), pairVector.end());

    addExpo.clear();
    addContr.clear();
    double expo = pairVector[0].first;
    double contr = pairVector[0].second;
    for (unsigned int i = 1; i < pairVector.size(); i++) {
      if (pairVector[i].first == expo) {
        contr += pairVector[i].second;
      }
      else {
        addExpo.push_back(expo);
        addContr.push_back(contr);
        expo = pairVector[i].first;
        contr = pairVector[i].second;
      }
    }
    addExpo.push_back(expo);
    addContr.push_back(contr);

    // add new complete shell
    std::shared_ptr<CombinedShellPair> addShellPair =
        std::make_shared<CombinedShellPair>(shellI, shellJ, addExpo, addContr, am, _settings.basis.makeSphericalBasis);
    completeShell.push_back(addShellPair);

    if (_settings.basis.extendSphericalACDShells == Options::EXTEND_ACD::FIRST)
      expandExponents = false;

    lastContr = addContr;
    lastExpo = addExpo;
  }

  return completeShell;
}

std::vector<std::shared_ptr<const CombinedShellPair>>
CDIntegralController::splitPrimitives(std::shared_ptr<CombinedShellPair> tmpShellPair) {
  std::vector<std::shared_ptr<const CombinedShellPair>> splitShells;

  unsigned int nPrim = tmpShellPair->getNPrimitives();

  // Hardcoded this to 20 because changes to this number in Libint.h caused the LRSCFTaskCC2Test.CD_ADC2 to fail, which
  // was traced back to different shell splitting behaviour of the aCD basis (Niklas Göllmann, Oct. 2024)

  // unsigned int nPrimMax = Libint::getNPrimMax() - (Libint::getNPrimMax() % 3);
  unsigned int nPrimMax = 20 - (20 % 3);

  auto shellI = tmpShellPair->getBaseShellA();
  auto shellJ = tmpShellPair->getBaseShellB();

  if (nPrim > nPrimMax) { // split the shell if nPrim is larger than libint->N_PRIM_MAX
    libint2::svector<double> contr = tmpShellPair->getNormContractions();
    libint2::svector<double> expo = tmpShellPair->getExponents();
    unsigned int i = 0;
    unsigned int step = nPrimMax;

    while (i < nPrim) {
      unsigned int am = tmpShellPair->getAngularMomentum();
      if (i + step >= nPrim)
        step = nPrim - i;
      std::vector<double> split_contr(contr.begin() + i, contr.begin() + i + step);
      std::vector<double> split_expo(expo.begin() + i, expo.begin() + i + step);

      std::shared_ptr<const CombinedShellPair> splitShellPair = std::make_shared<const CombinedShellPair>(
          shellI, shellJ, split_expo, split_contr, am, _settings.basis.makeSphericalBasis);
      splitShells.push_back(splitShellPair);

      i += step;
    }
  }
  else {
    splitShells.push_back(tmpShellPair);
  }
  return splitShells;
}

/*------------------
 * getter and setter
 *-----------------*/

double CDIntegralController::getDecompositionThreshold() {
  return _decompositionThreshold;
}

std::shared_ptr<CDStorageController> CDIntegralController::getStorageController(std::string label) {
  if (_cdStorageController.find(label) == _cdStorageController.end())
    this->addMatrix(label);
  return _cdStorageController[label];
}

void CDIntegralController::cleanup() {
  for (std::map<std::string, std::shared_ptr<CDStorageController>>::iterator it = _cdStorageController.begin();
       it != _cdStorageController.end(); it++) {
    if (!it->second.unique())
      throw SerenityError("shared_ptr not unique while trying to delete completely!");
    it->second.reset();
  }
  _cdStorageController.clear();
  unsetDiskMode();
}

void CDIntegralController::addMatrix(std::string label) {
  _cdStorageController.emplace(
      label, std::make_shared<CDStorageController>(_settings.path + _settings.name, label, this->shared_from_this()));
}

std::shared_ptr<bool> CDIntegralController::getDiskMode() {
  return _diskMode;
}

void CDIntegralController::setDiskMode() {
  if (*_diskMode)
    return;
  for (std::map<std::string, std::shared_ptr<CDStorageController>>::iterator it = _cdStorageController.begin();
       it != _cdStorageController.end(); ++it) {
    it->second->setDiskMode();
  }
  (*_diskMode) = true;
}

void CDIntegralController::unsetDiskMode() {
  (*_diskMode) = false;
}

} /* namespace Serenity */
