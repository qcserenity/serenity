/**
 * @file   Functional.h
 * @author Thomas Dresselhaus <t.dresselhaus at wwu.de>
 *
 * @date   12. Juli 2015, 18:27
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
#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H
/* Include Std and External Headers */
#include <vector>

namespace Serenity {
/* Forward Declarations */
namespace BasicFunctionals {
enum class BASIC_FUNCTIONALS;
}
namespace CompositeFunctionals {
enum class CLASSES;
}
namespace CompositeFunctionals {
enum class IMPLEMENTATIONS;
}
struct CUSTOMFUNCTIONAL;

/**
 * @class Functional Functional.h
 * @brief A DFT functional in the way we typically think of it. Typically composed of BASIC_FUNCTIONALS.
 */
class Functional {
 public:
  /**
   * @param impl              The implementation used to evaluate the functional.
   * @param basicFunctionals  The basic functionals making up this composite.
   * @param mixingFactors     The prefactors for each basic functional.
   * @param hfExchangeRatio   The amount exact HF exchange to be added (1.0 = full HF exchange). In case of a
   *                          CAM style functional this parameter is often called alpha. Values other than 0.0
   *                          will mark the functional as a hybrid. Default 0.0.
   * @param hfCorrelRatio    The amount of MP2 correlation to be added. Values other than 0.0 will mark the
   *                          functional as a double-hybrid. Default 0.0.
   * @param lrExchangeRatio   The amount of long-range exchange (HF like) to be added into the evaluation of
   *                          the functional. In CAM type functionals this value is often called beta.
   *                          Values other than 0.0 will mark the functional as a
   *                          range-separated hybrid. Default 0.0.
   * @param mu                The long-range short-range separation parameter mu, present in srDFT and CAM
   *                          style functionals. Default 0.0.
   * @param ssScaling         The same spin scaling factor employed in double hybrid functionals. Default 1.0.
   * @param osScaling         The opposite spin scaling factor employed in double hybrid functionals. Default 1.0.
   */
  Functional(CompositeFunctionals::IMPLEMENTATIONS impl, const std::vector<BasicFunctionals::BASIC_FUNCTIONALS>& basicFunctionals,
             const std::vector<double>& mixingFactors, double hfExchangeRatio = 0.0, double hfCorrelRatio = 0.0,
             double lrExchangeRatio = 0.0, double mu = 0.0, double ssScaling = 1.0, double osScaling = 1.0);
  /// Copy constructor.
  Functional(const Functional&) = default;
  /**
   * @brief Constructor using a CUSTOMFUNCTIONAL object which contains all settings necessary to define a composite
   * functional.
   */
  Functional(CUSTOMFUNCTIONAL customfunc);
  virtual ~Functional() = default;
  /**
   * @returns the underlying basic functionals which are used e.g. in calls to libraries.
   */
  inline const std::vector<BasicFunctionals::BASIC_FUNCTIONALS>& getBasicFunctionals() const {
    return _basicFunctionals;
  }
  /**
   * @returns the number of underlying basic functionals.
   */
  inline unsigned int getNBasicFunctionals() const {
    return _basicFunctionals.size();
  }
  /**
   * @returns the actual functional class to determine what the program needs to calculate
   *          (derivatives, Hartree--Fock exchange, ...)
   */
  inline CompositeFunctionals::CLASSES getFunctionalClass() const {
    return _functionalClass;
  }
  /**
   * @returns the prefactor for the Hartree--Fock exchange to add. Only works for hybrid functionals.
   */
  inline double getHfExchangeRatio() const {
    return _hfExchangeRatio;
  }
  /**
   * @returns the prefactor for the Hartree--Fock correlation (MP2) to add. Only works for double-hybrid functionals.
   */
  inline double getHfCorrelRatio() const {
    return _hfCorrelRatio;
  }
  /**
   * @returns the scaling for same spin contributions of the Hartree--Fock correlation (MP2). Only works for
   * double-hybrid functionals.
   */
  inline double getssScaling() const {
    return _ssScaling;
  }
  /**
   * @returns the scaling for opposite spin contributions of the Hartree--Fock correlation (MP2). Only works for
   * double-hybrid functionals.
   */
  inline double getosScaling() const {
    return _osScaling;
  }
  /**
   * @returns the prefactor for the LR Hartree--Fock exchange to add. Only works for hybrid functionals.
   */
  inline double getLRExchangeRatio() const {
    return _lrExchangeRatio;
  }
  /**
   * @returns the Range separation parameter mu. Only works for hybrid functionals.
   */
  inline double getRangeSeparationParameter() const {
    return _mu;
  }
  /**
   * @returns prefactors for the used BASIC_FUNCTIONALS. For LDAs/GGAs typically just 1.0s, for
   *          hybrids typically different.
   */
  inline std::vector<double> getMixingFactors() const {
    return _mixingFactors;
  }
  /**
   * @returns true if the functional is a hybrid (or meta-hybrid) functional, i.e. whether
   *          Hartree--Fock exchange needs to be calculated.
   */
  inline bool isHybrid() const {
    return _hfExchangeRatio != 0.0;
  }
  /**
   * @returns true if the functional is a double-hybrid functional, i.e. whether
   *          Hartree--Fock exchange and MP2 needs to be calculated.
   */
  inline bool isDoubleHybrid() const {
    return _hfCorrelRatio != 0.0;
  }
  /**
   * @returns true iff the functional is a range separated hybrid (or meta-hybrid) functional, i.e. whether
   *          LR Hartree--Fock exchange needs to be calculated.
   */
  inline bool isRSHybrid() const {
    return _lrExchangeRatio != 0.0;
  }
  /**
   * @returns The implementation to be used to evaluate the functional.
   */
  inline CompositeFunctionals::IMPLEMENTATIONS implementation() const {
    return _impl;
  }

  void print();

 private:
  // Which wrappers can calculate this functional
  CompositeFunctionals::IMPLEMENTATIONS _impl;
  // The list of basic fuctionals that make up the composite
  std::vector<BasicFunctionals::BASIC_FUNCTIONALS> _basicFunctionals;
  // In case some part of the energy is a mixture, e.g. in hybrid functionals.
  std::vector<double> _mixingFactors;
  /**
   * Typically the highest class of the included BASIC_FUNCTIONALS, but may also be hybrid or
   * something else (e.g. range-separated)
   */
  CompositeFunctionals::CLASSES _functionalClass;
  /**
   * For hybrid functionals: A prefactor for the amount of Hartree--Fock exchange that is mixed in.
   */
  double _hfExchangeRatio;
  /**
   * For double hybrid functionals: A prefactor for the amount of Hartree--Fock correlation (MP2) that is mixed in.
   */
  double _hfCorrelRatio;
  /**
   * For range separated hybrid functionals: A prefactor for the amount of Hartree--Fock exchange that is mixed in.
   */
  double _lrExchangeRatio;
  /**
   * For range separated hybrid functionals: The range separation parameter for ther erf
   */
  double _mu;
  /**
   * For double hybrid functionals: Same spin scaling factor for SCS-MP2
   */
  double _ssScaling;
  /**
   * For double hybrid functionals: Opposite spin scaling factor for SCS-MP2 and SOS-MP2
   */
  double _osScaling;
};

inline bool operator==(const Functional& lhs, const Functional& rhs) {
  bool check = (lhs.getFunctionalClass() == rhs.getFunctionalClass());
  check = (check && (lhs.getMixingFactors() == rhs.getMixingFactors()));
  check = (check && (lhs.getBasicFunctionals() == rhs.getBasicFunctionals()));
  check = (check && (lhs.implementation() == rhs.implementation()));
  check = (check && (lhs.getLRExchangeRatio() == rhs.getLRExchangeRatio()));
  check = (check && (lhs.getosScaling() == rhs.getosScaling()));
  check = (check && (lhs.getssScaling() == rhs.getssScaling()));
  check = (check && (lhs.getHfExchangeRatio() == rhs.getHfExchangeRatio()));
  check = (check && (lhs.getHfCorrelRatio() == rhs.getHfCorrelRatio()));
  return check;
}

} /* namespace Serenity */
#endif /* FUNCTIONAL_H */
