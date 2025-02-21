/**
 * @file AtomCenteredBasisController.h
 *
 * @date Jul 30, 2015
 * @author Thomas Dresselhaus
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
#ifndef ATOMCENTEREDBASISCONTROLLER_H
#define ATOMCENTEREDBASISCONTROLLER_H
/* Include Serenity Internal Headers */
#include "basis/BasisController.h"
/* Include Std and External Headers */
#include <memory>
#include <utility>
#include <vector>

namespace Serenity {
/* Forward Declarations */
class Geometry;
class AtomCenteredBasisFactory;
/**
 * @class AtomCenteredBasisController AtomCenteredBasisController.h
 * @brief In most cases atom-centered basis functions are used. This controller manages such a Basis.
 *
 * If basis functions are centered on a certain atom, they are typically associated with that atom.
 * This information is e.g. used in a Mulliken population analysis, or when calculating gradients
 * (with respect to the atomic coordinates).
 */
class AtomCenteredBasisController : public BasisController, public ObjectSensitiveClass<Geometry> {
  friend class BasisController__TEST_SUPPLY;

  struct BasisLoadingData {
    std::string label;
    Eigen::VectorXi shells;
    bool ecp;
  };

 public:
  /**
   * @param geometry with its atoms the basis functions will be associated
   * @param basisLibrary          Path to the basis set library (folder)
   * @param basisLabel            The name of the basis e.g. "6-31Gs"
   * @param makeSphericalBasis    Whether or not to transform the basis in the spherical space.
   * @param makePrimary           This switch will set the basis as the primary one in each atom.
   * @param firstECP              Nuclear charge threshold at which ECPs are used if available.
   */
  AtomCenteredBasisController(std::shared_ptr<Geometry> geometry, const std::string basisLibrary,
                              bool makeSphericalBasis, bool makePrimary, const std::string basisLabel = "",
                              int firstECP = 37, std::vector<Eigen::VectorXi> importantShells = {});
  virtual ~AtomCenteredBasisController() = default;
  /**
   *  @returns the vector of index pairs determining which basis functions belong to an atom. Don't
   *           confuse with getBasisIndicesRed()! Each entry in the vector holds the index pair for
   *           one atom. The first index inside that pair corresponds to the first index in the basis
   *           which belongs to a basis function centered on that atom. The second index is one after
   *           the last index in the basis which belongs to a basis function centered on that atom.
   *           This construction is consistent with what is often used in the standard library, e.g.
   *           the begin() and end() of a std::vector. To loop over basis functions centered on a
   *           particular atom you can construct a loop structure in this way:\n
   *           for (unsigned int mu=indices[atomIndex].first; mu<indices[atomIndex].second; ++mu)
   *  @see Basis
   */
  inline const std::vector<std::pair<unsigned int, unsigned int>>& getBasisIndices() {
    if (!_basis)
      produceBasis();
    return _basisIndicesOfAtom;
  }
  /**
   *  @returns the vector of index pairs determining which basis functions belong to an atom,
   *           counting each shell of functions as one (e.g. a set of px, py and pz increases the
   *           indices by 1 instead of 3). In all other respects this is the same as getBasisIndices().
   *  @see Basis
   */
  inline const std::vector<std::pair<unsigned int, unsigned int>>& getBasisIndicesRed() {
    if (!_basis)
      produceBasis();
    return _basisIndicesRedOfAtom;
  }
  /**
   * @brief Getter for the mapping of basis function to atom index.
   * @return A std vector of size nBasis containing the atom index that the basis function belongs to.
   */
  inline const std::vector<unsigned int>& getAtomIndicesOfBasis() {
    if (!_basis)
      produceBasis();
    return _atomIndicesOfBasisFunctions;
  }

  /**
   * @brief Getter for the mapping of basis function to atom index, but as a projection matrix (nAtoms x nBasis). It
   * contains a single 1 in each column (with the row indicating the atom on which the basis function is centered).
   * @return The (nAtoms x nBasis) projection matrix.
   */
  inline const Eigen::MatrixXd& getAtomBasisProjection() {
    if (!_basis)
      produceBasis();
    return _atomBasisProjection;
  }
  /**
   * @brief Getter for the mapping of basis shell to atom index.
   * @return A std vector of size nShells containing the atom index that the shell belongs to.
   */
  inline const std::vector<unsigned int>& getAtomIndicesOfBasisShells() {
    if (!_basis)
      produceBasis();
    return _atomIndicesOfShells;
  }
  /**
   *  @returns the name of the basis
   *  @see Basis
   */
  inline const std::string getBasisLabel() const {
    return _basisLabel;
  }

  /**
   * @return The geometry with respect to which the AtomCenteredBasisController is defined.
   */
  inline const Geometry& getGeometry() {
    return *_geometry;
  }

  inline const std::string& getBasisLibrary() {
    return _basisLibrary;
  }

  inline bool getMakeSphericalBasis() {
    return _makeSphericalBasis;
  }

  inline bool getMakePrimary() {
    return _makePrimary;
  }

  void toHDF5(std::string fBaseName, std::string id);

  void fromHDF5(std::string fBaseName, std::string id);

  std::vector<BasisLoadingData> getBasisLoadingData() {
    return _basisLoadingData;
  }

  void notify() override;

 private:
  std::unique_ptr<Basis> produceBasisFunctionVector() override final;
  void postConstruction() override final;

  std::shared_ptr<Geometry> _geometry;
  const std::string _basisLibrary;
  const std::string _basisLabel;
  bool _makeSphericalBasis;
  bool _makePrimary;
  std::vector<BasisLoadingData> _basisLoadingData;
  int _firstECP;

  std::vector<std::pair<unsigned int, unsigned int>> _basisIndicesOfAtom;
  std::vector<std::pair<unsigned int, unsigned int>> _basisIndicesRedOfAtom;
  std::vector<unsigned int> _atomIndicesOfBasisFunctions;
  Eigen::MatrixXd _atomBasisProjection;
  std::vector<unsigned int> _atomIndicesOfShells;

  void resetBasisLoadingData(const std::vector<Eigen::VectorXi>& importantShells = {});

  /*
   * This private constructor is for testing purposes only
   */
  AtomCenteredBasisController() : BasisController("TestBasis") {
  }
};

} /* namespace Serenity */
#endif /* ATOMCENTEREDBASISCONTROLLER_H */
