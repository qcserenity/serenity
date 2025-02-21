/**
 * @file   FCIDumpFileWriter.h
 *
 * @date   Jan 30, 2024
 * @author Moritz Bensberg
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

#ifndef SERENITY_FCIDUMPFILEWRITER_H
#define SERENITY_FCIDUMPFILEWRITER_H

/* Include Serenity Internal Headers */
#include "data/SpinPolarizedData.h"
#include "data/matrices/DensityMatrix.h"
#include "data/matrices/FockMatrix.h"
#include "tasks/FCIDumpFileWriterTask.h"
/* Include Std and External Headers */
#include <Eigen/Dense>
#include <memory>
#include <ostream>
#include <vector>

namespace Serenity {
class SystemController;
struct LocalCorrelationSettings;
class OrbitalPair;

/**
 * @class FCIDumpFileWriter FCIDumpFileWriter.h
 * @tparam SCFMode Template for the RESTRICTED or UNRESTRICTED implementation
 * @brief This class writes MOLPRO FCI dump files. This class uses the RI approximation.
 * The integral transformation is based on the local correlation code.
 *
 * An FCI dump file provides all two electron and one electron integrals in MO basis.
 * The orbital index counting starts at 1. The third and fourth index for the one
 * particle integrals are zero. The file start with (exactly) 4 lines of header, providing
 * the number of orbitals, electrons, the number of excess alpha electron (MS2), orbital
 * symmetry, and system symmetry. After the header the integrals are encoded through their
 * value and their orbital indices.
 *
 * A typical FCI dump file may look like this:
 * &FCI NORB= 8, NELEC= 10, MS2= 0,
 * ORBSYM= 1, 1, 1, 1, 1, 1, 1, 1,
 * ISYM= 1
 * &END
 *   8.2371169364e-01     1     1     1     1
 *   6.5011589197e-02     2     1     1     1
 *   6.5841092835e-01     2     2     1     1
 *   1.0997174278e-01     2     1     2     1
 *  -6.9437938103e-02     2     1     2     2
 *  ...
 *   7.7071181164e-01     8     6     0     0
 *  -6.6624329107e-16     8     7     0     0
 *  -6.2473727922e+00     8     8     0     0
 *
 */
template<Options::SCF_MODES SCFMode>
class FCIDumpFileWriter {
 public:
  /**
   * @brief Constructor
   * @param activeSystem          The active system controller.
   * @param environmentSystems    Environment systems which may contribute to the one particle integrals.
   * @param settings              The task settings.
   * @param oneParticleIntegrals  If already present, the one particle AO integrals can be provided to avoid
   * recalculating them.
   */
  FCIDumpFileWriter(std::shared_ptr<SystemController> activeSystem,
                    std::vector<std::shared_ptr<SystemController>> environmentSystems,
                    const FCIDumpFileWriterTaskSettings& settings,
                    std::shared_ptr<FockMatrix<SCFMode>> oneParticleIntegrals = nullptr);
  /**
   * @brief Write the FCI dump file to the output stream.
   * @param outputStream The output stream.
   * @param orbitalRange The orbital ranges for which to write the file.
   */
  void writeFCIDumpFile(std::ostream& outputStream, const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);
  /**
   * @brief Getter for the one particle integrals (everything but the active space two electorn contribution to the Fock
   * operator).
   * @param orbitalRange The orbital range defining the active space.
   * @return The one particle integrals in AO basis.
   */
  FockMatrix<SCFMode> getOneParticleIntegrals(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);
  /**
   * @brief Getter for the total uncorrelated energy.
   * @return The total energy.
   */
  double getTotalUncorrelatedEnergy();
  /**
   * @brief Getter for the uncorrelated energy of the active space.
   * @param orbitalRange The orbital range defining the active space.
   * @return The energy.
   */
  double getUncorrelatedActiveSpaceEnergy(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);
  /**
   * @brief Getter for the uncorrelated energy of the core/inactive orbitals.
   * @param orbitalRange The orbital range defining the active space.
   * @return The energy.
   */
  double getTotalCoreEnergy(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);

 private:
  /**
   * @brief Write the one particle part.
   */
  void writeCorePart(std::ostream& outputStream, const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);
  /**
   * @brief Write the two electron part.
   */
  void writeERIPart(std::ostream& outputStream, const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);
  /**
   * @brief Write the file header.
   */
  void writeHeader(std::ostream& outputStream, const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);
  /**
   * @brief Write the energy contribution from all orbitals not in the active space/orbital ranges.
   */
  void writeCoreEnergy(std::ostream& outputStream);
  /**
   * @brief Get the one particle integrals in the MO basis.
   * @return The one particle integrals in MO basis.
   */
  SPMatrix<SCFMode> getCore(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);
  /**
   * @brief Build the prescreening map for significant orbitals pairs ik and jk in the integral (ik|jl) using
   * differential overlap prescreening.
   * @param coefficients The orbital coefficients.
   * @return A sparse matrix for every orbital-orbital combination with a non-zero entry if the pair ik has significant
   * diff. overlap.
   */
  SparseMap getOrbital2OrbitalMap(std::shared_ptr<Eigen::MatrixXd> coefficients);
  /**
   * @brief Construct the orbital pair object which will be used to run the integral transformation.
   * @param orb2OrbMap The orbital to orbital map. Not that all pairs ij will be constructed since the integral (ik|jl)
   * does not vanish if ij have vanishing differential overlap.
   * @return The orbital pairs.
   */
  std::vector<std::shared_ptr<OrbitalPair>> getPseudoOrbitalPairs(std::shared_ptr<SparseMap> orb2OrbMap);
  /**
   * @brief Getter for the (embedded) Fock matrix.
   * @return The Fock matrix.
   */
  FockMatrix<SCFMode> getFockMatrix();
  /**
   * @brief Calculate the one particle integrals.
   */
  void updateOneParticleIntegrals(const SpinPolarizedData<SCFMode, std::vector<unsigned int>>& orbitalRange);

  ///@brief The settings.
  const FCIDumpFileWriterTaskSettings& _settings;
  ///@brief The active system controller.
  std::shared_ptr<SystemController> _activeSystem;
  ///@brief The environment system controller (if any).
  std::vector<std::shared_ptr<SystemController>> _environmentSystems;
  ///@brief The one particle integrals in AO basis.
  std::unique_ptr<FockMatrix<SCFMode>> _oneParticleIntegrals;
  ///@brief The energy of the inactive orbitals.
  std::unique_ptr<double> _coreEnergy = std::make_unique<double>(0.0);
  ///@brief The energy of the active orbitals evaluated with the uncorrelated density.
  std::unique_ptr<double> _totalInitialActiveSpaceEnergy = std::make_unique<double>(0.0);
  ///@brief Total uncorrelated energy.
  std::unique_ptr<double> _totalUncorrelatedEnergy = std::make_unique<double>(0.0);
  ///@brief The orbital range of the last active space.
  std::unique_ptr<SpinPolarizedData<SCFMode, std::vector<unsigned int>>> _lastActiveSpace;
};

} // namespace Serenity

#endif // SERENITY_FCIDUMPFILEWRITER_H
