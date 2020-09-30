Changelog
===============================

Release 1.3.1 (30.09.2020)
-------------------------------

### Technical Features
 - Allow compilation using Clang on both OSX and Linux
 - A few smaller technical bugs
 - Update Libecpint to v1.0.0

Release 1.3.0 (16.09.2020)
-------------------------------

### Functionalities

 - Added SystemSplittingTask and SystemAdditionTask to allow for modular
   system combining and splitting (Moritz Bensberg)
 - Added ElectronicStructureCopyTask to copy the orbitals between systems
   while taking care of displacement and rotation of the molecules
   (only implemented for spherical basis functions) (Moritz Bensberg)
 - Double hybrid functional support for FDE-type calculations (Moritz Bensberg)
 - Off-resonant Response Solver for TDDFT (standard and damped)
   (Niklas Niemeyer)
 - Response Properties from TDDFT (Niklas Niemeyer)
   - Dynamic Polarizabilities (and Linear-Absorption Cross Section)
   - Optical Rotation (and Electronic Circular Dichroism)
 - Added new functionals such as wB97, wB97X, wB97X-D, wB97X-V that became
   available with LibXC (Jan Unsleber)
 - Added x-only and lr-x gradients, enabling range-separated DFT gradient
   calculations (Jan Unsleber)
 - Continuum solvation (IEFPCM, CPCM) is now supported (Moritz Bensberg)
 - DLPNO-based methods are now available (DLPNO-(SCS-)MP2, DLPNO-CCSD(T0))
   (Moritz Bensberg)
 - The direct orbital selection scheme for embedding calculations is now
   available (Moritz Bensberg)
 - DLPNO-MP2 can now be used for double hybrid functionals (Moritz Bensberg)
 - Core and valence orbitals can now be localized independently (Moritz Bensberg)
 - The CubeFileTask is now the PlotTask and can also plot 2D heat-maps (Anja Massolle)

### Technical Features

 - Upgrade XCFun dependency to v2.0.2 (Jan Unsleber)
 - Added option to compile and use LibXC v5.0.0 (Jan Unsleber)
   - Both XCFun and LibXC can be present, default usage is an option at
     compile time.
   - Unittests require XCFun as default.
 - Upgrade Libint2 dependency to v2.7.0.beta6 (Jan Unsleber)
 - Allow linkage of parallel BLAS or Lapack to speed up Eigen3 (Jan Unsleber)
 - Remove `ext/` folder style external projects in favor of CMake submodules
   (Jan Unsleber)
 - XCFun and LibECPint are now cloned from mirrors located publicly at
   https://github.com/qcserenity/xcfun and
   https://github.com/qcserenity/libecpint (Jan Unsleber)
 - Separate evaluation of Coulomb and exchange when using RI
 - Streamlining the keywords used in various embedding tasks by adding
   input-blocks (Moritz Bensberg)
 - Added print-levels to every task (Moritz Bensberg)
 - Energy output files are now encoded as plain ascii files (Moritz Bensberg)
 - Rework of some integral contraction routines (Niklas Niemeyer, Johannes Tölle)
 - Incremental Fock matrix build in the SCF (Johannes Tölle, Moritz Bensberg)
 - Bugfix for range-separate hybrids for Hoffmann and Huzinaga operator
 - Bugfix exact exchange evaluation TDDFT for non-hybrid Nadd-XC
 - Updated density-initial guess files (Patrick Eschenbach).
 - Various smaller technical bugs

Release 1.2.2 (31.10.2019)
-------------------------------
 - Bug Fixes
   - Missing embedding settings in the Python wrapper
   - Generating directories in parallel runs

Release 1.2.1 (27.09.2019)
-------------------------------
 - Bug Fixes

Release 1.2.0 (13.09.2019)
-------------------------------
 - Various small improvements and unit tests
 - TDDFT rework (Michael Boeckers, Johannes Toelle, Niklas Niemeyer)
   - Rework of the eigenvalue solver (Niklas Niemeyer)
   - Rework numerical integration (Johannes Toelle)
   - Sigma Vector rework and RI implementation (Johannes Toelle)
   - Coupled TDDFT calculation with root-following (Michael Boeckers)
   - Exact subsystem TDDFT with root-following (Johannes Toelle,
     Michael Boeckers)
   - Various orbital space selection tools (Johannes Toelle, Niklas Niemeyer)
   - LMO - TDDFT (Johannes Toelle)
   - Rotatory strengths, analytical electric (velocity-gauge) and magnetic
     dipole integrals, manually settable gauge-origin (Niklas Niemeyer)
   - Added unit tests and stability improvements (Johannes Toelle,
     Niklas Niemeyer)
 - Huzinaga/Hoffmann projection operator rework, Fermi-shifted Huzinaga operator
   (Moritz Bensberg)
 - Rework of task input structure (Moritz Bensberg)
 - Speed up basis function in real space evaluation using sparse matrices
   (Moritz Bensberg)
 - Added superposition of atomic potentials as initial guess option
   (Jan Unsleber)

Release 1.1.0 (05.08.2019)
-------------------------------
 - Various small improvements and unit tests
 - Complete rework of the CMake system (Jan Unsleber)
 - Rework Python wrapper to use Pybind11 (Jan Unsleber)
   - Wrapped Loopers (single thread only)
   - Includes automated conversion from Eigen3 to NumPy objects and vice versa
 - Added first unittests for the Python wrapper (Jan Unsleber)
 - Added modified Zhang-Carter reconstruction (David Schnieders)
 - Refactored ProjectionBasedEmbeddingTask to TDEmbeddingTask, now featuring
   reconstruction techniques (David Schnieders, Jan Unsleber)
 - Added Huzinaga and Hoffmann projection operators (Moritz Bensberg)
 - Added various basis truncation techniques (Moritz Bensberg)
 - Added option to relax with respect to precalculated environment in
   TDEmbeddingTask (David Schnieders)
 - Added support for vectors in text input (David Schnieders)
 - Recalculated atom densities for initial guess (David Schnieders)
 - Added ECP support using Libecpint (Moritz Bensberg)
 - Enabled restarts with truncated and extended basis sets (David Schnieders)
 - Added double hybrid functionals (Lars Hellmann)
 - Added SOS/SCS-MP2 (Lars Hellmann)

Release 1.0.0 (29.03.2018)
-------------------------------
 - Various small improvements, bug fixes and code cleaning

Release 1.0.0.RC3 (21.12.2017)
-------------------------------
 - Various small improvements and unit tests
 - Added SAOP model potential (Moritz Bensberg)
 - Added ICC/ICPC 2017 support (Jan Unsleber)
 - Cleaned code for GCC/G++ 6.x and 7.x (Jan Unsleber)
 - Added Intel MKL support via Eigen3 (Jan Unsleber)
 - Added 'diskmode' for data to free memory (Kevin Klahr, Jan Unsleber)
 - Added grid localization using Hilbert R-tree (Jan Unsleber)
 - Added grid block vs. basis function prescreening (Jan Unsleber)
 - Added task for thermal corrections to the energy (Kevin Klahr)
 - Updated XCFun to own branch adding version checks (Moritz Bensberg)
 - Added LLP91/LLP91s,PBE2/PBE2S,PBE3,PBE4,E00 functionals (Moritz Bensberg)
 - Added option for a small supersystem grid in FDE/FaT calculations
   (Jan Unsleber)
 - Added 'SSF' scheme to replace Beckes scheme as default for the grid
   construction (Jan Unsleber)
 - Added 'signed density' to CubeFileTask (Jan Unsleber)
 - Moved energy evaluation to end FDE/FAT in order to allow subsystem grids
   during SCF (David Schnieders, Jan Unsleber)

Release 1.0.0.RC2 (25.09.2017)
-------------------------------
 - Various small improvements and unit tests

Release 1.0.0.RC1 (22.09.2017)
-------------------------------
 - Various bug fixes and unit tests
 - Added Mac OS support (Jan Unsleber)
 - Added Clang/Clang++ support (Jan Unsleber)
 - Added LRSCF (Michael Boeckers)
 - Optimized Grid (Jan Unsleber)
 - Optimized RI Integrals (Jan Unsleber)
 - Rework Grid and MatrixInBasis data objects (Jan Unsleber)
 - Optimized and extended density based initial guesses
   (David Schnieders, Jan Unsleber)
 - Reintroduced final grid (Jan Unsleber)
 - Rotational and translation invariance of gradients and normal modes
   (Kevin Klahr)
 - Switched from libxc to XCFun (Michael Boeckers)
 - Added ADIIS (Jan Unsleber)
 - Added full support for spherical basis functions (Jan Unsleber)
 - Added (R/U)-RI-MP2 (Jan Unsleber)

Beta 0.2.0 (18.03.2017)
-------------------------------
 - Various bug fixes and unit tests
 - Tracking warnings in separate file (David Schnieders)
 - Added HDF5 support for energies (David Schnieders)
 - Rewrote contribution guide (Jan Unsleber)
 - Added semi-numerical Hessian and frequency analysis (Kevin Klahr)
 - Integrated freeze-and-thaw optimization into optimization task (Kevin Klahr)

Beta 0.1.1 (04.02.2017)
-------------------------------
 - Clean GCC-5 warnings (Jan Unsleber)
 - Add basis and .h5 test files to build artifacts (Jan Unsleber)
 - Update to new repository location (Jan Unsleber)
 - Added support for a .pdf manual to CI (Jan Unsleber)
 - Added license to each source file (Jan Unsleber)
 - Added manual draft (Jan Unsleber)

Beta 0.1.0 (03.02.2017)
-------------------------------
(a short list of initial features)

### Electronic Structure Methods
 - Hartree-Fock (incl. gradients)
 - Density Functional Theory (incl. gradients and RI)
 - R-MP2/R-CCSD/R-CCSD(T) (only small systems)

### Initial Guesses
- hCore Guess (zero electron density)
- Extended Hueckel Guess
- Atomic Density Guess (simple version, but works well)

### Embedding Techniques
- FDE (incl. gradients and RI-J)
- Freeze and Thaw (incl. geo. opt.)
- Potential Reconstruction
- Projection Based Embedding (restricted/no basis set truncation)

### Further Tools
- Mulliken Population Analysis
- Orbital Localizations (PM, FB, IBO, EM)
- Output of electron density, MOs and other data in cube format

### Technical
- libint2 will be pre-generated and shipped with Serenity
- libint 2.3.0-beta3 is linked
