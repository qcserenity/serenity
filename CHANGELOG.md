Changelog
===============================

Release 1.3.0 (under development)
-------------------------------
 - Added SCINE Interface (Jan Unsleber)
 - Added MCSCF solver(s) (Jan Unsleber, Stefan Knecht)
 - Added DMRG-CI and DMRG-SCF capabilities via a QCMaquis Interface (Jan Unsleber)

Release 1.2.0 (13.09.2019)
-------------------------------
 - Various small improvements and unit tests
 - TDDFT rework (Michael Boeckers, Johannes Toelle, Niklas Niemeyer)
   - Rework of the eigenvalue solver (Niklas Niemeyer)
   - Rework numerical integration (Johannes Toelle)
   - Sigma Vector rework and RI implementation (Johannes Toelle)
   - Coupled TDDFT calculation with root-following (Michael Boeckers)
   - Exact subsystem TDDFT with root-following (Johannes Toelle, Michael Boeckers)
   - Various orbital space selection tools (Johannes Toelle, Niklas Niemeyer)
   - LMO - TDDFT (Johannes Toelle)
   - Rotatory strengths, analytical electric (velocity-gauge) and magnetic dipole integrals, manually settable gauge-origin (Niklas Niemeyer)
   - Added unit tests and stability improvements (Johannes Toelle, Niklas Niemeyer)
 - Huzinaga/Hoffmann projection operator rework, Fermi-shifted Huzinaga operator (Moritz Bensberg)
 - Rework of task input structure (Moritz Bensberg)
 - Speed up basis function in real space evaluation using sparse matrices (Moritz Bensberg)
 - Added superposition of atomic potentials as initial guess option (Jan Unsleber)

Release 1.1.0 (05.08.2019)
-------------------------------
 - Various small improvements and unit tests
 - Complete rework of the CMake system (Jan Unsleber)
 - Rework Python wrapper to use Pybind11 (Jan Unsleber)
   - Wrapped Loopers (single thread only)
   - Includes automated conversion from Eigen3 to NumPy objects and vice versa
 - Added first unittests for the Python wrapper (Jan Unsleber)
 - Added modified Zhang-Carter reconstruction (David Schnieders)
 - Refactored ProjectionBasedEmbeddingTask to TDEmbeddingTask, now featuring reconstruction techniques (David Schnieders, Jan Unsleber)
 - Added Huzinaga and Hoffmann projection operators (Moritz Bensberg)
 - Added various basis truncation techniques (Moritz Bensberg)
 - Added option to relax with respect to precalculated environment in TDEmbeddingTask (David Schnieders)
 - Added support for vectors in text input (David Schnieders)
 - Recalculated atom densities for initial guess (David Schnieders)
 - Added ECP support using Libecpint (Moritz Bensberg)
 - Enabled restarts with truncated and extended basis sets (David Schnieders)

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
 - Added task for thermal corections to the energy (Kevin Klahr)
 - Updated XCFun to own branch adding verison checks (Moritz Bensberg)
 - Added LLP91/LLP91s,PBE2/PBE2S,PBE3,PBE4,E00 functionals (Moritz Bensberg)
 - Added option for a small supersystem grid in FDE/FaT calculations (Jan Unsleber)
 - Added 'SSF' scheme to replace Beckes scheme as default for the grid construction (Jan Unsleber)
 - Added 'signed density' to CubeFileTask (Jan Unsleber)
 - Moved energy evaluation to end FDE/FAT in order to allow subsystem grids during SCF (David Schnieders, Jan Unsleber)
 
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
 - Optimized and extended density based initial guesses (David Schnieders, Jan Unsleber)
 - Reintroduced final grid (Jan Unsleber)
 - Rotational and translation invariance of gradients and normal modes (Kevin Klahr)
 - Switched from libxc to XCFun (Michael Boeckers)
 - Added ADIIS (Jan Unsleber)
 - Added full support for spherical basis functions (Jan Unsleber)
 - Added (R/U)-RI-MP2 (Jan Unsleber)

Beta 0.2.0 (18.03.2017)
-------------------------------
 - Various bug fixes and unit tests
 - Tracking warnings in seperate file (David Schnieders)
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
- libint2 will be pre-genrated and shiped with Serenity
- libint 2.3.0-beta3 is linked 
