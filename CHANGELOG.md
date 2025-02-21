Changelog
=========

Release 1.6.2 (21.02.2025)
--------------------------

### Functionalities

- Designable density functionals in the input by mixing basic functionals (Anton Rikus)
- Stabilized Quasi-Newton Method (SQNM) optimizer for minima (Thorben Wiegmann)
- Added WF in DFT geometry optimizations for HF (Thorben Wiegmann)
- CHELPG and CM5 partial charges (Thorben Wiegmann)
- Approximate embedding electrostatics via partial charges (Lars Hellmann)
- Approximate embedding via the Loewdin expansion of kinetic energy expectation values
  computed from non-orthogonal Slater determinants (Denis G. Artiukhin)
- TDDFT Gradients: restricted/unrestricted, LDA/GGA/hybrid and range-separated functionals, CIS/TDA/TDHF as well,
  RI-J possible, but so far only for an isolated system (Anton Rikus)
- Added tasks for exporting and importing solvation models with corresponding solvent cavities and charges
  to/from files that can be used by Serenity (Lukas Paetow)
- Transition, particle and hole densities can also be plotted for subsystem TDDFT (Anton Rikus)
- Transition and excited state densities from CC2/ADC(2) can be plotted (Anton Rikus)
- Read in external grid potential (Leon Fischer)
- Added a task that writes FCI dump files (Moritz Bensberg).
- Added a task that runs top-down embedding calculations without SCF-based orbital relaxation after subsystem
  partitioning (Moritz Bensberg).
- Added a new embedding flavor, ALMO-MSDFT (Lukas Lampe).

### Technical Features
- Added classes for easier calculation of gradient contributions from two-electron
  integrals using the RI approximation and from one-electron integrals (Anton Rikus)
- Basis files can be directly read in Turbomole format (Anton Rikus)
- LRSCFTaskSettings are written on disk (Anton Rikus)
- Added test systems with converged LRSCF excitation vectors (Anton Rikus)
- Python wrapper published on PyPI, allowing `pip install qcserenity` (Anton Rikus)

Release 1.6.1 (19.03.2024)
--------------------------

### Functionalities

- Added a task that provides direct access to integral files such as the core Hamiltonian (Moritz Bensberg)
- External charges may now be used as an additional potential and read from file (Moritz Bensberg)

### Technical Features

- Negative numbers as input for unsigned variables are now taken as their absolute value and a warning is issued (Niklas Göllmann)
- SCF-Damping reworked internally (Lukas Paetow)

Release 1.6.0 (16.11.2023)
--------------------------

### Functionalities

### Technical Features

- Updated ATOM_SCF initial guess atom densities, now BHLYP in a MINAO basis (Nadim Ramez)
- Removed deprecated ATOM_DENS initial guess (Niklas Niemeyer)
- Moved fractional occupancy keyword to the system block (Niklas Niemeyer)

#### Linear-Response Framework

- CC2/ADC(2) ground- and excited-state densities and dipole moments (Niklas Niemeyer)
- CC2 dynamic polarizabilities and optical rotation (Niklas Niemeyer)
- Triplet excitation energies for CC2/ADC(2) (Niklas Niemeyer)
- Rework Kernel sigmavector (Niklas Niemeyer)
- "Monomer-RI" Coulomb interaction subsystem TDDFT (Niklas Niemeyer)
- TDDFT-ris (one aux. basis function per atom for TDDFT, Niklas Niemeyer)
- Experimental:
  * Coupled CC2/ADC(2) excitation energies, transition moments, excited-state densities and response properties (Niklas Niemeyer)

Release 1.5.3 (25.10.2023)
-------------------------------

### Functionalities

- Added two flavors of restricted open-shell HF and KS for the ground-state (Niklas Niemeyer)
- Fermi-shifted Huzinaga EO Kernel for subsystem TDDFT (Niklas Niemeyer)
- Laplace-Transform GW (Johannes Tölle, Niklas Niemeyer)
- Renamed ReadOrbitalsTask to OrbitalsIOTask (Niklas Göllmann)
- Added the functionality to write Turbomole files (Niklas Göllmann)
- Added the functionality to write Molden files for both spherical and cartesian harmonics (Niklas Göllmann)
- Added three schemes to generate complete basis function products for the Cholesky 
  decomposition framework: Simple, First, Complete (Lars Hellmann) 
- Added the functionality to control density fitting for individual 
  contributions (Coulomb, exchange, long-range exchange, correlation)


Release 1.5.2 (22.03.2023)
-------------------------------

### Functionalities

- Added MOM and IMOM DeltaScf methods (Niklas Niemeyer, Niklas Göllmann)

#### Linear-Response Framework

- Added triplet exctations for TDHF/TDDFT (Niklas Niemeyer)
- Added the following stability analyses for SCF wavefunctions and instability root following (Niklas Niemeyer)
  - Real RHF -> Real RHF
  - Real RHF -> Real UHF
  - Real RHF -> Complex RHF
  - Real UHF -> Real UHF
  - Real UHF -> Complex UHF
- Added spin-flip TDHF/TDDFT (Niklas Niemeyer)

#### Bug Fixes

- Fixed a bug where the T0-correction failed for only 2 electrons.
- Fixed various incorrect settings files in the test resources.
- Fixed an error in FXDTask.cpp.
- Added a factor of one half for the restricted Levelshift potential to be consistent with the other EO potentials
- Serenity is now compilable on macOS, functioning memory management (Apple M1 Pro)


Release 1.5.1 (14.02.2023)
-------------------------------

#### Bug Fixes

- Delete removed libxc functional from Serenity


Release 1.5.0 (13.02.2023)
-------------------------------

### Technical Features

- CMake: changed "native" to "x86-64" as the default option for the march compile flag

#### Dependencies

- Updates the default Libxc library to libxc v6.1.0
- Updates the default ECP library to libecpint v1.0.7
- Updates the default GTest version to v1.13.0
- Updates the default Pybind11 version to v2.10.3
- Allow compilation without any downloads (SERENITY_DOWNLOAD_DEPENDENCIES=OFF)

#### Bug Fixes

- It is now possible to print GEPOL cavities to file.
- Correction to the environmental screening in subsystem-based GW/BSE
- Shifting procedure for not-included orbitals in G0W0/evGW

### Functionalities

#### Linear-Response Framework

- Gauge-origin invariant electronic circular dichroism in the length gauge (Niklas Niemeyer)
- Simplified subsystem TDDFT (Niklas Niemeyer)
- Frozen-virtual, frozen-core and core-only approximations for LR methods (Niklas Niemeyer)
- Interface to the laplace-minimax library (Niklas Niemeyer)
- Laplace-transformation for N4-scaling spin-opposite scaled MP2/ADC(2)/CC2 (Niklas Niemeyer)
- Double-hybrid TDDFT (CIS(D) correction) (Niklas Niemeyer)
- Integral-direct TDDFT sigma vector rework (Niklas Niemeyer)
- Arbitrary combination of couplings (tools/couple.py): FDEc, transition charges, dipole-dipole (Niklas Niemeyer)
- Some performance improvements
  * Adaptive prescreening based on residual norms
  * Exchange and LR-exchange sigmavector contraction symmetry
  * Numerical integration XC potential
  * Numerical integration and kernel contraction
- Experimental:
  * Laplace-transform GW
  * FDEc-BSE calculations possible without TDA

#### General

- The default for implicit solvation is now CPCM instead of IEF-PCM.
- The ReadOrbitalsTask is now able to read Molpro-xml orbital files and
  Molcas-HDF5 orbital files (Moritz Bensberg).
- The ReadOrbitalsTask may now replace the orbital definition in a Molcas-HDF5
  file by Serenity's orbitals (Moritz Bensberg).
- The unrelaxed density is now available for RI-MP2 and DLPNO-MP2 and can 
  be used in embedding calculations (Lukas Lampe).
- Valence virtual orbitals may now be mapped between structures with the DOS algorithm (Moritz Bensberg).
- Valence virtual orbitals may now be localized with the IBO and orbital alignment schemes (Moritz Bensberg).
- The DOS selection threshold may now be optimized automatically to provide a qualitative orbital map (Moritz Bensberg).

Release 1.4.0 (21.10.2021)
-------------------------------

### Functionalities

#### General/Other Features
- SCF convergence thresholds were changed! The new defaults are
  * energy convergence threshold:   5e-8 (old: 1e-8)
  * density convergence threshold:  1e-8 (old: 1e-8)
  * max(FP-PF) threshold:      5e-7 (old: 1e-7)
- Add Broken-Symmetry calculations via KS-DFT and sDFT (Anja Massolle).
- Add a task that orthogonalizes orbitals between subsystems (Anja Massolle).
- The EnergyTask can now evaluate the non-additive kinetic energy contribution
  from orthogonalized subsystem orbitals (Anja Massolle).
- Add ECP gradients (Jan Unsleber).
- Add multi-state FDE Electron Transfer (FDE-ET) and FDE-diab (Patrick Eschenbach).
- Add a task that allows reading of orbitals from other programs.
  Currently, only the ASCII format from turbomole and Serenity's own format are
  supported (Moritz Bensberg).
- Add calculation of quasi-restricted orbitals (Moritz Bensberg).
- Makes Serenity compatible with the MoViPac program (Moritz Bensberg).

##### Local Correlation
- Add occupied orbital partitioning into an arbitrary number of subsystems
  by the generalized direct orbital selection procedure (Moritz Bensberg).
- Add input simplification tasks for local correlation calculations
  (LocalCorrelationTask) and DFT-embedded local correlation calculations
  (DFTEmbeddedLocalCorrelationTask) (Moritz Bensberg).
- Add a task for coupled-cluster-in-coupled-cluster embedding by adjusting
  the DLPNO-thresholds for each region [see JCTC 13, 3198-3207 (2017)]
  (Moritz Bensberg).
- Added a task that allows the fully automatized calculations of relative energies
  form multi-level DLPNO-CC (DOSCCTask) (Moritz Bensberg).
- Core orbitals may be specified in the orbital localization task either by an
  energy cut-off, by tabulated, element-specific numbers, or by explicitly
  giving a number of core orbitals (Moritz Bensberg).

#### Polarizable Continuum Model
- Add a task to calculate the PCM energy contributions for a given
  subsystem density (Jan Unsleber, Moritz Bensberg).
- Add CPCM gradients (Moritz Bensberg).
- Add cavity creation energy calculation from scaled particle
  theory (Moritz Bensberg).
- Changed the default for "minDistance" in the PCM-input block from 0.1 to 0.2.

#### Response Calculations
- Restricted/unrestricted CC2/CIS(Dinf)/ADC(2) excitation energies
  and transition moments from the ground state (Niklas Niemeyer).
- Spin-component and spin-opposite scaled CC2/CIS(Dinf)/ADC(2) (Niklas Niemeyer).
- Quasi-linear and DIIS nonlinear eigenvalue solver (Niklas Niemeyer).
- Natural auxiliary functions (NAFs) for GW/BSE/CC2/CIS(Dinf)/ADC(2) (Niklas Niemeyer).
- Non-orthonormal eigenvalue subspace solver (Niklas Niemeyer).
- Restart system of non-converged eigenpairs in the iterative eigenvalue solvers (Niklas Niemeyer).
- Gauge-origin invariant optical rotation in the length gauge (Niklas Niemeyer).
- Virtual orbital space selection [tested for GW/BSE/TDDFT/TDA/CIS/TDHF/CC2/CIS(Dinf)/ADC(2)/MP2] (Johannes Tölle).
- Diabitazation procedures (multistate FXD, FED, FCD) (Johannes Tölle).
- GW and BSE (with and without environmental screening) (Johannes Tölle).
- Partial response-matrix construction (TDA, TDDFT) (Johannes Tölle, Niklas Niemeyer).
- LibXC support for TDDFT/TDA-Kernel evaluation (Johannes Tölle).
- Mixed exact-approximate embedding schemes for ground and excited states (Johannes Tölle).
- Reimplementation of natural transition orbitals and support for coupled TDDFT (Johannes Tölle).
- Grimme's simplified TDA and TDDFT (Niklas Niemeyer).
- Sigmavector for Exchange contribution using RI, support for long-range exchange and coupled sTDDFT support (Niklas Niemeyer, Johannes Tölle).
- Löwdin transition, hole, and particle charges for response calculations (Anton Rikus, Niklas Niemeyer).
- Transition densities, hole densities, and particle densities can be plotted with the PlotTask (Anton Rikus).
- Natural Response Orbitals can now be plotted (Anton Rikus).

#### Cholesky Decomposition Techniques
- Added Cholesky decomposition techniques (full Cholesky decomposition,
  atomic Cholesky decomposition, atomic-compact Cholesky decomposition) for the evaluation
  of Coulomb and exchange contributions (Lars Hellmann).
- Added atomic and atomic-compact Cholesky basis sets to be used in place of the auxiliary
  basis sets used in the RI formalism (Lars Hellmann).
- Added atomic and atomic-compact Cholesky basis sets to fit integrals in the range-separation
  approach (Lars Hellmann).

#### Electric Fields
- Numerical external electric fields can now be included through point charges arranged in circular 
  capacitor plates around a molecule (Niklas Niemeyer, Patrick Eschenbach).
- Analytical external electric fields and corresponding geometry gradients can now be included through dipole integrals 
  and their derivatives. (Niklas Niemeyer, Patrick Eschenbach).
- Finite-Field Task for (FDE-embedded) numerical and semi-numerical
  calculation of (hyper) polarizabilities (Niklas Niemeyer, Patrick Eschenbach).

### Technical Features

- Update Libecpint to v1.0.4.
- Rework of Libint precision handling.
- Output modifications for simplified handling with MoViPac.
- The MultipoleMomentTask now accepts multiple systems and is able to print
  their total multipole moments.
- The GradientTask may now print the gradient for all atoms in all systems
  in one table.
- Removed outdated keyword "dispersion" from GradientTask, GeometryOptimizationTask
  and HessianTask.
- All basis-set files have been updated to the latest version available
  on www.basissetexchange.org.
- Errors in the def2-series RI MP2 basis sets have been fixed. The old versions were
  actually the MP2 fitting-basis sets of the def-series.
- Rework of DLPNO-MP2/CCSD/CCSD(T).
  Now significantly faster, linear scaling, and caches integrals on disk.
- Fixed an error where the tabulated probe radii for the PCM cavity construction
  where given in Bohr instead of angstrom.
- The Schwarz-prescreening threshold is now by default tied to the basis set size.
  It is calculated as 1e-8/(3M), where M is the number of Cartesian basis
  functions.
- The settings of other tasks may now be forwarded with the block-input system.

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

### Bug Fixes

  - Missing embedding settings in the Python wrapper
  - Generating directories in parallel runs

Release 1.2.1 (27.09.2019)
-------------------------------

### Bug Fixes

- Various smaller Bug Fixes

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
