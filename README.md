![](doc/SerenityLogo.png)

**Table of Contents** 
- [Synopsis](#synopsis)
- [License and Copyright Info](#license-and-copyright-info)
- [Download](#download)  
- [Install](#install)
    - [Prerequisites](#prerequisites)
    - [Install Using CMake and Make](#install-using-cmake-and-make)
    - [Install Including Python Interface](#install-including-python-interface)
    - [Tests](#tests)
    - [Documentation](#documentation)
    - [Troubleshooting](#troubleshooting) 
- [User Manual](#user-manual)
- [How to Cite Serenity](#how-to-cite-serenity)
- [Contact](#contact)
    - [Bugs and Feature Requests](#bugs-and-feature-requests)
    - [Other](#other)

## Synopsis
Serenity is a quantum chemistry code originally
developed in the workgroup of Johannes Neugebauer
at the WWU Münster.       
Serenity has a strong focus on quantum chemical subsystem/embedding methods.

## License and Copyright Info

Serenity is free software: you can redistribute it and/or modify
it under the terms of the LGNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the LGNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Download 

In order to download the source files of the 
latest Serenity release please visit:  
https://github.com/qcserenity/serenity

## Install

Please read the following instructions carefully.

### Prerequisites
The code has been tested and compiled with GCC/G++ (Versions 7.X and 8.X)
and ICC/ICPC (Versions 17.X) compilation with other compilers and complier versions might still
be possible though.

The compilation of the code is supported on Linux systems, and should also be possible on MacOS.
However, the latter will only be fully supporten in future releases.

The following programs/libraries must be available on your system:
 - CMake (Version >= 3.9)
 - Boost (for package managers: including boost-devel)
 - OpenMP
 - Eigen3 
 - HDF5 (Version >= 1.10.1; including header files and cmake files)
 - A recent GMP version, including C++ support (for libint2)
 - The standard GNU toolchain (make, tar, autoconf, libtool)  
 
The following libraries may be preinstalled and will otherwise be automatically 
downloaded and installed locally:
 - xcfun (The Serenity verison is forked to: https://github.com/moritzBens/xcfun.git)  
 - libecpint (The Serenity version is forked to: https://github.com/moritzBens/libecpint.git)  
The following libraries will be downloaded and installed locally regardless of their earlier
presence:
 - GTest (Google Test and Google Mock)
 - libint2 (Version 2.2.0-beta3, pre-configured and hosted at: https://thclab.uni-muenster.de/serenity/ext-libint)

The following libraries are optional and needed for additional features:
 - Intel MKL (for SMP parallel Eigen3 eigenvalue solvers)
 - Doxygen (for the documentation)
 - Python-devel (for the python wrapper)   
 - pybind11 (for the python wrapper)

### Install Using CMake and Make
Extract or pull the source code, then create a build directory:
> cd serenity
> mkdir build  
> cd build  

Then run cmake:
> cmake ..    
(If the build folder is not located inside the main directory of Serenity 
please adapt the path accordingly.)

Finally run make and make install to build the program:
> make         

Please source serenity.sh located in the main folder to set all necessary environment 
variables:
> cd ..  
> source serenity.sh

### Install Including Python Interface
In order to activate the compilation of a Python interface to the code
a flag can be set as follows
> cmake -DSERENITY_PYTHON_BINDINGS=ON ..

additionally this and other flags can be toggled using `ccmake`
> ccmake ..

The wrapper will be build for the Python version that CMake finds first
In order to point CMake to a specific version of Python the following 
option can be used:
> cmake -DPYTHON_EXECUTABLE=/usr/bin/python3 ..

The wrapper is shipped in form of a shared library (serenipy.so).
In order for Python to find this package the library folder has to be present 
in the `PYTHONPATH` environment variable and the Serenity library has to be 
present in a path searchd by the system for shared libraries.
The latter is done when sourceing the `serenith.sh` script, the former requires 
to uncomment one line in this file.
Afterwards the interface should importable as follows:
> python  
> import serenipy as spy  

### Tests
Serenity comes with a decent set of unittests in order to run them source 
the `serenity.sh` script and run
> PATH_TO_BUILD_DIR/bin/serenity_tests

Please use the complete path, or else GTEST might run into problems.  
  
The Python wrapper comes with its own set of tests, these can be run using
> python -m unittest discover PATH_TO_SERENITY/src/python/tests

### Documentation
After configuring the project using CMake it is possible to create the documentation 
using:
> make doc  

Then you can open up the doc/html/index.html in a browser.  
The python wrapper is not featured inside the documentation.
Its documentation is available via pythons help() function,
which displays the build-in doc-strings.

### Troubleshooting:
If you run into trouble during the compilation please make sure all
required libraries are present and read the output carfully.
If you can make sure that your problems are not due to lacking requirements
feel free to contact the main developers or preferably open an issue in the 
repository.

## User Manual
The user manual resides in a seperate folder (manual) in this repository.

## How to Cite Serenity
Serenity is published as:  
  
J. P. Unsleber, T. Dresselhaus, K. Klahr, D. Schnieders, M. Böckers, D. Barton and J. Neugebauer,   
Serenity: A Subsystem Quantum Chemistry Program,   
*J. Comput. Chem.*, **39**, 788--798, (2018).  
  
The BibTeX code would thus be:  
>@article{serenity_pub,  
> title = {Serenity: A Subsystem Quantum Chemistry Program},  
> author = {Jan P. Unsleber and Thomas Dresselhaus and Kevin Klahr
>             and David Schnieders and Michael B{\"o}ckers and
>             Dennis Barton and Johannes Neugebauer},  
> journal = {J. Comput. Chem.},  
> volume = {39},  
> pages = {788--798},  
> year = {2018}  
>}  
  
In order to allow others to reproduce your data please also reference the Verison of Serenity used like this:  
  
Serenity Version: 1.1.0, https://github.com/qcserenity/serenity (2019)  
  
or as BibTeX code:  
>@misc{serenity_version,  
>  title = {Serenity},  
>  notes = {Serenity Version: 1.0.0, \url{{https://github.com/qcserenity/serenity}}},  
>  year  = {2019}  
>  }  

## Contact

### Bugs and Feature Requests
For both bugs and feature requests please use the issue tracker.
For specific parts of the code feel free to tag the one or more of the 
developers that contributed to it (as seen in the doxygen header).

### Other
For other question, requests or simply to give some feedback feel free to send an e-mail
to: 
serenity@uni-muenster.de
