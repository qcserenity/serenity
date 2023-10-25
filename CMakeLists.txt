cmake_minimum_required(VERSION 3.9)

project(serenity
        VERSION 1.5.3
        DESCRIPTION "Serenity: A subsystem quantum chemistry program.")

set(CMAKE_MODULE_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

#########################
## Buildtype and Flags ##
#########################

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type specified. Will build Release")
  set(CMAKE_BUILD_TYPE "RELEASE")
endif(NOT CMAKE_BUILD_TYPE)


# Optional flags
set(SERENITY_MARCH "x86-64" CACHE STRING "Compile with -march=<FLAG>")
option(WERROR "Compile with warnings as errors" OFF)
option(GCC_PROFILE "Compile with profile flags" OFF)
option(GCC_COVERAGE "Compile with coverage flags" OFF)
if (GCC_COVERAGE)
  set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
endif()
option(SERENITY_PREFER_XCFUN "If both XCFun and LibXC are present, prefer XCFun per default" ON)
option(SERENITY_USE_XCFUN "Include XCFun as DFT functional library" ON)
option(SERENITY_USE_LIBXC "Include LibXC as DFT functional library" ON)
option(SERENITY_USE_LAPLACE_MINIMAX "Do not include laplace-minimax by default" OFF)
if ((NOT SERENITY_USE_LIBXC) AND (NOT SERENITY_USE_XCFUN))
  message(FATAL_ERROR "At least one of SERENITY_USE_LIBXC and SERENITY_USE_XCFUN has to be enabled.")
endif()
option(SERENITY_DOWNLOAD_DEPENDENCIES "If OFF, will require the operator to provide all dependencies." ON)

# Get some architecture info
get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)
set(ARCH_LIB_PATH "lib")
if ("${LIB64}" STREQUAL "TRUE")
  set(ARCH_LIB_PATH "lib64")
endif()


##########################################################
## Set up lists of directories that need to be included ##
##########################################################

# Create the dummy task and dummy input file if not present
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/src/tasks/DummyTask.h)
  file(COPY ${PROJECT_SOURCE_DIR}/dev/templates/DummyTask.h DESTINATION ${PROJECT_SOURCE_DIR}/src/tasks/)
endif()
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/src/tasks/DummyTask.cpp)
  file(COPY ${PROJECT_SOURCE_DIR}/dev/templates/DummyTask.cpp DESTINATION ${PROJECT_SOURCE_DIR}/src/tasks/)
endif()
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/test.input)
  file(COPY ${PROJECT_SOURCE_DIR}/dev/templates/test.input DESTINATION ${PROJECT_SOURCE_DIR}/)
endif()

# Set up lists of files
include(src/Files.cmake)

##################
## Main Targets ##
##################

add_library(serenity ${SERENITY_CPPS} ${SERENITY_HEADERS})
set_property(TARGET serenity PROPERTY POSITION_INDEPENDENT_CODE ON)
add_executable(serenity_exe ${PROJECT_SOURCE_DIR}/src/serenity.cpp)

##################
## Dependencies ##
##################
message("#=========================#")
message("|  External Dependencies  |")
message("#=========================#")

# Eigen3
# Intel MKL
message("-----Checking Eigen3-----")
message("-----Checking Intel MKL-----")
include(AddEigen)
add_eigen(serenity PUBLIC)

# OpenMP
message("-----Checking OpenMP-----")
find_package(OpenMP REQUIRED)

# Laplace Minimax
if (SERENITY_USE_LAPLACE_MINIMAX)
  message("-----Checking Laplace_Minimax-----")
  if (SERENITY_DOWNLOAD_DEPENDENCIES)
    include(ImportLaplaceMinimax)
    import_laplace_minimax()
  else()
    message(FATAL_ERROR "Laplace-Minimax cannot be detected via CMake please enable download or disable this dependency")
  endif()
endif(SERENITY_USE_LAPLACE_MINIMAX)

# Libint2
message("-----Checking Libint2-----")
if (SERENITY_DOWNLOAD_DEPENDENCIES)
  include(ImportLibint)
  import_libint()
else()
  find_package(Libint2 CONFIG REQUIRED)
  add_library(libint2-static INTERFACE IMPORTED)
  target_link_libraries(libint2-static INTERFACE Libint2::libint2)
endif()

if (SERENITY_USE_XCFUN)
  # xcfun
  message("-----Checking XCFun----")
  if (SERENITY_DOWNLOAD_DEPENDENCIES)
    include(ImportXCFun)
    import_xcfun()
  else()
    find_package(XCFun CONFIG REQUIRED)
    add_library(xcfun INTERFACE IMPORTED)
    target_link_libraries(xcfun INTERFACE XCFun::xcfun)
  endif()
endif()

if (SERENITY_USE_LIBXC)
  # libxc
  message("-----Checking Libxc----")
  if (SERENITY_DOWNLOAD_DEPENDENCIES)
    include(ImportLibxc)
    import_libxc()
  else()
    find_package(PkgConfig QUIET)
    pkg_check_modules(pc_libxc libxc)
    if(pc_libxc_FOUND)
      add_library(xc INTERFACE IMPORTED)
      target_link_libraries(xc INTERFACE ${pc_libxc_LINK_LIBRARIES})
      target_include_directories(xc INTERFACE ${pc_libxc_INCLUDE_DIRS})
    endif()
  endif()
endif()

# HDF5
message("-----Checking HDF5-----")
if (NOT DEFINED HDF5_USE_STATIC_LIBRARIES)
  set(HDF5_USE_STATIC_LIBRARIES OFF)
endif()
find_package(HDF5 REQUIRED COMPONENTS HL CXX)
# The following if is a small helper for systems
#  that have HDF5 installed but do not have the
#  cmake files installed (this seems common for
#  many distributions)
if (HDF5_USE_STATIC_LIBRARIES)
  set(HDF5_TARGET hdf5_cpp-static)
else()
  set(HDF5_TARGET hdf5_cpp-shared)
endif()
if(NOT TARGET ${HDF5_TARGET})
  add_library(${HDF5_TARGET} INTERFACE IMPORTED)
  set_target_properties(${HDF5_TARGET} PROPERTIES
    INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
  )
endif()

# Ensure that the include directory set by HDF5 contains hdf5.h
get_target_property(_hdf5_include_dirs ${HDF5_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
set(_has_hdf5_h FALSE)
foreach(_hdf5_include_dir ${_hdf5_include_dirs})
  if(EXISTS ${_hdf5_include_dir}/hdf5.h)
    set(_has_hdf5_h TRUE)
    break()
  endif()

  # Try looking for a hdf5 directory within the include directory
  if(EXISTS ${_hdf5_include_dir}/hdf5/hdf5.h)
    set_target_properties(${HDF5_TARGET} PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${_hdf5_include_dir}/hdf5
    )
    set(_has_hdf5_h TRUE)
    break()
  endif()
endforeach()
if(NOT _has_hdf5_h)
  message(FATAL_ERROR "Could not find hdf5.h in HDF's include directorie: ${_hdf5_include_dirs}")
endif()
unset(_has_hdf5_h)
unset(_hdf5_include_dirs)

# Boost
message("-----Checking Boost-----")
find_package(Boost REQUIRED)

# Libecpint
message("-----Checking Libecpint-----")
if (SERENITY_DOWNLOAD_DEPENDENCIES)
  include(ImportLibecpint)
  import_libecpint()
else()
  find_package(ecpint CONFIG REQUIRED)
  add_library(ecpint INTERFACE IMPORTED)
endif()

#################
## Definitions ##
#################
find_package(Git)
if(Git_FOUND)
  message(STATUS "Found git executable: ${GIT_EXECUTABLE}")
  # In case no git repository is found, only redirect stderr to /dev/null.
  execute_process ( COMMAND bash -c "git rev-parse --abbrev-ref HEAD 2> /dev/null" OUTPUT_VARIABLE GIT_BRANCH )
  if (GIT_BRANCH STREQUAL "")
    message(STATUS "No git repository found.")
    set(GIT_BRANCH "UNKNOWN")
  else ()
    # Remove new line character.
    string ( STRIP ${GIT_BRANCH} GIT_BRANCH )
    message(STATUS "Found git repository: " ${GIT_BRANCH})
  endif()
else ()
  set(GIT_BRANCH "UNKNOWN")
endif()
add_definitions (-DGIT_BRANCH="${GIT_BRANCH}")
# configure_file (${PROJECT_SOURCE_DIR}/src/io/FormattedOutput.cpp ${PROJECT_SOURCE_DIR}/src/io/FormattedOutput.cpp)

#################################
## Add dependencies to targets ##
#################################

# Shared library
target_include_directories(serenity PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/serenity>
  $<$<BOOL:${MKL_FOUND}>:${MKL_INCLUDE_DIRS}>
)
target_compile_options(serenity PRIVATE ${OpenMP_CXX_FLAGS})
set_target_properties(serenity PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  CXX_STANDARD 14
)
target_link_libraries(serenity
  PUBLIC
    Boost::boost
    OpenMP::OpenMP_CXX
    $<$<BOOL:${GCC_COVERAGE}>:-lgcov>
    $<$<BOOL:${GCC_PROFILE}>:-pg>
    libint2-static
  PRIVATE
    xc
    xcfun
    ecpint
    $<$<BOOL:${SERENITY_USE_LAPLACE_MINIMAX}>:laplace-minimax>
    ${HDF5_TARGET}
)
target_compile_definitions(serenity PUBLIC
  $<$<BOOL:${SERENITY_PREFER_XCFUN}>:SERENITY_PREFER_XCFUN>
  $<$<BOOL:${SERENITY_USE_XCFUN}>:SERENITY_USE_XCFUN>
  $<$<BOOL:${SERENITY_USE_LIBXC}>:SERENITY_USE_LIBXC>
  $<$<BOOL:${SERENITY_USE_LAPLACE_MINIMAX}>:SERENITY_USE_LAPLACE_MINIMAX>
  $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
  $<$<BOOL:${MKL_FOUND}>:EIGEN_USE_MKL_ALL>
)
target_compile_options(serenity PUBLIC
  -Wall
  -Wextra
  -Wno-comment
  $<$<BOOL:${WERROR}>:-Werror>
  $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
  $<$<BOOL:${GCC_PROFILE}>:-pg>
  $<$<AND:$<BOOL:${MKL_FOUND}>,$<STREQUAL:"${CMAKE_CXX_COMPILER_ID}","Intel">>:-DMKL_LP64>
)
if(NOT "${SERENITY_MARCH}" STREQUAL "" AND NOT MSVC)
  target_compile_options(serenity PUBLIC -march=${SERENITY_MARCH})
endif()


# Executable
set_target_properties(serenity_exe PROPERTIES OUTPUT_NAME serenity)
target_link_libraries(serenity_exe
  PUBLIC
    $<$<BOOL:${GCC_COVERAGE}>:-lgcov>
    $<$<BOOL:${GCC_PROFILE}>:-pg>
  PRIVATE
    serenity
)
set_target_properties(serenity_exe PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  CXX_STANDARD 14
)
target_compile_options(serenity_exe
  PUBLIC
    -Wall
    -Wextra
    -Wno-comment
    $<$<BOOL:${WERROR}>:-Werror>
    $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
    $<$<BOOL:${GCC_PROFILE}>:-pg>
)
if(NOT "${SERENITY_MARCH}" STREQUAL "" AND NOT MSVC)
  target_compile_options(serenity_exe PUBLIC -march=${SERENITY_MARCH})
endif()

#############
## Install ##
#############

# Targets
install(TARGETS serenity serenity_exe
        EXPORT serenityTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/
  DESTINATION include/serenity
  FILES_MATCHING PATTERN "*.h"
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/basis
        ${CMAKE_CURRENT_SOURCE_DIR}/data/initialGuess
        DESTINATION share/serenity/data
)



###############
## Unittests ##
###############

option(SERENITY_ENABLE_TESTS "Compiles Serenity's unittests." ON)
if(SERENITY_ENABLE_TESTS)
  message("#================#")
  message("| Test Framework |")
  message("#================#")
  enable_testing()
  message("-----Checking GTest-----")
  # Load gtest dependency
  if (SERENITY_DOWNLOAD_DEPENDENCIES)
    include(ImportGTest)
    import_gtest()
  else()
    find_package(GTest REQUIRED)
  endif()
  message("-----Adding Test Executable-----")
  # Test exe
  add_executable(serenity_tests ${SERENITY_TEST_FILES})
  target_link_libraries(serenity_tests
    PUBLIC
      $<$<BOOL:${GCC_COVERAGE}>:-lgcov>
      $<$<BOOL:${GCC_PROFILE}>:-pg>
    PRIVATE
      GTest::Main
      GMock::GMock
      serenity
      xcfun
      ecpint
      ${HDF5_TARGET}
  )
  set_target_properties(serenity_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    CXX_STANDARD 14
  )
  target_compile_options(serenity_tests
    PUBLIC
      -Wall
      -Wextra
      -Wno-comment
      $<$<BOOL:${WERROR}>:-Werror>
      $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
      $<$<BOOL:${GCC_PROFILE}>:-pg>
  )
  if(NOT "${SERENITY_MARCH}" STREQUAL "" AND NOT MSVC)
    target_compile_options(serenity_tests PUBLIC -march=${SERENITY_MARCH})
  endif()
  add_test(NAME serenity COMMAND serenity_tests)
endif()


####################
## Python Wrapper ##
####################

option(SERENITY_PYTHON_BINDINGS "Enable Python interface" OFF)
if (SERENITY_PYTHON_BINDINGS)
  message("#===============================#")
  message("| Python Interface Dependencies |")
  message("#===============================#")
  if(SERENITY_DOWNLOAD_DEPENDENCIES)
    include(ImportPybind11)
    import_pybind11()
  else()
    find_package(pybind11 REQUIRED)
  endif()
  # target
  #  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})
  pybind11_add_module(serenipy
    ${PROJECT_SOURCE_DIR}/src/python/serenipy.cpp
    ${SERENITY_PYTHON_FILES}
  )
  set_target_properties(serenipy
    PROPERTIES
      SUFFIX ".so"
      CXX_STANDARD 14
  )
  target_link_libraries(serenipy
    PUBLIC
      $<$<BOOL:${GCC_COVERAGE}>:-lgcov>
      $<$<BOOL:${GCC_PROFILE}>:-pg>
    PRIVATE
      serenity
  )
  set_target_properties(serenipy PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  )
  target_compile_options(serenipy
    PUBLIC
      $<$<BOOL:${WERROR}>:-Werror>
      $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
      $<$<BOOL:${GCC_PROFILE}>:-pg>
  )
  if(NOT "${SERENITY_MARCH}" STREQUAL "" AND NOT MSVC)
    target_compile_options(serenipy PUBLIC -march=${SERENITY_MARCH})
  endif()

  # Install
  install(TARGETS serenipy
    EXPORT serenityTargets
    LIBRARY DESTINATION lib
    )
endif(SERENITY_PYTHON_BINDINGS)


###############################
## Define 'make doc' command ##
###############################

message("#============================#")
message("| Documentation Dependencies |")
message("#============================#")
message("-----Doxygen-----")
find_package(Doxygen)
if(DOXYGEN_FOUND)
  add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${PROJECT_SOURCE_DIR}/doc/serenity.doxyfile
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc
  COMMENT "Generating API documentation with Doxygen" VERBATIM)
else (DOXYGEN_FOUND)
  add_custom_target(doc
  COMMENT "Doxygen needed in order to create documentation" VERBATIM)
endif(DOXYGEN_FOUND)


#################################
## Allow usage without install ##
#################################

option(SERENITY_USAGE_FROM_SOURCE "Generate files for easy usage without invoking 'make install'." ON)
if(SERENITY_USAGE_FROM_SOURCE)
  if (NOT EXISTS ${PROJECT_SOURCE_DIR}/serenity.sh)
    file(COPY ${PROJECT_SOURCE_DIR}/dev/templates/serenity.sh DESTINATION ${PROJECT_SOURCE_DIR})
  endif()
  # Create links to bin, lib and include
  add_custom_command(TARGET serenity POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_BINARY_DIR}/lib" "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  add_custom_command(TARGET serenity POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_BINARY_DIR}/bin" "${CMAKE_CURRENT_SOURCE_DIR}/bin")
endif(SERENITY_USAGE_FROM_SOURCE)


###########################
## CMake files for users ##
###########################

include(CMakePackageConfigHelpers)
# Config Version file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/serenity-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Config file
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/serenity-config.cmake"
  INSTALL_DESTINATION "lib/cmake/serenity"
)

# Install serenity-config.cmake and serenity-config-version.cmake
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/serenity-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/serenity-config-version.cmake"
  DESTINATION lib/cmake/serenity
)

# Add all targets to the build-tree export set
export(
  TARGETS
    xc
    xcfun
    libint2-static
    serenity
    serenity_exe
  FILE "${PROJECT_BINARY_DIR}/serenity-targets.cmake"
)

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE Serenity)

# Install the export set for use with the install-tree
install(EXPORT serenityTargets
  FILE serenity-targets.cmake
  DESTINATION lib/cmake/serenity
)
