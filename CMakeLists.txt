cmake_minimum_required(VERSION 3.9)

project(serenity
        VERSION 1.1.0
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
option(NATIVE "Compile with -march=native" ON)
option(WERROR "Compile with warnings as errors" OFF)
option(GCC_PROFILE "Compile with profile flags" OFF)
option(GCC_COVERAGE "Compile with coverage flags" OFF)
if (GCC_COVERAGE)
  set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
endif()
option(SERENITY_PREFER_XCFUN "If both XCFun and LibXC are present, prefer XCFun per default" ON)
option(SERENITY_USE_XCFUN "Include XCFun as DFT functional library" ON)
option(SERENITY_USE_LIBXC "Include LibXC as DFT functional library" ON)
if ((NOT SERENITY_USE_LIBXC) AND (NOT SERENITY_USE_XCFUN))
  message(FATAL_ERROR "At least one of SERENITY_USE_LIBXC and SERENITY_USE_XCFUN has to be true.")
endif()

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

add_library(${PROJECT_NAME} ${SERENITY_CPPS} ${SERENITY_HEADERS})
set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
add_executable(${PROJECT_NAME}_exe ${PROJECT_SOURCE_DIR}/src/serenity.cpp)

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
add_eigen(${PROJECT_NAME} PUBLIC)

# OpenMP
message("-----Checking OpenMP-----")
find_package(OpenMP REQUIRED)

# Libint2
message("-----Checking Libint2-----")
include(ImportLibint)
import_libint()

if (SERENITY_USE_XCFUN)
  # xcfun
  message("-----Checking XCFun----")
  include(ImportXCFun)
  import_xcfun()
endif()

if (SERENITY_USE_LIBXC)
  # libxc
  message("-----Checking Libxc----")
  include(ImportLibxc)
  import_libxc()
endif()

# HDF5
message("-----Checking HDF5-----")
if (NOT DEFINED HDF5_USE_STATIC_LIBRARIES)
  set(HDF5_USE_STATIC_LIBRARIES OFF)
endif()
find_package(HDF5 REQUIRED COMPONENTS HL CXX)
if (HDF5_USE_STATIC_LIBRARIES)
  # The following if is a small helper for systems
  #  that have HDF5 installed but do not have the
  #  cmake files installed (this seems common for
  #  many distributions)
  if(NOT (TARGET hdf5::hdf5_cpp-static))
    add_library(hdf5::hdf5_cpp-static INTERFACE IMPORTED)
    set_property(TARGET hdf5::hdf5_cpp-static PROPERTY
      INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES}"
    )
  endif()
  set(HDF5_TARGET hdf5::hdf5_cpp-static)
else()
  # The following if is a small helper for systems
  #  that have HDF5 installed but do not have the
  #  cmake files installed (this seems common for
  #  many distributions)
  if(NOT (TARGET hdf5::hdf5_cpp-shared))
    add_library(hdf5::hdf5_cpp-shared INTERFACE IMPORTED)
    set_property(TARGET hdf5::hdf5_cpp-shared PROPERTY
      INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES}"
    )
  endif()
  set(HDF5_TARGET hdf5::hdf5_cpp-shared)
endif()

# Boost
message("-----Checking Boost-----")
find_package(Boost REQUIRED)

# Libecpint
message("-----Checking Libecpint-----")
include(ImportLibecpint)
import_libecpint()

#################################
## Add dependencies to targets ##
#################################

# Shared library
set_target_properties(${PROJECT_NAME}
  PROPERTIES PUBLIC_HEADER
    "${SERENITY_HEADERS}"
)
target_include_directories(${PROJECT_NAME}
  PUBLIC
    #$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/${PROJECT_NAME}>
    ${HDF5_INCLUDE_DIRS}
    $<$<BOOL:${MKL_FOUND}>:${MKL_INCLUDE_DIRS}>
)
target_compile_options(${PROJECT_NAME} PRIVATE
  ${OpenMP_CXX_FLAGS}
)
set_target_properties(${PROJECT_NAME} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  CXX_STANDARD 14
)
target_link_libraries(${PROJECT_NAME}
  PUBLIC
    Boost::boost
    OpenMP::OpenMP_CXX
    $<$<BOOL:${GCC_COVERAGE}>:-lgcov>
    $<$<BOOL:${GCC_PROFILE}>:-pg>
  PRIVATE
    libint2-static
    xc
    xcfun
    ecpint
    ${HDF5_TARGET}
)
target_compile_definitions(${PROJECT_NAME}
  PUBLIC
    $<$<BOOL:${SERENITY_PREFER_XCFUN}>:SERENITY_PREFER_XCFUN>
    $<$<BOOL:${SERENITY_USE_XCFUN}>:SERENITY_USE_XCFUN>
    $<$<BOOL:${SERENITY_USE_LIBXC}>:SERENITY_USE_LIBXC>
    $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
    $<$<BOOL:${MKL_FOUND}>:EIGEN_USE_MKL_ALL>
)
target_compile_options(${PROJECT_NAME}
  PUBLIC
    -Wall
    -Wextra
    -Wno-comment
    $<$<BOOL:${NATIVE}>:-march=native>
    $<$<BOOL:${WERROR}>:-Werror>
    $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
    $<$<BOOL:${GCC_PROFILE}>:-pg>
    $<$<AND:$<BOOL:${MKL_FOUND}>,$<STREQUAL:"${CMAKE_CXX_COMPILER_ID}","Intel">>:-DMKL_LP64>
)

# Executable
set_target_properties(${PROJECT_NAME}_exe PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_exe
  PUBLIC
    $<$<BOOL:${GCC_COVERAGE}>:-lgcov>
    $<$<BOOL:${GCC_PROFILE}>:-pg>
  PRIVATE
    ${PROJECT_NAME}
    libint2-static
)
set_target_properties(${PROJECT_NAME}_exe PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  CXX_STANDARD 14
)
target_compile_options(${PROJECT_NAME}_exe
  PUBLIC
    -Wall
    -Wextra
    -Wno-comment
    $<$<BOOL:${NATIVE}>:-march=native>
    $<$<BOOL:${WERROR}>:-Werror>
    $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
    $<$<BOOL:${GCC_PROFILE}>:-pg>
)

#############
## Install ##
#############

# Targets
install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_exe
        EXPORT ${PROJECT_NAME}Targets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include/${PROJECT_NAME}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data/basis
        ${CMAKE_CURRENT_SOURCE_DIR}/data/initialGuess
        DESTINATION share/${PROJECT_NAME}/data
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
  include(ImportGTest)
  import_gtest()
  message("-----Adding Test Executable-----")
  # Test exe
  add_executable(${PROJECT_NAME}_tests ${SERENITY_TEST_FILES})
  target_link_libraries(${PROJECT_NAME}_tests
    PUBLIC
      $<$<BOOL:${GCC_COVERAGE}>:-lgcov>
      $<$<BOOL:${GCC_PROFILE}>:-pg>
    PRIVATE
      GTest::Main
      GMock::GMock
      ${PROJECT_NAME}
      libint2-static
      xcfun
      ecpint
      ${HDF5_TARGET}
  )
  set_target_properties(${PROJECT_NAME}_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    CXX_STANDARD 14
  )
  target_compile_options(${PROJECT_NAME}_tests
    PUBLIC
      -Wall
      -Wextra
      -Wno-comment
      $<$<BOOL:${NATIVE}>:-march=native>
      $<$<BOOL:${WERROR}>:-Werror>
      $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
      $<$<BOOL:${GCC_PROFILE}>:-pg>
  )
  add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME}_tests)
endif()


####################
## Python Wrapper ##
####################

option(SERENITY_PYTHON_BINDINGS "Enable Python interface" OFF)
if (SERENITY_PYTHON_BINDINGS)
  message("#===============================#")
  message("| Python Interface Dependencies |")
  message("#===============================#")
  include(ImportPybind11)
  import_pybind11()

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
      ${PROJECT_NAME}
      libint2-static
  )
  set_target_properties(serenipy PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  )
  target_compile_options(serenipy
    PUBLIC
      $<$<BOOL:${NATIVE}>:-march=native>
      $<$<BOOL:${WERROR}>:-Werror>
      $<$<BOOL:${GCC_COVERAGE}>:-fprofile-arcs -ftest-coverage>
      $<$<BOOL:${GCC_PROFILE}>:-pg>
  )

  # Install
  install(TARGETS serenipy
    EXPORT ${PROJECT_NAME}Targets
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
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_BINARY_DIR}/lib" "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_BINARY_DIR}/bin" "${CMAKE_CURRENT_SOURCE_DIR}/bin")
endif(SERENITY_USAGE_FROM_SOURCE)


###########################
## CMake files for users ##
###########################

include(CMakePackageConfigHelpers)
# Config Version file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Config file
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
)

# Install serenity-config.cmake and serenity-config-version.cmake
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION lib/cmake/${PROJECT_NAME}
)

# Add all targets to the build-tree export set
export(
  TARGETS
    ecpint
    xc
    xcfun
    libint2-static
    ${PROJECT_NAME}
    ${PROJECT_NAME}_exe
  FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake"
)

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE Serenity)

# Install the export set for use with the install-tree
install(EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION lib/cmake/${PROJECT_NAME}
)

