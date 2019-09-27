# - Find Intel MKL
# Try to find the MKL libraries as needed for Serenity
#

include(FindPackageHandleStandardArgs)

set(INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains intel libs")
set(MKL_ROOT $ENV{MKLROOT} CACHE PATH "Folder contains MKL")

# Find include dir
find_path(MKL_INCLUDE_DIRS mkl.h
          PATHS ${MKL_ROOT}/include)

# Set path according to architecture
execute_process(COMMAND getconf LONG_BIT OUTPUT_VARIABLE SYSTEM_BIT)
string(STRIP ${SYSTEM_BIT} SYSTEM_BIT)
if ("${SYSTEM_BIT}" STREQUAL "64")
  set(MKL_LIB_ARCH lib/intel64/)
else()
  set(MKL_LIB_ARCH lib/ia32/)
endif()

# Find libraries
if ("${SYSTEM_BIT}" STREQUAL "64")
  find_library(MKL_INTERFACE_LIBRARY libmkl_intel_lp64.so  PATHS ${MKL_ROOT}/lib/intel64/)
else()
  find_library(MKL_INTERFACE_LIBRARY libmkl_intel.so  PATHS ${MKL_ROOT}/lib/ia32/)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
  find_library(MKL_THREADING_LIBRARY libmkl_intel_thread.so PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  find_library(MKL_THREADING_LIBRARY libmkl_gnu_thread.so PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
else()
  unset(MKL_FOUND)
endif()

find_library(MKL_CORE_LIBRARY libmkl_core.so PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
find_library(MKL_AVX_LIBRARY libmkl_avx2.so PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
find_library(MKL_DEF_LIBRARY libmkl_def.so PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})

set(MKL_LIBRARIES ${MKL_AVX_LIBRARY} ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY})

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIRS MKL_LIBRARIES)

if(MKL_FOUND)
  message("-- Will use Eigen3 with interface to SMP parallel MKL libraries, using:")
  message("-- ${MKL_INCLUDE_DIRS} ")
  message("-- ${MKL_INTERFACE_LIBRARY} ")
  message("-- ${MKL_THREADING_LIBRARY} ")
  message("-- ${MKL_CORE_LIBRARY} ")
  message("-- ${MKL_AVX_LIBRARY} ")
  message("-- ${MKL_DEF_LIBRARY} ")
else()
  set(MKL_LIBRARIES "MKL_LIBRARIES-NOTFOUND")
endif()
