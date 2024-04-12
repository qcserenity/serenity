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
  find_library(MKL_INTERFACE_LIBRARY NAMES libmkl_intel_lp64.so libmkl_intel_lp64.so.1 libmkl_intel_lp64.so.2 PATHS ${MKL_ROOT}/lib/intel64/)
else()
  find_library(MKL_INTERFACE_LIBRARY NAMES libmkl_intel.so libmkl_intel.so.1 libmkl_intel.so.2 PATHS ${MKL_ROOT}/lib/ia32/)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
  find_library(MKL_THREADING_LIBRARY NAMES libmkl_intel_thread.so libmkl_intel_thread.so.1 libmkl_intel_thread.so.2 PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  find_library(MKL_THREADING_LIBRARY NAMES libmkl_gnu_thread.so libmkl_gnu_thread.so.1 libmkl_gnu_thread.so.2 PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  find_library(MKL_THREADING_LIBRARY NAMES libmkl_gnu_thread.so libmkl_gnu_thread.so.1 libmkl_gnu_thread.so.2 PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
else()
  unset(MKL_FOUND)
endif()

find_library(MKL_CORE_LIBRARY NAMES libmkl_core.so libmkl_core.so.1 libmkl_core.so.2 PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
find_library(MKL_AVX2_LIBRARY NAMES libmkl_avx2.so libmkl_avx2.so.1 libmkl_avx2.so.2 PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})
find_library(MKL_VML_AVX2_LIBRARY NAMES libmkl_vml_avx2.so libmkl_vml_avx2.so.1 libmkl_vml_avx2.so.2 PATHS ${MKL_ROOT}/${MKL_LIB_ARCH})

set(MKL_LIBRARIES ${MKL_AVX2_LIBRARY} ${MKL_VML_AVX2_LIBRARY} ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY})

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIRS MKL_LIBRARIES)

if(MKL_FOUND)
  message("-- Will use Eigen3 with interface to SMP parallel MKL libraries, using:")
  message("-- ${MKL_INCLUDE_DIRS} ")
  message("-- ${MKL_AVX2_LIBRARY} ")
  message("-- ${MKL_VML_AVX2_LIBRARY} ")
  message("-- ${MKL_CORE_LIBRARY} ")
  message("-- ${MKL_INTERFACE_LIBRARY} ")
  message("-- ${MKL_THREADING_LIBRARY} ")
else()
  set(MKL_LIBRARIES "MKL_LIBRARIES-NOTFOUND")
endif()
