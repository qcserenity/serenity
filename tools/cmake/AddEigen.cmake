# If the target already exists, do nothing
if(NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 3.3.3 REQUIRED)
endif()

option(SERENITY_USE_INTEL_MKL "Use the Intel MKL libraries with Eigen" ON)
option(SERENITY_USE_LAPACK "Use a LAPACK library with Eigen" ON)
option(SERENITY_USE_BLAS "Use a BLAS library with Eigen" ON)

# Attempt to find external linalg libraries that accelerate basic calls in
# the following order:
#
# 1. Intel MKL
# 2. LAPACK (brings BLAS too)
# 3. BLAS
#
if(NOT ADD_EIGEN_SEARCHED_EXTERNAL_LINALG_LIBRARIES)
  if(SERENITY_USE_INTEL_MKL)
    include(FindMKL)
  endif()

  if(MKL_FOUND)
    find_package(OpenMP REQUIRED)
    message(STATUS "Found MKL for use with Eigen3")
  else()
    if(SERENITY_USE_LAPACK)
      include(FindLAPACK)
      find_package(LAPACK QUIET)
    endif()

    if(LAPACK_FOUND)
      message(STATUS "Found LAPACK/BLAS for use with Eigen3")
    else()
      if(SERENITY_USE_BLAS)
        include(FindBLAS)
        find_package(BLAS QUIET)
      endif()

      if(BLAS_FOUND)
        message(STATUS "Found BLAS for use with Eigen3")
      endif()
    endif()
  endif()

  set(ADD_EIGEN_SEARCHED_EXTERNAL_LINALG_LIBRARIES TRUE)
endif()

function(add_eigen target_name mode)
  # Everything relating to the special case of object libraries can be
  # simplified down to the target_link_libraries command when the minimum cmake
  # version is bumped past 3.12.0
  if("${CMAKE_VERSION}" VERSION_GREATER_EQUAL "3.12.0")
    set(_CAN_LINK_OBJECT_LIBRARIES TRUE)
  endif()
  get_target_property(target_type ${target_name} TYPE)
  if("${target_type}" STREQUAL "OBJECT_LIBRARY")
    set(_IS_OBJECT_LIBRARY TRUE)
  endif()

  # Append the required properties to the passed target
  if(MKL_FOUND AND SERENITY_USE_INTEL_MKL)
    # MKL AND EIGEN_USE_MKL_ALL
    if(_CAN_LINK_OBJECT_LIBRARIES OR NOT _IS_OBJECT_LIBRARY)
      target_link_libraries(${target_name} ${mode} Eigen3::Eigen ${MKL_LIBRARIES} OpenMP::OpenMP_CXX)
    else()
      target_include_directories(${target_name} SYSTEM ${mode} $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>)
    endif()
    target_include_directories(${target_name} SYSTEM ${mode} ${MKL_INCLUDE_DIRS})
    target_compile_definitions(${target_name} ${mode} EIGEN_USE_MKL_ALL)
    target_compile_options(${target_name} ${mode} $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>)
  else()
    if(LAPACK_FOUND AND SERENITY_USE_LAPACK)
      # LAPACK and EIGEN_USE_LAPACK / EIGEN_USE_BLAS
      if(_CAN_LINK_OBJECT_LIBRARIES OR NOT _IS_OBJECT_LIBRARY)
        target_link_libraries(${target_name} ${mode} Eigen3::Eigen ${LAPACK_LIBRARIES})
      endif()
      target_compile_definitions(${target_name} ${mode} EIGEN_USE_LAPACK EIGEN_USE_BLAS)
    else()
      if(BLAS_FOUND AND SERENITY_USE_BLAS)
        # Blas and EIGEN_USE_BLAS
        if(_CAN_LINK_OBJECT_LIBRARIES OR NOT _IS_OBJECT_LIBRARY)
          target_link_libraries(${target_name} ${mode} Eigen3::Eigen ${BLAS_LIBRARIES})
        else()
          target_include_directories(${target_name} SYSTEM ${mode} $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>)
        endif()
        target_compile_definitions(${target_name} ${mode} EIGEN_USE_BLAS)
      else()
        # Just Eigen, no definitions
        if(_CAN_LINK_OBJECT_LIBRARIES OR NOT _IS_OBJECT_LIBRARY)
          target_link_libraries(${target_name} ${mode} Eigen3::Eigen)
        else()
          target_include_directories(${target_name} SYSTEM ${mode} $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>)
        endif()
      endif()
    endif()
  endif()
endfunction()
