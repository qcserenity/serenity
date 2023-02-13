function(import_laplace_minimax)
  # If the target already exists, do nothing
  if(TARGET laplace-minimax)
    return()
  endif()

  include(ExternalProject)

  set (LAPLACE_COMPILER "")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Using gfortran for laplace-minimax.")
    set(LAPLACE_COMPILER "gfortran")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Using ifort for laplace-minimax.")
    set(LAPLACE_COMPILER "ifort")
  else ()
    message(STATUS "Compiler family used is not support for laplace-minimax. Will try to compile and link with gfortran.")
    set(LAPLACE_COMPILER "gfortran")
  endif()

  ExternalProject_Add(laplace-minimax-static
    GIT_REPOSITORY https://github.com/bhelmichparis/laplace-minimax.git
    GIT_TAG 55414f3
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ext-laplace-minimax
    CONFIGURE_COMMAND cmake ../laplace-minimax-static -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_Fortran_COMPILER=${LAPLACE_COMPILER}
    BUILD_COMMAND make VERBOSE=1
    UPDATE_COMMAND ""
    INSTALL_COMMAND cp <BINARY_DIR>/liblaplace-minimax.a ${CMAKE_CURRENT_BINARY_DIR}/lib/liblaplace-minimax.a
  )

  add_library(laplace-minimax STATIC IMPORTED)
  set_property(TARGET laplace-minimax PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib/liblaplace-minimax.a)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set_property(TARGET laplace-minimax PROPERTY INTERFACE_LINK_LIBRARIES libgfortran.so)
    set(LAPLACE_COMPILER "gfortran")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set_property(TARGET laplace-minimax PROPERTY INTERFACE_LINK_LIBRARIES ifcore.so ifport.so)
  else ()
    set_property(TARGET laplace-minimax PROPERTY INTERFACE_LINK_LIBRARIES libgfortran.so)
  endif()
  set_property(TARGET laplace-minimax PROPERTY POSITION_INDEPENDENT_CODE ON)
  add_dependencies(laplace-minimax laplace-minimax-static)

  # Final check if all went well
  if(NOT TARGET laplace-minimax)
    string(CONCAT error_msg
      "Laplace-Minimax was not found and could not be established through a download."
    )
    message(FATAL_ERROR ${error_msg})
   endif()
endfunction()
