function(import_libecpint)
  # If the target already exists, do nothing
  if(TARGET Libecpint::Libecpint)
    return()
  endif()

  add_library(Libecpint::Libecpint SHARED IMPORTED)
  set_property(TARGET Libecpint::Libecpint PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/ext/libecpint/lib/libecpint.so)
  set_property(TARGET Libecpint::Libecpint PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/ext/libecpint/include)

  # Download it instead
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/ext/libecpint/lib/libecpint.so")
    file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/ext/libecpint/include)
    include(ExternalProject)
    ExternalProject_Add(Libecpint
    GIT_REPOSITORY "https://github.com/moritzBens/libecpint.git"
    CMAKE_ARGS -DENABLE_FORTRAN_INTERFACE=OFF -DENABLE_TESTALL=OFF -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/ext/libecpint)
    add_dependencies(Libecpint::Libecpint Libecpint)
  else()
    message(STATUS "Libecpint was found the 'ext/' folder.")
    return()    
  endif()

  # Final check if all went well
  if(TARGET Libecpint::Libecpint)
    message(STATUS
      "Libecpint was not found in your PATH, so it was downloaded."
    )
  else()
    string(CONCAT error_msg
      "Libecpint was not found in your PATH and could not be established "
      "through a download. Try specifying the Libecpint_DIR variable or " 
      "altering the CMAKE_PREFIX_PATH."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
