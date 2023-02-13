function(import_libecpint)
  # If the target already exists, do nothing
  if(TARGET ecpint)
    return()
  endif()

  message(STATUS
    "Checking Libecpint source files."
  )
  include(DownloadProject)
  download_project(
    PROJ ext-ecpint
    GIT_REPOSITORY "https://github.com/robashaw/libecpint.git"
    GIT_TAG v1.0.7 
    QUIET
  )
  set(LIBECPINT_BUILD_TESTS OFF CACHE BOOL "Disable tests")
  set(LIBECPINT_BUILD_DOCS OFF CACHE BOOL "Disable docs")
  set(_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(BUILD_SHARED_LIBS OFF)
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/ext-ecpint-src ${CMAKE_CURRENT_BINARY_DIR}/ext-ecpint-build)
  set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})


  # Final check if all went well
  if(NOT TARGET ecpint)
    string(CONCAT error_msg
      "Libecpint was not found and could not be established through a download."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
