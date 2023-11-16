function(import_libint)
  # If the target already exists, do nothing
  if(TARGET libint2-static)
    return()
  endif()

  message(STATUS
    "Checking Libint2 source files."
  )
  include(DownloadProject)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ext-libint-src)
  download_project(
    PROJ ext-libint
    #GIT_REPOSITORY https://thclab.uni-muenster.de/serenity/libint.git
    #GIT_SHALLOW TRUE
    #GIT_PROGRESS TRUE
    #SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/ext-libint-git
    #PATCH_COMMAND tar -xzf libint-2.7.0-beta.6.tgz --directory=${CMAKE_CURRENT_BINARY_DIR}/ext-libint-src --strip-component=1
    URL https://www.uni-muenster.de/Chemie.oc/THCLAB/libint/libint-2.7.0-beta.6.tgz
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/ext-libint-src
    QUIET
  )
  set(LIBINT2_BUILD_SHARED_AND_STATIC_LIBS ON CACHE BOOL "" FORCE)
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/ext-libint-src ${CMAKE_CURRENT_BINARY_DIR}/ext-libint-build)

  # Final check if all went well
  if(NOT TARGET libint2-static)
    string(CONCAT error_msg
      "Libint was not found and could not be established through a download."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
