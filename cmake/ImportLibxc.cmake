function(import_libxc)
  # If the target already exists, do nothing
  if(TARGET xc)
    return()
  endif()

  message(STATUS
    "Checking LibXC source files."
  )
  include(DownloadProject)
  download_project(
    PROJ ext-libxc
    GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
    GIT_TAG 6.1.0
    QUIET
  )

  set(_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(_BUILD_TESTING ${BUILD_TESTING})
  set(BUILD_SHARED_LIBS OFF)
  set(BUILD_TESTING OFF)
  set(DISABLE_KXC OFF CACHE BOOL "Libxc: Don't compile third derivative code")
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/ext-libxc-src ${CMAKE_CURRENT_BINARY_DIR}/ext-libxc-build)
  set(BUILD_TESTING ${_BUILD_TESTING})
  set(BUILD_SHARED_LIBS ${_BUILD_SHARED_LIBS})
  install(TARGETS xc EXPORT serenityTargets DESTINATION lib)
  target_include_directories(xc
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ext-libxc-build>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ext-libxc-build/src>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ext-libxc-src/src>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ext-libxc-build/gen_funcidx>
  )

  # Final check if all went well
  if(NOT TARGET xc)
    string(CONCAT error_msg
      "LibXC was not found and could not be established through a download."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
