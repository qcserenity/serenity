function(import_xcfun)
  # If the target already exists, do nothing
  if(TARGET xcfun)
    return()
  endif()

  message(STATUS
    "Checking XCFun source files."
  )
  include(DownloadProject)
  download_project(
    PROJ ext-xcfun
    GIT_REPOSITORY https://github.com/qcserenity/xcfun.git
    QUIET
  )

  set(_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(_ENABLE_TESTALL ${ENABLE_TESTALL})
  set(ENABLE_TESTALL OFF CACHE BOOL "XCFun: Enable all tests")
  set(BUILD_SHARED_LIBS OFF)
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/ext-xcfun-src ${CMAKE_CURRENT_BINARY_DIR}/ext-xcfun-build)
  set(XCFUN_MAX_ORDER 3 CACHE STRING "XCFun: Maximum order of derivatives of the exchange-correlation kernel" FORCE)
  mark_as_advanced(XCFUN_PYTHON_INTERFACE XCFUN_MAX_ORDER PYMOD_INSTALL_LIBDIR)
  set(ENABLE_TESTALL ${_ENABLE_TESTALL})
  set(BUILD_SHARED_LIBS ${_BUILD_SHARED_LIBS})
  install(TARGETS xcfun EXPORT serenityTargets DESTINATION lib)

  # Final check if all went well
  if(NOT TARGET xcfun)
    string(CONCAT error_msg
      "XCFun was not found and could not be established through a download."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
