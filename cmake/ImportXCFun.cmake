function(import_xcfun)
  # If the target already exists, do nothing
  if(TARGET XCFun::XCFun)
    return()
  endif()

  add_library(XCFun::XCFun SHARED IMPORTED)
  set_property(TARGET XCFun::XCFun PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/ext/xcfun/lib/libxcfun.so)
  set_property(TARGET XCFun::XCFun PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/ext/xcfun/include)

  # Download it instead
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/ext/xcfun/lib/libxcfun.so")
    file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/ext/xcfun/include/)
    include(ExternalProject)
    ExternalProject_Add(XCFun
     GIT_REPOSITORY "https://github.com/moritzBens/xcfun.git"
     GIT_TAG "b2bfa8e92fde74897f08a5f412196df2677af956"
     CMAKE_ARGS -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=ON -DENABLE_FORTRAN_INTERFACE=OFF -DENABLE_TESTALL=OFF -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/ext/xcfun)
    add_dependencies(XCFun::XCFun XCFun)
  else()
    message(STATUS "XCFun was found the 'ext/' folder.")
    return()
  endif()

  # Final check if all went well
  if(TARGET XCFun::XCFun)
    message(STATUS
      "XCFun was not found in your PATH, so it was downloaded."
    )
  else()
    string(CONCAT error_msg
      "XCFun was not found in your PATH and could not be established "
      "through a download. Try specifying the XCFun_DIR variable or " 
      "altering the CMAKE_PREFIX_PATH"
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
