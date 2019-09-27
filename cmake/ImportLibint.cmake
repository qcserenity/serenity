function(import_libint)
  # If the target already exists, do nothing
  if(TARGET Libint::Libint)
    return()
  endif()

  if (APPLE)
    set(LIBINT_LIBRARY_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}int2-beta.3.2${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(LIBINT_LIBRARY_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}int2-beta.3${CMAKE_SHARED_LIBRARY_SUFFIX}.2)
  endif()

  add_library(Libint::Libint SHARED IMPORTED)
  set_property(TARGET Libint::Libint PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/ext/libint/${ARCH_LIB_PATH}/${LIBINT_LIBRARY_NAME})
  set_property(TARGET Libint::Libint PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                ${PROJECT_SOURCE_DIR}/ext/libint/include
                ${PROJECT_SOURCE_DIR}/ext/libint/include/libint2
  )
  install(FILES
    ${PROJECT_SOURCE_DIR}/ext/libint/${ARCH_LIB_PATH}/${LIBINT_LIBRARY_NAME}
    DESTINATION lib
  )

  # Download it instead
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/ext/libint/${ARCH_LIB_PATH}/${LIBINT_LIBRARY_NAME}")
    option(SERENITY_DOWNLOAD_LIBINT "Download a compiled version of libint rather than compiling, if possible." ON)
    set(LIBINT_DOWNLOAD_UNSUPPORTED OFF)
    # Download binary if requested and possible
    if (SERENITY_DOWNLOAD_LIBINT)
      file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/ext)
      # Check if the binary version is supported
      if (NOT (CMAKE_SIZEOF_VOID_P EQUAL 8))
	      # x32 unsupported
        set(LIBINT_DOWNLOAD_UNSUPPORTED ON)
      endif()
      # Download
      if (NOT LIBINT_DOWNLOAD_UNSUPPORTED)
        include(DownloadProject)
        download_project(
          PROJ libint
          GIT_REPOSITORY https://thclab.uni-muenster.de/serenity/ext-libint.git
          GIT_TAG linux-x64
          GIT_SHALLOW TRUE
          GIT_PROGRESS TRUE
          PATCH_COMMAND tar -xf libint.tgz --directory ${PROJECT_SOURCE_DIR}/ext/
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ""
          QUIET
        )
        # Ubuntu has only a lib/ folder even if x64 this patches the path in the downloaded archive
        if(NOT (${ARCH_LIB_PATH} STREQUAL "lib64"))
          file(RENAME ${PROJECT_SOURCE_DIR}/ext/libint/lib64 ${PROJECT_SOURCE_DIR}/ext/libint/${ARCH_LIB_PATH})
        endif()
      endif()
    endif()
    # Compile fresh libint
    if (NOT SERENITY_DOWNLOAD_LIBINT OR LIBINT_DOWNLOAD_UNSUPPORTED)
      include(ExternalProject)
      ExternalProject_Add(libint
        GIT_REPOSITORY https://thclab.uni-muenster.de/serenity/ext-libint.git
        GIT_TAG master
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        PATCH_COMMAND  tar -xvf libint-2.3.0-beta.3.tgz --directory <PROJECT_SOURCE_DIR/src> --strip-components=1
        CONFIGURE_COMMAND ./configure "CXXFLAGS=${CMAKE_CXX_FLAGS}" --enable-shared  --prefix=${PROJECT_SOURCE_DIR}/ext/libint  --libdir=${PROJECT_SOURCE_DIR}/ext/libint/${ARCH_LIB_PATH}
        BUILD_IN_SOURCE 1
        BUILD_COMMAND $(MAKE)
        INSTALL_COMMAND make install
      )
    endif()
    file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/ext/libint/include/libint2)
    add_dependencies(Libint::Libint libint)
  else()
    message(STATUS "Libint was found the 'ext/' folder.")
    return()
  endif()

  # Final check if all went well
  if(TARGET Libint::Libint)
    message(STATUS
      "Libint was not found in your PATH, so it was downloaded."
    )
  else()
    string(CONCAT error_msg
      "Libint was not found in your PATH and could not be established "
      "through a download. Try specifying the Libint_DIR variable or "
      "altering the CMAKE_PREFIX_PATH."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
