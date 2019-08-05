function(import_pybind11)
  # If the target already exists, do nothing
  if(TARGET pybind11)
      return()
  endif()

  # Download it instead
  include(DownloadProject)
  download_project(
    PROJ pybind11
    GIT_REPOSITORY      https://github.com/pybind/pybind11.git
    GIT_TAG             v2.2.4
    QUIET
  )

  add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})

  # Final check if all went well
  if(EXISTS "${pybind11_SOURCE_DIR}/CMakeLists.txt")
    message(STATUS
      "Pybind11 was not found in your PATH, so it was downloaded."
    )
  else()
    string(CONCAT error_msg
      "Pybind11 was not be established through a download."
    )
    message(FATAL_ERROR ${error_msg})
  endif()
endfunction()
