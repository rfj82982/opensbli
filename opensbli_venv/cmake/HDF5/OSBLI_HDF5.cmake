if (OSBLI_BUILD_HDF5)
  message(STATUS "Downloading and Building HDF5")
  enable_language(C)
  configure_file(${CMAKE_SOURCE_DIR}/cmake/HDF5/downloadinstallHDF5.cmake.in hdf5-build/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
           RESULT_VARIABLE result
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5-build )
  if(result)
	  message(FATAL_ERROR "CMake step for HDF5 failed: ${result}")
  else()
	  message(STATUS "CMake step for OPS completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5-build )
  if(result)
	  message(FATAL_ERROR "Build step for HDF5 failed: ${result}")
  endif()
  set(HDF5_INSTALL_PATH ${CMAKE_SOURCE_DIR}/osbli_opt CACHE PATH "Env path to HDF5")
  # This is set only in the process and won't be available in the wider enviroment
  set(ENV{HDF5_INSTALL_PATH} ${HDF5_INSTALL_PATH})
  message(STATUS "HDF5 location $ENV{HDF5_INSTALL_PATH}")
endif()

#find_package(HDF5 QUIET COMPONENTS C HL PATHS ${HDF5_INSTALL_PATH})
set(ENV{HDF5_ROOT} ${HDF5_INSTALL_PATH})
find_package(HDF5 QUIET COMPONENTS C HL)
if(NOT HDF5_FOUND)
  message(WARNING "HDF5 support NOT FOUND. CHECK THE PATH ${HDF5_INSTALL_PATH}")
else()
  set(HDF5_PREFER_PARALLEL true)
endif()
