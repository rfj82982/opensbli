if (OSBLI_BUILD_OPS)
  message(STATUS "Downloading and Building OPS")
  configure_file(${CMAKE_SOURCE_DIR}/cmake/OPS/downloadinstallOPS.cmake.in ops-build/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
           RESULT_VARIABLE result
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ops-build )
  if(result)
	  message(FATAL_ERROR "CMake step for OPS failed: ${result}")
  else()
	  message(STATUS "CMake step for OPS completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ops-build )
  if(result)
	  message(FATAL_ERROR "Build step for OPS failed: ${result}")
  endif()
  set(OPS_INSTALL_DIR ${CMAKE_SOURCE_DIR}/osbli_opt CACHE PATH "Env path to OPS")
endif()
#set(OPS_TRANSLATOR ${OPS_INSTALL_DIR}/translator/ops_translator/ops-translator CACHE PATH "Env path to OPS Translator")
set(OPS_TRANSLATOR ${OPS_INSTALL_DIR}/ops_translator/ops-translator CACHE PATH "Env path to OPS Translator")
# This is set only in the process and won't be available in the wider enviroment
set(ENV{OPS_INSTALL_DIR} ${OPS_INSTALL_DIR})
set(ENV{OPS_TRANSLATOR} ${OPS_TRANSLATOR})
message(STATUS "OPS location $ENV{OPS_INSTALL_DIR}")
message(STATUS "OPS translator $ENV{OPS_TRANSLATOR}")

