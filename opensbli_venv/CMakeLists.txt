# To specifiy the cmake version range Currenty work under 3.18
cmake_minimum_required(VERSION 3.22)
project(OSBLI_VENV CXX)

if(DEFINED ENV{HDF5_INSTALL_PATH})
  message(STATUS "HDF5 location $ENV{HDF5_INSTALL_PATH}")
  set(HDF5_INSTALL_PATH $ENV{HDF5_INSTALL_PATH} CACHE PATH "Env path to HDF5")
  option(OSBLI_BUILD_HDF5 "Build the HDF5 library with the current loaded compiler" OFF)
else ()
  option(OSBLI_BUILD_HDF5 "Build the HDF5 library with the current loaded compiler" ON)
endif()

if(DEFINED ENV{OPS_INSTALL_DIR})
  message(STATUS "OPS location $ENV{OPS_INSTALL_DIR}")
  set(OPS_INSTALL_DIR $ENV{OPS_INSTALL_DIR} CACHE PATH "Env path to OPS")
  option(OSBLI_BUILD_OPS "Build the OPS library with the current loaded compiler" OFF)
else ()
  option(OSBLI_BUILD_OPS "Build the OPS library with the current loaded compiler" ON)
endif()

# Check is HDF5 and OPS needs to be buind and installed
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/HDF5" "${CMAKE_SOURCE_DIR}/cmake/OPS") 

include(OSBLI_FindPackages)


