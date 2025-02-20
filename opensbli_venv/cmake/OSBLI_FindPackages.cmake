# MPI
find_package(MPI QUIET)
if(NOT MPI_FOUND)
  message(WARNING "MPI environment NOT FOUND! Only sequential code will be compiled!")
endif()
message(STATUS "MPICXX FOUND ${MPI_CXX_FOUND}")

# HDF5
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/HDF5")
include(OSBLI_HDF5)

# OPS
include(OSBLI_OPS)

# OSBLI
set(OSBLI_PYTHON_PATH ${CMAKE_SOURCE_DIR}/../ CACHE PATH "Env path to OSBLI")
# This is set only in the process and won't be available in the wider enviroment
set(ENV{OSBLI_PYTHON_PATH} ${OSBLI_PYTHON_PATH})

# Write of the configuration file
file(WRITE ${CMAKE_SOURCE_DIR}/osbli_env.sh
"#!/bin/bash
export PYTHONPATH=$PYTHONPATH:${OSBLI_PYTHON_PATH}
export HDF5_INSTALL_PATH=${HDF5_INSTALL_PATH}
export OPS_INSTALL_DIR=${OPS_INSTALL_DIR} 
export OPS_TRANSLATOR=${OPS_TRANSLATOR}
")

