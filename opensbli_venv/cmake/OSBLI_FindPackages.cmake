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
set(osbli_venv_dir "$ENV{osbli_venv_dir}")
set(osbli_venv_activate "$ENV{osbli_venv_activate}")

# Write of the configuration file
file(WRITE ${CMAKE_SOURCE_DIR}/osbli_env.sh
"#!/bin/bash
export PYTHONPATH=$PYTHONPATH:${OSBLI_PYTHON_PATH}
export HDF5_INSTALL_PATH=${HDF5_INSTALL_PATH}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF5_INSTALL_PATH/lib
export OPS_INSTALL_DIR=${OPS_INSTALL_DIR} 
export OPS_TRANSLATOR=${OPS_TRANSLATOR}
alias ops_translation=\"source ${OPS_INSTALL_DIR}/translator/ops_translator/ops_venv/bin/activate && python ${OPS_TRANSLATOR}/ops.py opensbli.cpp\"
#
osbli_code_generation (){
  if [ -z \"$1\" ];then
    echo \"One input file needs to be provided \"
    return
  elif [ -z \"$2\" ];then
    input_file=$1
  else
    echo \"No single input file given \"
    return
  fi
  source ${osbli_venv_activate} && python $input_file  
}

#
osbli_cmake_command (){
  cp ${osbli_venv_dir}/../../apps/CMakeLists.txt . 
  source ${OPS_INSTALL_DIR}/translator/ops_translator/ops_venv/bin/activate && cmake -S . -B test-build -DOPS_INSTALL_DIR=${OPS_INSTALL_DIR} -DCMAKE_BUILD_TYPE=Release -DLEGACY_CODEGEN=OFF -DHDF5_ROOT=$HDF5_INSTALL_PATH  
}
")

