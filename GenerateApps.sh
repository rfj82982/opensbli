#!/bin/bash
##
###############
if [ $# = 2 ]
then
 path_hdf5=$1
 path_ops=$2
else
  echo "This Scripts is building the Apps for the OpenSBLI framework"
  echo "The Build is performed by executing the python script apps/generate_all_applications.py"
  echo "A local Python virtual enviroment is created to have all correct dependancies and particularly Python3.7 and SymPy==1.1"
  echo "This script mush be excuted from the root directory of the OpenSBLI repository"
  echo "Two inputs are necessary"
  echo "  1=> Path to HDF5 installation (i.e. HDF5 bin location)"
  echo "  2=> Path to OPS installation  (i.e. OPS bin location)"
  exit
fi
path_sbli=`pwd`
export PYTHONPATH=$PYTHONPATH:${path_sbli}
export HDF5_INSTALL_PATH=${path_hdf5}
export OPS_INSTALL_DIR=${path_ops} 
export OPS_TRANSLATOR=$OPS_INSTALL_DIR/translator/ops_translator/ops-translator

echo "HDF5 $HDF5_INSTALL_PATH"
echo "OPS  $OPS_INSTALL_DIR"

cd opensbli_venv/
source setup_venv.sh 
echo "My activate command ${python_activate}"
cd ../apps
python generate_all_applications.py
cd ..

