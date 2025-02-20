#!/bin/bash
##
###############
echo "This Scripts is setting up the enviroment to build/test OpenSBLI Apps"
echo "A local Python virtual enviroment is created to have all correct dependancies and particularly Python3.7 and SymPy==1.1"
echo "This script mush be excuted from the root directory of the OpenSBLI repository"
echo "The script can be take the following optional inputs   "
echo "  1=> Set up OSBLI mode:                               "
echo "      - ENV : Set up the OSBLI enviroment only (default)"
echo "      - APP : Set up the OSBLI enviroment and build apps"
echo "      - TEST: Set up the OSBLI enviroment and run tests "
echo "  2=> Path to HDF5 installation (i.e. HDF5 bin location)"
echo "  3=> Path to OPS installation  (i.e. OPS bin location)"
echo "If inputs are not provided HDF5 and OPS are installed"
echo "In case the virtual enviroment has already been set up HDF5 and OPS are only activated"
echo "------------------------------------------------------------------------------------- "
if [ -z "$1" ];
then
  echo "No inputs provided HDF5 and OPS will be installed if not already present"
  mode=ENV
elif [ -z "$2" ];
then
  mode=$1
elif [ -z "$3" ];
then
  mode=$1
  path_hdf5=$2
  export HDF5_INSTALL_PATH=${path_hdf5}
  echo "HDF5 is provided and this will be used $HDF5_INSTALL_PATH"
elif [ -z "$4" ];	
then
  mode=$1
  path_hdf5=$2
  path_ops=$3
  export HDF5_INSTALL_PATH=${path_hdf5}
  export OPS_INSTALL_DIR=${path_ops} 
  echo "HDF5 is provided and this will be used $HDF5_INSTALL_PATH"
  echo "OPS is provided and this will be used $OPS_INSTALL_DIR"
else
  echo "Too many inputs we exit"
 exit
fi

cd opensbli_venv/
source setup_venv.sh 
echo "My activate command ${python_activate}"
echo "My HDF5 path ${HDF5_INSTALL_PATH}"
echo "My OPS  path ${OPS_INSTALL_DIR}"
cd ../
if [ "$mode" = "APP" ];
then	
  cd apps
  python generate_all_applications.py
  cd ..
elif [ "$mode" = "TEST" ];
then	
  cd tests
  python test_opensbli.py --verif-only 
  cd ..
fi 
  

