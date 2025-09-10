#!/bin/bash
local_dir=`pwd`
# Check version of available Python3
PYTHON_VERSION=$(python3 --version | cut -d ' ' -f 2)
# Display the version
echo "Python version is: $PYTHON_VERSION"
# Example usage of the variable in a conditional statement
if [[ "$PYTHON_VERSION" =~ ^3\.[0-9]+\.[0-9]+ ]]; then
  echo "You are using Python 3"
  # Extract major, minor, and patch versions
  MAJOR_VERSION=$(echo $PYTHON_VERSION | cut -d. -f1)
  MINOR_VERSION=$(echo $PYTHON_VERSION | cut -d. -f2)
  PATCH_VERSION=$(echo $PYTHON_VERSION | cut -d. -f3)
  echo "Major version: $MAJOR_VERSION"
  echo "Minor version: $MINOR_VERSION"
  echo "Patch version: $PATCH_VERSION"
  if [ "$MINOR_VERSION" -lt "8" ]; then
    echo "Python3 must be version 8 or above for OPS to work"
    echo "We are exiting "
    return 1
  fi
else
  echo "You are not using Python 3. needed for OSBLI to work"
  echo "We are exiting "
  return 1
fi
# Test for python 3.7
PY37=$(python3.7 --version | cut -d ' ' -f 2)
if [ -z "${PY37}" ]; then
  py37install=${local_dir}/python37/opt
  py37=${py37install}/bin/python3.7
  if ! [ -f ${py37} ]; then
    echo "No Python3.7 a Miniconda version will be used"
    mkdir python37
    cd python37
    PackageName="Miniconda3-py37_22.11.1-1-Linux-x86_64.sh" 
    wget -c https://repo.anaconda.com/miniconda/${PackageName}
    chmod a+x ./${PackageName}
    ./${PackageName} -b -p ${py37install}
    cd ../
  else
    echo "Python3.7 is installed locally via Miniconda"
  fi
else
  echo "Python3.7 is available"
  py37="python3.7"
fi
# Now we need to create the virtual enviroment
export osbli_venv_dir=${local_dir}/osbli_opt
export osbli_venv_activate=${osbli_venv_dir}/bin/activate
#if ! [ -f ${osbli_venv_activate} ]; then
#  echo "Virtual Enviroment do not exists we need to create it"
#  python3 -m pip install --user virtualenv
#  virtualenv -p ${py37} ${osbli_venv_dir}
#fi
# Activate the Virtual Enviroment
echo "My activate command call ${osbli_venv_activate}"
source ${osbli_venv_activate}
# Update Python3.7 packages only in case of non Conda download
#if [ -z "${py37install}/bin/python3.7" ]; then
# Not strictly necessary but we can make sure that all libs are available
python3 -m pip install --upgrade pip
python3 -m pip install -r ${local_dir}/requirements.txt
#fi
if [ ! -f osbli_env.sh ]; then
  echo "Myosbli_venv_activate HDF5 ENV PATH ${HDF5_INSTALL_PATH}" 
  cmake -S . -B build 
fi
echo "SOURCE OSBLI ENV"
source osbli_env.sh
