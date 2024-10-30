#!/bin/bash

#Do not run below command if virtualenv is already installed
local_dir=`pwd`
echo $local_dir
python_dir=${local_dir}/py37_opt
export python_activate=${python_dir}/bin/activate
if ! [ -f ${python_activate} ]; then
  echo "Python3.7 Virtual Enviroment do not exists we need to create it"
  python3 -m pip install --user virtualenv
  virtualenv -p python3.7 ${python_dir}
fi
# Activate the Virtual Enviroment
source ${python_activate}
# Not strictly necessary but we can make sure that all libs are available
python3 -m pip install --upgrade pip
python3 -m pip install -r ${local_dir}/requirements.txt
echo "My activate command call ${python_activate}"

