# Installation 
This guide provides instructions on how to install and use the OpenSBLI code generator framework to build up APPS. 
The following dependencies are required for generating a code from Python into C++ based on the 
OPS (Oxford Parallel library for Structured mesh solvers) domain specific language (DSL) framework 

* Python 3.7
* Numpy
* Python-utils
* Cached-property
* Pytest-shutil
* Sympy 1.1
* h5py
* Scipy
* Matplotlib

A Python virtual enviroment has been created under [opensbli_venv](opensbli_venv) 
and It can be activated by sourcing the script 
```
$ cd opensbli_venv
$ source setup_venv.sh
```
All relevant libraries will be installed under `opensbli_venv/py37_opt`.
To be sure that you are working under the `python_venv` the propt of your terminal 
should  looks like
```
(py37_opt) $  cd ../
``` 
OpenSBLI requires also a working version of the HDF5 library installed with MPI support 
and a working version of the development branch of OPS defined as
```
(py37_opt) $ export OPS_INSTALL_DIR=path/to/OPS/opt 
(py37_opt) $ export OPS_TRANSLATOR=$OPS_INSTALL_DIR/translator/ops_translator/ops-translator
(py37_opt) $ export PYTHONPATH=$PYTHONPATH:path/to/opensbli/repository
(py37_opt) $ export HDF5_INSTALL_PATH=path/to/hdf5/mpi
```
A set-up-script to install, activate the virtual enviroment and set up the enviromental variables 
for OpenSBLI is available as [SetUpEnviroment.sh](SetUpEnviroment.sh) and can be activated by
```
(py37_opt) $ source SetUpEnviroment.sh path/to/hdf5/mpi path/to/OPS/opt   
``` 
where the first input is the path to the HDF5 install directory and the second is the path 
to the OPS installation. 

Work is on progress to have an automatic installation of HDF5 and OPS requirements in case the
requirements are not found. 

## Automatic Apps Installation
It is possible to generate and compile all [apps](apps) using the python script 
[generate_all_applications.py](apps/generate_all_applications.py) as:
```
(py37_opt) $ cd ../apps
(py37_opt) $ python generate_all_applications.py 
``` 
that can take the following inputs:
```
(py37_opt) rfj82982@ccp-gpu2:apps$ python generate_all_applications.py --help 
usage: generate_all_applications.py [-h] [--legacy-translator] [--verbose]
                                    [--generate] [--app APP] [--target TARGET]

optional arguments:
  -h, --help           show this help message and exit
  --legacy-translator
  --verbose
  --generate
  --app APP
  --target TARGET
``` 
where `--app`  triggers the build of a single app like [wave](app/wave) 
and `--target` restricts the build of few targets and can take the 
values `mpi`, `CUDA`, `seq` or `all`. 
The default name for the workspace directory is **apps/_opensbli-build-workspace** where all
apps are build. 
Inside every *app* directory there is the following struture: 

* the *app.py* original OpenSBLI python script
* the *opensbli.cpp* code generated file together with all relevant include files
* the *CMakeList.txt* for app build
* the **test-build** folder with the app  build done with the CMake system

**P.S.** A *Makefile* is also available in the working directory, however it is not guarantee 
to work unless OPS has been built using the traditional make procedure. The procedure here 
described assumes a CMake build for both OPS and OpenSBLI app. 

A script is available to perform the activation of the python virtual enviroment, set up of the
OpenSBLI variables and run of the python app generation and can be used as 
```
(py37_opt) $ source GenerateApps.sh path/to/hdf5/mpi path/to/OPS/opt   
``` 

### App modification and rebuild
It is possible to modify and regenerate/rebuild a specific app. To do so two possibilities are available
and are here described. In both cases it is assumed that the worrking enviroment is the generated 
working directory located 
```
(py37_opt) $ cd apps/_opensbli-build-workspace**/app   
``` 
where **app** is the app that we are trying to build like [wave](app/wave)

1. Direct modification of the *opensbli.cpp* source file. In this case only the OPS translation needs to be 
redone followed by the build. This can be done as follows: 
```
(py37_opt) $ cd apps/_opensbli-build-workspace**/app/test-build
(py37_opt) $ cmake --build . --target clean ! remove the already generated targets
(py37_opt) $ cmake ../
(py37_opt) $ cmake --build . -j 2 ! or more for parallel build   
``` 
1. Modification of the OpenSBLI python file *app.py*. In this case am OpenSBLI code translations needs to done by
```
(py37_opt) $ cd apps/_opensbli-build-workspace**/app
(py37_opt) $ python app.py
``` 
**PS** The OpenSBLI python code generation is guarantee to work only within the OpenSBLI python virtual enviroment. 
Please see section [install](#Install) on how to activate it. 


















