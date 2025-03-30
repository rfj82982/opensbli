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

The code translation to enable the automatic parallelisation on different hybrid computing architectures
is performed using OPS with the possibility to have parallel HDF5 for data postprocessing. 

To have all requirements together a virtual enviroment is created and can be activated using [SetUpOSBLI.sh](SetUpOSBLI.sh). 
The script is creating and activating the correct Python virtual enviroment. 
HDF5 and OPS are optional inputs and, if not provided, are installed in the same virtual enviroment 
that is located at *opensbli_venv/opensbli_opt* 
The virtual enviroment can be activated by sourcing the script 
```
$ source SetUpOSBLI.sh
```
To be sure that you are working under the virtual enviroment is using the correct Python 3.7
the terminal prompt should looks like
```
(osbli_opt) $  python --version
Python 3.7.17
``` 
The result of the set-up bash script is the installation of the required libraries and the
set up the enviromental variables necessary for OpenSBLI to do the code generation 
and OPS to do the code translation. 
The necessary enviromental variables are:
```
(osbli_opt) $ export PYTHONPATH=$PYTHONPATH:path/to/opensbli/repository
(osbli_opt) $ export HDF5_INSTALL_PATH=path/to/hdf5/mpi
(osbli_opt) $ export OPS_INSTALL_DIR=path/to/OPS/opt 
(osbli_opt) $ export OPS_TRANSLATOR=$OPS_INSTALL_DIR/translator/ops_translator/ops-translator
```
These are stored in the file *opensbli_venv/osbli_env.sh*. The set-up script is always checking 
at every invocation if the OpenSBLI enviromental variables are defined to avoid constant 
reinstall HDF5 and OPS. 
The final structure of the OpenSBLI virtual enviroment is 
```
(osbli_opt) $ ls opensbli_venv/
CMakeLists.txt  build  cmake  osbli_env.sh  osbli_opt  requirements.txt  setup_venv.sh
``` 
where 

* *CMakeLists.txt*: CMake main file for download and build HDF5 and OPS
* **build**: working directory for the additional library build
* **cmake**: support CMake files for configure and build
* *osbli_env.sh*: file with enviromental variables for OpenSBLI and OPS
* **osbli_opt**: installtion directory for the virtual enviroment
* *requirements.txt*: list of requirements for the python virtual enviroment
* *setup.venv.sh*: set-up script for the virtual enviroment 

## About Python requirements
The OpenSBLI workflow requires two different Python versions: 

* **Python3.7** for the OpenSBLI code generation
* **Python3.8** or above for the OPS code translation

The handling of Python3.7 is done via the OpenSBLI virtual enviroment. 
If a Python3.7 is available on the system this will be used and the   
additional library requirements will be installed. If no Python3.7 is available
a MiniConda is downlaoded and installed. 
The requirement for Python3.8 is instead managed by the OPS virtual enviroment. 

## Automatic Apps Installation
It is possible to generate and compile all [apps](apps) using the python script 
[generate_all_applications.py](apps/generate_all_applications.py) as:
```
(osbli_opt) $ cd ../apps
(osbli_opt) $ python generate_all_applications.py 
``` 
that can take the following inputs:
```
(osbli_opt) $ python generate_all_applications.py --help 
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
(osbli_opt) $ source SetUpOSBLI.sh TEST  
``` 

### App modification and rebuild
It is possible to modify and regenerate/rebuild a specific app. To do so two possibilities are available
and they are here described. In both cases it is assumed that the worrking enviroment is the generated 
working directory located 
```
(osbli_opt) $ cd apps/_opensbli-build-workspace**/app   
``` 
where **app** is the app that we are trying to build like [wave](app/wave)

1. Direct modification of the *opensbli.cpp* source file. In this case only the OPS translation step needs to be 
performed, which is followed by the build step. This can be done as follows: 
```
(osbli_opt) $ cd apps/_opensbli-build-workspace**/app/test-build
(osbli_opt) $ cmake --build . --target clean ! remove the already generated targets
(osbli_opt) $ cmake ../
(osbli_opt) $ cmake --build . -j 2 ! or more for parallel build   
``` 
1. Modification of the OpenSBLI python file *app.py*. In this case an OpenSBLI code translations needs to be
performed 
```
(osbli_opt) $ cd apps/_opensbli-build-workspace**/app
(osbli_opt) $ osbli_code_generation app.py
``` 
and followed by the build as decribed in the above point.  

### App build from scratch 
The OpenSBLI virtual enviroment provides also commands to build an apps in a step to step fashion. 
This is useful in case of building new cases starting from a new *newcase.py* python file. 

1. *osbli_code_generation*: command to generate the OPS C++ file from the file *newcase.py*. Remember that 
   the command always requires am open_sbli python input file. The result is the creation of a *opensbli.cpp*
   file and relative includes files.
1. *ops_translation*: command to run the OPS python code translation. The input of the command is always the file *opensbli.cpp*
1. *osbli_cmake_command*: this command creates the CMake configure to compile the fifferent backends for the app. 
   The CMake build system will be created under the folder **test-build**

Please remember that afterwards the apps needs to be build by using: 
```
(osbli_opt) $ cmake --build test-build -j ...
```      

An example of how to build an app step-by-step is the following: 
```
(osbli_opt) $ mkdir mytestapp
(osbli_opt) $ cd mytestapp
(osbli_opt) $ cp path_to_opensbli/apps/wave/* . 
(osbli_opt) $ osbli_code_generation wave.py
(ops_venv) $ ops_translation
(ops_venv) $ osbli_cmake_command
(ops_venv) $ cmake --build test-build -j 4 
(ops_venv) $ ls test-build/OpenSBLI_* 
test-build/OpenSBLI_cuda     test-build/OpenSBLI_mpi             test-build/OpenSBLI_mpi_openmp  test-build/OpenSBLI_seq
test-build/OpenSBLI_dev_mpi  test-build/OpenSBLI_mpi_cuda        test-build/OpenSBLI_mpi_tiled   test-build/OpenSBLI_tiled
test-build/OpenSBLI_dev_seq  test-build/OpenSBLI_mpi_cuda_tiled  test-build/OpenSBLI_openmp
```
If you would like to build less targets you can use the CMake variable `OPS_TARGET`. This can be modified using *ccmake* as
```
cmake --build test-build --target clean 
ccmake test-build
```

## Testing
A testing framework is also available for OpenSBLI. This is based on the [test_opensbli.py](tests/test_opensbli.py) 
python script. 
The code generation, OPS translation and CMake build follow the structure given in the above sections. 
In this case the workspace is located under **tests/_opensbli-test-workspac** and the full output is 
recorded in *tests/test.log*. In case of a fully successful test run the working directory is deleted. 
To run the full test from scratch with the activation of the OpenSBLI virtual enviroment use: 
```
$ source SetUpOSBLI.sh TEST
```

## Usage of external build of HDF5 and OPS
System installation of HDF5 and external OPS can be used and the location path can be passed to the 
set-up script as:
```
$ source SetUpOSBLI.sh ENV /path/to/root/HDF5 /path/to/root/OPS
```
with the caveat that the OPS has to be built with the same HDF5. 
 
