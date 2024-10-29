## Run the Tests
In order to perform the tests some python dependacies are needed. 
To take care of these a Python Virtual Enviroment has been set up under the 
folder `opensbli_venv/`.
To run the tests the bash script `RunTests.sh` located in the root directory of 
this git reporitory. 
The script is installing and activate the Python virtual enviroment, export 
the relevant enviromental variables to locate OpenSBLI and OPS and proceed 
with the code translation, compililation and running. 
To run the script some enviromental variables needs to be set together 
with two inputs with the location of HDF5 and OPS

```
$ export CXX=mpicxx  
$ export CC=mpicc
$ ./RunTests.sh /path/to/HDF5 /path/to/OPS
```  
