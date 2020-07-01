#MADELEINE Python Wrapper - The CYTHON guide

##Conda Environment
In yml file ***finaleEnvironment.yml*** in the Madeleine root folder all the conda and pip packages to run this project and GEMS are listed. Besically, new mandatory packages are ***lapacke*** and ***cblas***

##C++ Coupling library
It is important that the ***coupling.so*** library is compiled using the conda environment, to be compatible with GEMS. To have a coherent compilation, ***bitpit*** library, which coupling.so depends on, must be compiled using the same environment and as a shared library.

Once the right conda environment is activated and the bash environment is prepared (see below), cmake should be able to find the right compiler and libraries from the conda environment for both couplng.so and bitpit compilation. 

##Bash Environment
As usual, the bash environment must be compliant to the conda one. The user must export bash variables:
- PATH
- LD_LIBRARY_PATH
- PYTHONPATH

to have the correct paths to the conda folders where libraries, headers and python modules are.

Let's say:
 - the conda environment path is in bash variable ***CONDA_ENV_PATH***,
 - conda installation path is in ***CONDA_PATH***,
 - bitpit installation is in ***BITPIT_PATH*** and
 - the folder containing madeleine Python module coupling.so is in ***MADELEINE_WRAPPER_PATH***
 - the folder containing the C++ madeleine library is ***MADELEINE_PATH***
 - the GEMS folder is in ***GEMS_PATH***
 
 then:
- ***PATH*** should contain at first positions:
    - CONDA_ENV_PATH/bin and
    - CONDA_PATH/condabin
- ***LD_LIBRARY_PATH*** should contain at first positions:
    - CONDA_ENV_PATH/lib
    - BITPIT_PATH/lib
    - MADELEINE_PATH
- ***PYTHONPATH*** should contain at first positions:
    - GEMS_PATH
    - MADELEINE_WRAPPER_PATH
    
See exportVariable script for example on my machine. Customizing ***CONDA_ENV_PATH***, ***CONDA_PATH***, ***BITPIT_PATH***, ***MADELEINE_WRAPPER_PATH***, ***MADELEINE_PATH***, ***GEMS_PATH*** in *exportVariables* file and launching
 
`source ./exportVariables`

prepares the bash environment for python wrapper usage.

##Cython arguments
In order to produce the ***coupling*** python module, the cythonization of the python interface to the C++ madeleine library must be performed.

Basically, this is achieved by launching

`python setup.py build_ext --inplace [...] --extensions-source=coupling.pyx`

in ***$MADELEINE_WRAPPER_PATH*** folder, where `[...]` is a placeholder for the following arguments. If the bash environment has been prepared as above, exported variables can be used (**pay attention to the trailing slash**):

- *--bitpit-path*=$BITPIT_PATH/
- *--madeleine-path*=$MADELEINE_PATH/
- *--petsc-path*=$CONDA_ENV_PATH/lib/
- *--lapack-path*=$CONDA_ENV_PATH/lib/
- *--mpi-include-path*=$CONDA_ENV_PATH/include/
<br/><br/><br/><br/>
A special attention should be paid for
- *--metis-path*=[metis-installation-folder]

at the moment metis is not used but the library is ready to exploit its algorithms. The argument is needed and a working metis library is mandatory. In the future, some improvements will be introduced and metis should be installed into the conda environment. For now, any metis installation can be used for the sake of compilation.