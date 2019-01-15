# madeleine interpolation tool (MIT)

MIT is a simple C++ interpolation tool based on bitpit with Python bindings.
MIT provides user with a C++ library for interpolation between two meshes (defined as discipline meshes) via a buffer mesh (defined neutral mesh).
This library can be used in Python thanks to a Cythonization of the C++ library.

##### C++ library dependencies:
- cmake
- bitpit (compiled as a shared library with flag "-fPIC")
- MPI (optional)

##### Python wrapping dependencies:
- Cython
- Numpy
- other Python modules: distutils, libc, libcpp, traceback, imp, os, re, sys

##### Configuring C++ library compilation:
In MIT's root folder make a building folder, e.g. build
```bash
    madeleine$ mkdir build
```
Enter the `build` folder
```bash
    madeleine$ cd build
```
 In order to configure it, run:
```bash
    madeleine/build$ ccmake ../
```
By this way, MIT can be configured for C++ and/or Python usage.
Setting some variable in ccmake interface you can customize a bit your configuration.
The only mandatory cmake variable that has to be set is 
```bash
    BITPIT_DIR = [/bitpit/installation/dir]/lib/cmake/bitpit-{bitpit_version}
```
where bitpit-version should be 1.5 (previous versions have not been tested).
In order to use the Python wrappers the following cmake variable has to be set:
```bash
    BUILD_SHARED_LIBS = ON
```
C++ examples are compiled if the following bariable is set
```bash
    BUILD_EXAMPLES = ON
```
##### Building C++ library:
Once configured the C++ library can be compiled with
```bash
    make
```

##### C++ library example:
If ```BUILD_EXAMPLES``` variable has been set to ```ON```, in ```[build_folder]/example``` an exacutable named ```example_00002``` is produced.

Just move to the above folder and launch it as follows
```bash
    ./example_00002
```
By this way the example will use to default meshes and analtical data to perform interpolation from neutral mesh to discipline mesh, a very easy computation and an interpolation from discipline mesh to neutral mesh.
The example procudes vtu files to check the result.

The user can customize the example run by providing 2 stl files containing 2 different triagulations of a unit sphere.
```bash
    ./example_00002 /path/to/first/stl/file /path/to/second/stl/file
```

##### Building Python wrappers:
Once the MIT C++ library has been compiled, in order to build Cythonized Python wrappers, move to the
```bash
 python-wrapper
```
folder into the MIT's root folder and launch
```bash
    python setup.py build_ext --inplace --bitpit-path=[/path/to/bitpit/installation/folder] \
    --madeleine-path=[/path/to/madeleine/root/folder]/build/src/ --extensions-source=coupling.pyx
```
This string will produce a Python module named
```bash
    coupling.so
```
In order to be able to use it the following environment variable should contain at least
```bash
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:[/path/to/bitpit/dynamic/library]:[/path/madeleine/root/folder]/build/src
```
Optionally
```bash
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:[/path/to/bitpit/dynamic/library]:[/path/madeleine/root/folder]/build/src:[/path/to/MPI/library]
```

##### Running Python wrappers example:
In the ```python-wrapper``` folder, an example, named ```py_example_00002.py``` is provided. Just simply run
```bash
    python py_example_00002.py
```
to launch it. This script reproduces exactly the same computation of the C++ example (```example_00002```)

