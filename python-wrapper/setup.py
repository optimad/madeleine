# set tabstop=8 softtabstop=0 expandtab shiftwidth=4 smarttab
from Cython.Build import cythonize
from Cython.Distutils import Extension
# Change this commented line with the following one gives us the possibility
# to pass \"cython_compile_time_env\" at \"Extension\" because in line 154 we do
# not have more to call \"cythonize\" function before returning \"ext_modules\".
#from distutils.command.build_ext import build_ext as _build_ext
from Cython.Distutils import build_ext as _build_ext
from distutils.core import setup

from numpy import get_include as numpy_get_include

import imp
import os
import re
import sys

try:
    imp.find_module('mpi4py')
    ENABLE_MPI4PY = 1
    BITPIT_ENABLE_MPI = 1
    ENABLE_MPI = 1
except ImportError:
    ENABLE_MPI4PY = 0
    print("No module \"mpi4py\" found in your environment. We are " +
          "going serial in Python.")

class build_ext(_build_ext):
    description = ("Custom build_ext \"bitpit-include-path\" "   +
                   " \"mpi-include-path\" \"madeleine-include-path\" " +
                   " \"extension-source\" command for Cython")

    # Getting \"user_options\" from class \"build_ext\" imported as 
    # \"_build_ext\".
    user_options = _build_ext.user_options

    # Add new \"user_options\" \"bitpit-lib-path\", \"mpi-library-path\",
    # \"extensions-source\" and \"bitpit-lib-path\".
    user_options.append(("bitpit-path=", "P", "bitpit path"))
    user_options.append(("madeleine-path=", "I", "madeleine path"))
    user_options.append(("mpi-include-path=", "M", "mpi include path"))
    user_options.append(("metis-path=", "N", "metis path"))
    user_options.append(("petsc-path=", "W", "petsc path"))
    user_options.append(("lapack-path=", "L", "lapack path"))
    user_options.append(("extensions-source=", "E", "extensions source file"))

    def find_mpi_include_path(self):
        LIBRARY_PATHS = os.environ.get("LD_LIBRARY_PATH")
        mpi_checked = False
        MPI_INCLUDE_PATH = None

        for LIBRARY_PATH in LIBRARY_PATHS.split(":"):
            if (("mpi" in LIBRARY_PATH.lower()) and 
                (not mpi_checked)):
                MPI_INCLUDE_PATH = LIBRARY_PATH
                if not "/include/" in MPI_INCLUDE_PATH:
                    MPI_INCLUDE_PATH = re.sub("/lib?"    , 
                                              "/include/",
                                              MPI_INCLUDE_PATH)
                
                mpi_checked = True
                break
    
        if (not mpi_checked):
            print("No \"mpi-include-path\" found in your " +
                  "\"LD_LIBRARY_PATH\" environment variable. "   +
                  "We are going serial in C/C++.")
            
        return MPI_INCLUDE_PATH


    def find_bitpit_path(self):
        BITPIT_PATH = os.environ.get("BITPIT_PATH")

        if (BITPIT_PATH is None):
            print("No \"BITPIT_INCLUDE_PATH\" env variable found. " +
                  "Please, check this out or enter it via shell.")

            sys.exit(1)

        return BITPIT_PATH


    def find_madeleine_path(self):
        MADELEINE_PATH = os.environ.get("MADELEINE_PATH")

        if (MADELEINE_PATH is None):
            print("Dude, no \"MADELEINE_PATH\" env variable found. Please, " + 
                  "check this out or enter it via shell.")

            sys.exit(1)

        return MADELEINE_PATH

    def find_metis_path(self):
        METIS_PATH = os.environ.get("METIS_PATH")

        if (METIS_PATH is None):
            print("Dude, no \"METIS_PATH\" env variable found. Please, " + 
                  "check this out or enter it via shell.")

            sys.exit(1)

        return METIS_PATH

    def find_petsc_path(self):
        PETSC_PATH = os.environ.get("PETSC_PATH")

        if (PETSC_PATH is None):
            print("Dude, no \"PETSC_PATH\" env variable found. Please, " + 
                  "check this out or enter it via shell.")

            sys.exit(1)

        return PETSC_PATH

    def find_lapack_path(self):
        LAPACK_PATH = os.environ.get("LAPACK_PATH")

        if (LAPACK_PATH is None):
            print("Dude, no \"LAPACK_PATH\" env variable found. Please, " + 
                  "check this out or enter it via shell.")

            sys.exit(1)

        return LAPACK_PATH


    def check_extensions_source(self):
        if ((self.extensions_source is None) or 
            (not self.extensions_source.endswith(".pyx"))):
            print("Insert source \".pyx\" file to build.")
        
            sys.exit(1)


    def initialize_options(self):
        # Initialize father's \"user_options\".
        _build_ext.initialize_options(self)

        # Initializing own new \"user_options\".
        self.bitpit_path = None
        self.madeleine_path = None
        self.metis_path = None
        self.petsc_path = None
        self.lapack_path = None
        self.mpi_include_path = None
        self.extensions_source = None


    def finalize_options(self):
        # Finalizing father's \"user_options\".
        _build_ext.finalize_options(self)
        
        # If yet \"None\", finalize own new \"user_options\" searching their
        # values.
        if (self.mpi_include_path is None):
            self.mpi_include_path = self.find_mpi_include_path()
        if (self.bitpit_path is None):
            self.bitpit_path = self.find_bipit_path()
        if (self.madeleine_path is None):
            self.madeleine_path = self.find_madeleine_path()
        if (self.metis_path is None):
            self.metis_path = self.find_metis_path()
        if (self.petsc_path is None):
            self.petsc_path = self.find_petsc_path()
        if (self.lapack_path is None):
            self.lapack_path = self.find_lapack_path()

        # Check if the source to pass at the \"Extension\" class is present and
        # finishes with \".pyx\".
        self.check_extensions_source()
        # Define \"custom cython\" extensions.
        self.extensions = self.def_ext_modules()


    def def_ext_modules(self):
        # Define \"Extension\" being cythonized.
        # Overloading compilers.
        os.environ["CXX"] = "c++"
        os.environ["CC"] = "gcc"
        BITPIT_ENABLE_MPI = 0
        ENABLE_MPI = 1
        include_paths = [self.bitpit_path + "/include/bitpit/", self.madeleine_path + "../../src/"]

        mpi_lib = ""
        if ((not (not self.mpi_include_path)) and (ENABLE_MPI4PY)):
            BITPIT_ENABLE_MPI = 1
            ENABLE_MPI = 1
            include_paths.append(self.mpi_include_path)
            os.environ["CXX"] = "mpic++"
            os.environ["CC"] = "mpic++"
            mpi_lib = re.sub("/include/","/lib/",self.mpi_include_path) + "libmpi.so"
            #mpicxx_lib = re.sub("/include/","/lib/",self.mpi_include_path) + "libmpi_cxx.so"
            print(mpi_lib)

        _extra_compile_args = ["-std=c++11",
                               "-g"        ,
                               "-O0"       ,
                               "-fPIC"     ,
                               "-DVERBOSE=1",
                               #include_paths,
                               "-DBITPIT_ENABLE_MPI=" + str(BITPIT_ENABLE_MPI),
                               "-DENABLE_MPI=" + str(ENABLE_MPI)]
        _extra_link_args = ["-fPIC"] # Needed? We have already the same flag for 
                                     # the compiler args above.
        _cython_directives = {"boundscheck": False,
                              "wraparound": False,
                              # http://stackoverflow.com/questions/23351813/how-to-declare-an-ndarray-in-cython-with-a-general-floating-point-type
                              "nonecheck": False}
        _language = "c++"
        bitpit_lib = self.bitpit_path + "/lib/"
        if(BITPIT_ENABLE_MPI):
            bitpit_lib = bitpit_lib + "libbitpit_MPI_D.so"
        else:
            bitpit_lib = bitpit_lib + "libbitpit_D.so"
        madeleine_lib = self.madeleine_path + "libmadeleine_MPI_D.so"
        metis_lib = self.metis_path + "libmetis.a"
        petsc_lib = self.petsc_path + "libpetsc.so"
        lapack_lib = self.lapack_path + "liblapacke.so"
        _extra_objects = [madeleine_lib, bitpit_lib, petsc_lib, mpi_lib, metis_lib, lapack_lib, "-lxml2"]

##FINO QUI!
        #print(os.path.dirname(self.bitpit_include_path))
        #piercedVector_subdir = os.path.dirname(self.bitpit_include_path) + "/containers/"
        _include_dirs=["."                        , 
                       numpy_get_include()        ]
        _include_dirs.extend(include_paths) 
       
        # Cython compile time environment.
        _cc_time_env = {"BITPIT_ENABLE_MPI": BITPIT_ENABLE_MPI}
    
        ext_modules = [Extension(os.path.splitext(self.extensions_source)[0],
                                [self.extensions_source]                   ,              
                                extra_compile_args = _extra_compile_args   ,
                                extra_link_args = _extra_link_args         ,
                                cython_directives = _cython_directives     , 
                                language = 'c++'                       ,
                                extra_objects = _extra_objects             ,
                                include_dirs = _include_dirs               ,
                                cython_compile_time_env = _cc_time_env     ,
                    )]
        return cythonize(ext_modules,gdb_debug=True,verbose=True)
        #return ext_modules


    def run(self):
        _build_ext.run(self)
        

setup(cmdclass = {"build_ext" : build_ext})
