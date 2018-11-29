# set tabstop=8 softtabstop=0 expandtab shiftwidth=4 smarttab
from Cython.Build import cythonize
from Cython.Distutils import Extension
# Change this commented line with the following one gives us the possibility
# to pass \"cython_compile_time_env\" at \"Extension\" because in line 154 we do
# not have more to call \"cythonize\" function before returning \"ext_modules\".
#from distutils.command.build_ext import build_ext as _build_ext
from Cython.Distutils import build_ext as _build_ext
from distutils.core import setup

from numpy import import get_include as numpy_get_include

import imp
import os
import re
import sys

try:
    imp.find_module('mpi4py')
    ENABLE_MPI4PY = 1
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

    # Add new \"user_options\" \"BitP_Mesh-lib-path\", \"mpi-library-path\",
    # \"extensions-source\" and \"BitP_Base-lib-path\".
    user_options.append(("bitpit-include-path=", "P", "bitpit include path"))
    user_options.append(("madeleine-include-path=", "I", "madeleine include path"))
    user_options.append(("mpi-include-path=", "M", "mpi include path"))
    user_options.append(("extensions-source=", "E", "extensions source file"))
