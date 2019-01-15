import traceback
from libc.stdint cimport uintptr_t

cdef extern from "surfunstructured.hpp" namespace "bitpit":
    cdef cppclass SurfUnstructured:
        SurfUnstructured() except +
    
cdef class Py_SurfUnstructured:
    cdef SurfUnstructured* thisptr
    
    def __cinit__(self,*args):
        
        cdef uintptr_t int_ptr = 0
        n_args = len(args)
        
        if(n_args == 0):
            self.thisptr = new SurfUnstructured()
        else:
            print("Py_SurfUnstructured, wrong number of arguments!")

    def __dealloc__(self):
        del self.thisptr


cdef class Py_SharedSurfUnstructured:
    cdef SurfUnstructured* thisptr

    def __cinit__(self,*args):
        cdef uintptr_t int_ptr = 0
        n_args = len(args)
        if(n_args == 1):
            int_ptr = args[0]
            self.thisptr = (<SurfUnstructured*><void*>int_ptr)
        else:
            print("Py_SurfUnstructured, wrong number of arguments!")
