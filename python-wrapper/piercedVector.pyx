import traceback
from libc.stdint cimport uintptr_t

cdef extern from "piercedVector.hpp" namespace "bitpit":
    cdef cppclass PiercedVector[value_t,id_t]:
        PiercedVector() except +
        PiercedVector(const PiercedVector<value_t, id_t> &x) except +
        PiercedVector<value_t, id_t> & operator=(PiercedVector<value_t, id_t> x)
        
    
cdef class Py_PiercedVector:
    cdef PiercedVector[double,long]* thisptr
    
    def __cinit__(self,*args):
        
        cdef uintptr_t int_ptr = 0
        n_args = len(args)
        
        if(n_args == 0):
            self.thisptr = new PiercedVector[double,long]()
        elif:
            try:
                int_ptr = args[0]
                self.thisptr = new PiercedVector[double,long]((<PiercedVector[double]*><void*>int_ptr)[0])
            except Exception as e:
                traceback.print_exc()
        else:
            print("Py_PiercedVector, wrong number of arguments!")
            
    def __dealloc__(self):
        del self.thisptr