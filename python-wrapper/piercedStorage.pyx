import traceback
from libc.stdint cimport uintptr_t

cdef extern from "piercedStorage.hpp" namespace "bitpit":
    cdef cppclass PiercedStorage[value_t,id_t]:
        PiercedStorage() except +
        PiercedStorage(const PiercedStorage[value_t, id_t] &x) except +
        PiercedStorage[value_t, id_t] & operator=(PiercedStorage[value_t, id_t] x)
        
    
cdef class Py_PiercedStorage:
    cdef PiercedStorage[double,long]* thisptr
    
    def __cinit__(self,*args):
        
        cdef uintptr_t int_ptr = 0
        n_args = len(args)
        
        if(n_args == 0):
            self.thisptr = new PiercedStorage[double,long]()
        elif(n_args == 1):
            try:
                int_ptr = args[0]
                self.thisptr = new PiercedStorage[double,long]((<PiercedStorage[double,long]*><void*>int_ptr)[0])
            except Exception as e:
                traceback.print_exc()
        else:
            print("Py_PiercedStorage, wrong number of arguments!")
            
    def __dealloc__(self):
        del self.thisptr
