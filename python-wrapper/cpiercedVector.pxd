cdef extern from "piercedVector.hpp" namespace "bitpit":
    cdef cppclass PiercedVector[value_t,id_t]:
        PiercedVector() except +
        PiercedVector(const PiercedVector[value_t, id_t] &x) except +
        PiercedVector[value_t, id_t] & operator=(PiercedVector[value_t, id_t] x)
        
    
