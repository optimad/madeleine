cdef extern from "surfunstructured.hpp" namespace "bitpit":
    cdef cppclass SurfUnstructured:
        SurfUnstructured() except +
    
