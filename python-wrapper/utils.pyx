include "piercedVector.pyx"
include "surfUnstructured.pyx"
from libc.stdint cimport uintptr_t

cdef extern from "couplingUtils.hpp" namespace "coupling":
    cdef void initDoubleDataOnMesh(const SurfUnstructured * mesh, PiercedVector[double,long] * data)
    
def Py_initDoubleDataOnMesh(Py_SharedSurfUnstructured mesh, Py_PiercedVector data):
    initDoubleDataOnMesh(<SurfUnstructured*><void*>mesh.thisptr,<PiercedVector[double,long]*><void*>data.thisptr)
