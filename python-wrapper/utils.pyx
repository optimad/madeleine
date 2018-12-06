include "piercedVector.pyx"
include "surfUnstructured.pyx"
from libc.stdint cimport uintptr_t

cdef extern from "couplingUtils.hpp" namespace "coupling":
    cdef void initDoubleDataOnMesh(const SurfUnstructured * mesh, PiercedVector[double,long] * data)
    
def Py_initDoubleDataOnMesh(Py_SurfUnstructured mesh, Py_PiercedVector data):
    initDoubleDataOnMesh(mesh.thisptr,data.thisptr)
