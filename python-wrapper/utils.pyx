include "piercedVector.pyx"
include "surfUnstructured.pyx"

cdef extern from "utils.hpp" namespace "coupling":
    PiercedVector[double,long] initDoubleDataOnMesh(const SurfUnstructured * mesh)