#include "piercedVector.pyx"
#include "surfUnstructured.pyx"

cimport cpiercedVector
cimport csurfUnstructured

cdef extern from "utils.hpp":
    void initDoubleDataOnMesh(const csurfUnstructured.SurfUnstructured * mesh, cpiercedVector.PiercedVector[double,long] * data)