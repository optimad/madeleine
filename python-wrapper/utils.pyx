include "piercedVector.pyx"

cdef extern from "utils.hpp" namespace "coupling":
    PiercedVector<double> initDoubleDataOnMesh(const SurfUnstructured * mesh)