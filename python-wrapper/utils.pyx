include "piercedVector.pyx"
include "surfUnstructured.pyx"
from libc.stdint cimport uintptr_t
import numpy as np

cdef extern from "couplingUtils.hpp" namespace "coupling":
    cdef void initDoubleDataOnMesh(const SurfUnstructured * mesh, PiercedVector[double,long] * data)
    cdef void initDataOnMeshFromArray(const SurfUnstructured * mesh, PiercedVector[double,long] * data, double * arr, size_t arrSize)
    cdef void moveDataOnMeshToArray(const SurfUnstructured * mesh, PiercedVector[double,long] * data, double * arr, size_t arrSize)
    
def Py_initDoubleDataOnMesh(Py_SharedSurfUnstructured mesh, Py_PiercedVector data):
    initDoubleDataOnMesh(<SurfUnstructured*><void*>mesh.thisptr,<PiercedVector[double,long]*><void*>data.thisptr)

def Py_initDataOnMeshFromArray(Py_SharedSurfUnstructured mesh, Py_PiercedVector data, arr):
    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr)
    cdef double[::1] arr_memview = arr
    initDataOnMeshFromArray(<SurfUnstructured*><void*>mesh.thisptr,<PiercedVector[double,long]*><void*>data.thisptr,&arr_memview[0],arr.shape[0])
     
def Py_moveDataOnMeshToArray(Py_SharedSurfUnstructured mesh, Py_PiercedVector data, arr):
    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr)
    cdef double[::1] arr_memview = arr
    moveDataOnMeshToArray(<SurfUnstructured*><void*>mesh.thisptr,<PiercedVector[double,long]*><void*>data.thisptr,&arr_memview[0],arr.shape[0])
    