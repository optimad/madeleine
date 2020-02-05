include "piercedVector.pyx"
include "surfUnstructured.pyx"
from libc.stdint cimport uintptr_t
import numpy as np
cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "couplingUtils.hpp" namespace "coupling":
    cdef void initDoubleDataOnMesh(const SurfUnstructured * mesh, PiercedVector[double,long] * data)
    cdef void initDataOnMeshFromArray(const SurfUnstructured * mesh, PiercedVector[double,long] * data, double * arr, size_t arrSize)
    cdef void moveDataOnMeshToArray(const SurfUnstructured * mesh, PiercedVector[double,long] * data, double * arr, size_t arrSize)
    cdef vector[int] computeMeshFilePartitioning(const string meshFile,MPI_Comm);

    
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
    
def Py_computeMeshFilePartitioning(meshFile,MPI.Comm comm):
    cdef MPI_Comm mpi_comm
    cdef vector[int] cellIndicesPerRank
    mpi_comm = comm.ob_mpi    
    #filename = meshFile
    cellIndicesPerRank = computeMeshFilePartitioning(<string&> meshFile,mpi_comm)
    #print(cellIndicesPerRank)
    #cdef int[::1] view = <int[:cellIndicesPerRank.size()]> cellIndicesPerRank.data()
    ret = np.asarray(cellIndicesPerRank)
    return ret