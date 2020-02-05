include "surfUnstructured.pyx"
from libc.stdint cimport uintptr_t
import numpy as np
cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "couplingUtils.hpp" namespace "coupling":
    cdef vector[int] computeMeshFilePartitioning(const string meshFile,MPI_Comm);

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
