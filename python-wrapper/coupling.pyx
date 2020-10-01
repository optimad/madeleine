from Cython.Shadow import NULL
include "utils.pyx"

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

cdef extern from "coupling.hpp" namespace "coupling":
    # ctypedef vector[string] svector

    cdef cppclass MeshCoupling:
        MeshCoupling(string, MPI_Comm) except +
        MeshCoupling(vector[string], vector[string], string, MPI_Comm) except +
        void initialize(string, string, double, double, double, bool, double, vector[double], double, double, double, vector[int], int) except +
        void compute(double * arr, size_t arrSize, double, double)
        void close()
        const SurfUnstructured * getNeutralMesh()
        size_t getNeutralMeshSize()
        void extractOutputInputJacobianRow(int,int,vector[int],vector[double])
        void extractOutputControlJacobianRow(int,int,vector[int],vector[double])
        long getNeutralFirstCellId()
        long getNeutralGlobalConsecutiveId(long)

cdef class Py_MeshCoupling:
    cdef MeshCoupling * thisptr

    def __cinit__(self, *args):
        thisptr = NULL
        cdef vector[string] inputs
        cdef vector[string] outputs
        cdef MPI_Comm mpi_comm

        n_args = len(args)

        if(n_args == 0):
            name = "dummy"
            mpi_comm = ( < MPI.Comm > MPI_COMM_WORLD).ob_mpi
            self.thisptr = new MeshCoupling(name, mpi_comm)
        elif(n_args == 4):
            for str in args[0]:
                inputs.push_back(str)
            for str in args[1]:
                outputs.push_back(str)
            name = args[2]
            mpi_comm = ( < MPI.Comm > args[3]).ob_mpi
            self.thisptr = new MeshCoupling(inputs, outputs, name, mpi_comm)
        else:
            print("Py_MeshCoupling, wrong number of arguments!")

    def __dealloc__(self):
        del self.thisptr

    def initialize(self, unitDisciplineMeshFile, unitNeutralMeshFile, sphereRadius, sphereNeutralRadius,
    sphereThickness, innerSphere, sourceIntensity, sourceDirection, thermalDiffusivityCoefficient,
     emissivity, infinityTemperature, cellIndicesPerRank, kernel):
#        cdef vector[double] sDirection
#        sDirection.push_back(sourceDirection[0])
#        sDirection.push_back(sourceDirection[1])
#        sDirection.push_back(sourceDirection[2])
        self.thisptr.initialize( < string&> unitDisciplineMeshFile,
                                < string & > unitNeutralMeshFile,
                                < double > sphereRadius,
                                < double > sphereNeutralRadius,
                                < double > sphereThickness,
                                < bool > innerSphere,
                                < double > sourceIntensity,
                                < const vector[double] & > sourceDirection,
                                < double > thermalDiffusivityCoefficient,
                                < double > emissivity,
                                < double > infinityTemperature,
                                < const vector[int] & > cellIndicesPerRank,
                                < int > kernel)

    def compute(self, neutralData,double newRadius,double otherRadius):
        if not neutralData.flags['C_CONTIGUOUS']:
            neutralData = np.ascontiguousarray(neutralData)
        cdef double[::1] arr_memview = neutralData
        self.thisptr.compute(& arr_memview[0], neutralData.shape[0], newRadius, otherRadius)

    def close(self):
        self.thisptr.close()

    def getNeutralMesh(self):
        cdef uintptr_t int_ptr = <uintptr_t > ( < MeshCoupling*> < void*>self.thisptr)[0].getNeutralMesh()
        py_mesh = Py_SharedSurfUnstructured(int_ptr)
        return py_mesh

    def getNeutralMeshSize(self):
        return ( < MeshCoupling*> < void*>self.thisptr)[0].getNeutralMeshSize()

    def extractOutputInputJacobianRow(self,int cellId):
        cdef vector[int] columnIds
        cdef int cellGlobalId = 0
        cdef vector[double] columnValues
        self.thisptr.extractOutputInputJacobianRow(cellId,<int&> cellGlobalId,<vector[int]&> columnIds, <vector[double]&> columnValues)
        retIds = np.asarray(columnIds,dtype=np.int32)
        retValues = np.asarray(columnValues,dtype=np.float64)
        return cellGlobalId,retIds,retValues

    def extractOutputControlJacobianRow(self,int cellId):
        cdef vector[int] columnIds
        cdef int cellGlobalId = 0
        cdef vector[double] columnValues
        self.thisptr.extractOutputControlJacobianRow(cellId,<int&> cellGlobalId,<vector[int]&> columnIds, <vector[double]&> columnValues)
        retIds = np.asarray(columnIds,dtype=np.int32)
        retValues = np.asarray(columnValues,dtype=np.float64)
        return cellGlobalId,retIds,retValues

    def getNeutralFirstCellId(self):
        return self.thisptr.getNeutralFirstCellId()
    
    def getNeutralGlobalConsecutiveId(self,long cellId):
        cdef globalIndex = self.thisptr.getNeutralGlobalConsecutiveId(cellId)
        return globalIndex