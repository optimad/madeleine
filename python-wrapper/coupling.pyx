from Cython.Shadow import NULL
include "utils.pyx"

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "coupling.hpp" namespace "coupling":
    #ctypedef vector[string] svector
    
    cdef cppclass MeshCoupling:
        MeshCoupling() except +
        MeshCoupling(vector[string],vector[string]) except +
        void initialize(string,string,double)
        void compute(PiercedVector[double,long]*)
        void close()
        const SurfUnstructured * getNeutralMesh()
        size_t getNeutralMeshSize()
        
cdef class Py_MeshCoupling:
    cdef MeshCoupling* thisptr
    
    def __cinit__(self,*args):
        thisptr = NULL
        cdef vector[string] inputs
        cdef vector[string] outputs
        
        n_args =len(args)
        
        if(n_args == 0):
            self.thisptr = new MeshCoupling()
        elif(n_args == 2):
            for str in args[0]:
                inputs.push_back(str)
            for str in args[1]:
                outputs.push_back(str)
            self.thisptr = new MeshCoupling(inputs,outputs)
        else:
            print("Py_MeshCoupling, wrong number of arguments!")
            
    def __dealloc__(self):
        del self.thisptr
    
    def initialize(self,unitDisciplineMeshFile,unitNeutralMeshFile,sphereRadius):
        self.thisptr.initialize(<string&> unitDisciplineMeshFile, <string&> unitNeutralMeshFile,<double> sphereRadius)
    
    def compute(self,Py_PiercedVector neutralData):
        self.thisptr.compute(<PiercedVector[double,long]*><void*> neutralData.thisptr)

    def close(self):
        self.thisptr.close()

    def getNeutralMesh(self):
        cdef uintptr_t int_ptr = <uintptr_t>(<MeshCoupling*><void*>self.thisptr)[0].getNeutralMesh()
        py_mesh = Py_SharedSurfUnstructured(int_ptr)
        return py_mesh

    def getNeutralMeshSize(self):
        return (<MeshCoupling*><void*>self.thisptr)[0].getNeutralMeshSize()