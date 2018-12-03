from Cython.Shadow import NULL
include "piercedVector.pyx"
include "surfUnstructured.pyx"

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "coupling.hpp" namespace "coupling":
    #ctypedef vector[string] svector
    
    cdef cppclass MeshCoupling:
        MeshCoupling() except +
        MeshCoupling(vector[string],vector[string]) except +
        void initialize(string,string,double)
        void compute(PiercedVector[double,long]&)
        void close()
        
cdef class Py_MeshCoupling:
    cdef MeshCoupling* thisptr
    
    def __cinit__(self,*args):
        thisptr = NULL
        
        n_args =len(args)
        
        if(n_args == 0):
            self.thisptr = new MeshCoupling()
        elif(n_args == 1):
            self.thisptr = new MeshCoupling(args[0],args[1])
        else:
            print("Py_MeshCoupling, wrong numbr of arguments!")
            
    def __dealloc__(self):
        del self.thisptr
    
    def initialize(self,unitDisciplineMeshFile,unitNeutralMeshFile,sphereRadius):
        self.thisptr.initialize(<string&> unitDisciplineMeshFile, <string&> unitNeutralMeshFile,<double> sphereRadius)
    
    def compute(self,neutralData):
        self.thisptr.compute(<PiercedVector[double,long]&> neutralData)
    
    def close(self):
        self.thisptr.close()
    
    
    