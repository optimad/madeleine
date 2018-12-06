#def Py_initDoubleDataOnMesh(Py_SurfUnstructured mesh, Py_PiercedVector data):
#    initDoubleDataOnMesh(mesh.thisptr,data.thisptr)

cimport cpiercedVector
cimport csurfUnstructured
cimport cutils

cdef class Dummy:
    
    def __cinit__(self):
        print("Dummy created")
        
    def __dealloc__(self):
        print("Dummy destroyed")
        
    cdef test(self,const csurfUnstructured.SurfUnstructured * mesh, cpiercedVector.PiercedVector[double,long] * data):
        cutils.initDoubleDataOnMesh(mesh,data)
    