import coupling
import numpy as np
import mpi4py
from mpi4py import MPI
inputs=["Temperature"] #outer sphere
outputs=["Flux"] #outer sphere
comm = MPI.COMM_WORLD
mc=coupling.Py_MeshCoupling(inputs, outputs, "InnerSphere", comm)
disciplineFile = "../examples/data/unitSphere5.stl"
neutralFile = "../examples/data/unitSphere4.stl"
cellIndicesPerRank = coupling.Py_computeMeshFilePartitioning(neutralFile,comm)
print(cellIndicesPerRank)
sourceDirection = [1.0, 0.0, 0.0]
npSourceDirection = np.asarray(sourceDirection)
mc.initialize(disciplineFile, neutralFile, 2.0, 2.0, 1.0, False, 1.0, npSourceDirection, cellIndicesPerRank,0)
nofNeutralElements=mc.getNeutralMeshSize()
arr=np.ones(nofNeutralElements)
arr=arr * 274.0
mc.compute(arr,2.0,2.0)
print(arr)
mc.close()
