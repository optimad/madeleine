import coupling
import numpy as np
import mpi4py
from mpi4py import MPI
inputs=["Forces"]
outputs=["Pressure"]
comm = MPI.COMM_WORLD
mc=coupling.Py_MeshCoupling(inputs,outputs,"toy1",comm)
disciplineFile = "../examples/data/unitSphere5.stl"
neutralFile = "../examples/data/unitSphere4.stl"
cellIndicesPerRank = coupling.Py_computeMeshFilePartitioning(neutralFile,comm)
print(cellIndicesPerRank)
mc.initialize(disciplineFile,neutralFile,2.0,cellIndicesPerRank)
nofNeutralElements=mc.getNeutralMeshSize()
arr=np.ones(nofNeutralElements)
mc.compute(arr)
print(arr)
mc.close()
