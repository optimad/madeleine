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
thermalDiffusivityCoefficient = 0.001
emissivity = 0.001
infinityTemperature = 274.0
mc.initialize(disciplineFile, neutralFile, 2.0, 2.0, 1.0, True, 1.0, npSourceDirection, thermalDiffusivityCoefficient,
              emissivity, infinityTemperature, cellIndicesPerRank, 0)
nofNeutralElements=mc.getNeutralMeshSize()
arr=np.ones(nofNeutralElements)
arr=arr * 274.0
mc.compute(arr,2.0,2.0)
cellGlobalId,retIds,retValues = mc.extractOutputInputJacobianRow(0)
#print(retValues.shape)
print(arr)
mc.close()
