import coupling
inputs=["Forces"]
outputs=["Pressure"]
mc=coupling.Py_MeshCoupling(inputs,outputs)
disciplineFile = "../examples/data/unitSphere1.stl"
neutralFile = "../examples/data/unitSphere2.stl"
mc.initialize(disciplineFile,neutralFile,2.0)
nD=coupling.Py_PiercedVector()
#mesh=coupling.Py_SurfUnstructured()
mesh=mc.getNeutralMesh()
coupling.Py_initDoubleDataOnMesh(mesh,nD)
mc.compute(nD)
mc.close()
