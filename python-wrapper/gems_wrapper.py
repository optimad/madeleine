# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
# Copyright (c) 2019 IRT-AESE
# All rights reserved.
#
# Contributors:
#    INITIAL AUTHORS - initial API and implementation and/or
#                      initial documentation
#        :author:  Francois Gallard
#    OTHER AUTHORS   - MACROSCOPIC CHANGES
from numpy import ones
import numpy as np
from petsc4py import PETSc

from gems.core.discipline import MDODiscipline
from gems.parallel.api import create_execution_context, create_user_partition
from gems.parallel.core.mpi_manager import get_world_comm
import coupling


class ToySphereDiscipline(MDODiscipline):

    def __init__(self, name, inputs, outputs, mesh_file,
                 neutral_mesh_file, sphere_radius=1.0, sphere_neutral_radius=1.0,
                 sphere_thickness=0.001, is_inner_sphere=True,
                 source_intensity=1.0, source_direction=None,
                 n_cpus=1, kernel=1):

        if source_direction is None:
            source_direction = [1.0, 0.0, 0.0]
        np_source_direction = np.asarray(source_direction)
        self.inputs = inputs
        self.outputs = outputs
        self.mesh_file = mesh_file
        self.neutral_mesh_file = neutral_mesh_file
        self.sphere_radius = sphere_radius
        self.kernel = kernel

        super(ToySphereDiscipline, self).__init__(name)

        self.input_grammar.initialize_from_data_names([inputs[0], 'r'])
        self.output_grammar.initialize_from_data_names([outputs[0]])

        # Creation of the execution context and register
        self.execution_context = create_execution_context(self, n_cpus)

        comm = self.execution_context.comm

        cell_indices_per_rank = 0
        local_size = 0
        all_sizes = None
        if self.execution_context.is_alive():
            cell_indices_per_rank = coupling.Py_computeMeshFilePartitioning(neutral_mesh_file,
                                                                        comm)
            local_indices = [i for i in cell_indices_per_rank if i == comm.rank]
            local_size = len(local_indices)
            all_sizes = comm.allgather(local_size)

        # Broadcast all_sizes to all the ranks of comm_world
        comm_world = get_world_comm()
        root = self.execution_context.comm_world_ranks[0]
        all_sizes = comm_world.bcast(all_sizes, root=root)

        obj_variables = [0] * n_cpus
        obj_variables[0] = 1

        variables_sizes = {inputs[0]: all_sizes,
                           outputs[0]: all_sizes,
                           "r": obj_variables}

        partition = create_user_partition(self.execution_context,
                                          variables_sizes)

        # Init the discipline only if the ranks belong to the discipline
        # subcommunicator
        partition.initialize_default_inputs(1.)

        if self.execution_context.is_alive():
            if comm.rank == 0:
                self.default_inputs['r'] = ones(1)*sphere_radius

        if self.execution_context.is_alive():
            self.mesh_coupling = coupling.Py_MeshCoupling(inputs,
                                                          outputs,
                                                          name,
                                                          comm)

            self.mesh_coupling.initialize(mesh_file,
                                          neutral_mesh_file,
                                          sphere_radius,
                                          sphere_neutral_radius,
                                          sphere_thickness,
                                          is_inner_sphere,
                                          source_intensity,
                                          source_direction,
                                          cell_indices_per_rank,
                                          self.kernel)

    def _run(self):
        input_vector = self.local_data[self.inputs[0]]
        r = self.local_data['r']

        comm = self.execution_context.comm
        r = comm.bcast(r, root=0)

        # self.mesh_coupling.set_radius(r)
        self.mesh_coupling.compute(input_vector, r, r)
        output = {self.outputs[0]: input_vector}
        self.store_local_data(**output)

    def close(self):
        self.mesh_coupling.close()

    def _compute_jacobian(self, inputs=None, outputs=None):
        inputs = self.get_input_data_names()
        outputs = self.get_output_data_names()

        # Compute global number of neutral cell
        nofNeutralLocalCells = self.mesh_coupling.getNeutralMeshSize()
        comm = self.execution_context.comm
        nofNeutralGlobalCells = comm.allreduce(nofNeutralLocalCells)

        # Initialize Jacobian matrices
        self._init_jacobian(with_zeros=True)
        self.jac['Pressure'] = {}
        mat = PETSc.Mat().createAIJ(nofNeutralGlobalCells,nofNeutralGlobalCells,comm=comm)
        mat.setUp()

        for i in range(nofNeutralLocalCells):
            rowId,colIds,values = self.mesh_coupling.computeJacobianRow(i)
            rowIds = [rowId]
            #CAVEAT: petsc4py accepts only int32 by default. bitpit indices are long integers. Cast is possible but very large meshes are not feasible
            mat.setValues(rowIds,colIds,values, addv=1)

        mat.assemblyBegin()
        mat.assemblyEnd()
        self.jac['Pressure']['Forces'] = mat
        viewer = PETSc.Viewer().createASCII("ForcesJac.dat",mode=1,comm=comm)
        viewer.view(self.jac['Pressure']['Forces'])

        mat = PETSc.Mat().createAIJ(nofNeutralGlobalCells,1,comm=comm)
        mat.setUp()
        for i in range(beginId,beginId+nofNeutralLocalCells):
            mat.setValue(self.mesh_coupling.getNeutralGlobalConsecutiveId(i),0,i, addv=1)
        mat.assemblyBegin()
        mat.assemblyEnd()
        self.jac['Pressure']['r'] = mat

        viewer = PETSc.Viewer().createASCII("RadiusJac.dat",mode=1,comm=comm)
        viewer.view(self.jac['Pressure']['r'])
