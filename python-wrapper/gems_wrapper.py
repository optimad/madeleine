# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
# Copyright (c) 2019 IRT-AESE
# All rights reserved.
#
# Contributors:
#    INITIAL AUTHORS - initial API and implementation and/or
#                      initial documentation
#        :author:  Francois Gallard
#    OTHER AUTHORS   - MACROSCOPIC CHANGES
from __future__ import print_function
from numpy import ones
import numpy as np
from petsc4py import PETSc

from gems.core.discipline import MDODiscipline
from gems_mpi_plugins.api import create_execution_context, create_user_partition
from gems_mpi_plugins.core.mpi_manager import MPIManager
import coupling

PETSC_DETERMINE = PETSc.DETERMINE


class ToySphereDiscipline(MDODiscipline):
    def __init__(
        self,
        name,
        inputs,
        outputs,
        mesh_file,
        neutral_mesh_file,
        sphere_radius=1.0,
        sphere_neutral_radius=1.0,
        sphere_thickness=0.001,
        is_inner_sphere=True,
        source_intensity=10.0,
        source_direction=None,
        thermalDiffusivityCoefficient=0.5,
        emissivity=0.5,
        infinityTemperature=400.0,
        n_cpus=1,
        kernel=1,
    ):

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

        self.input_grammar.initialize_from_data_names([inputs[0], "r"])
        self.output_grammar.initialize_from_data_names([outputs[0]])

        # Creation of the execution context and register
        self.execution_context = create_execution_context(self, n_cpus)

        comm = self.execution_context.comm

        cell_indices_per_rank = 0
        local_size = 0
        all_sizes = None
        if self.execution_context.is_alive():
            cell_indices_per_rank = coupling.Py_computeMeshFilePartitioning(
                neutral_mesh_file.encode(), comm
            )
            local_indices = [i for i in cell_indices_per_rank if i == comm.rank]
            local_size = len(local_indices)
            all_sizes = comm.allgather(local_size)

        # Broadcast all_sizes to all the ranks of comm_world
        comm_world = MPIManager().main_comm
        root = self.execution_context.comm_world_ranks[0]
        all_sizes = comm_world.bcast(all_sizes, root=root)

        obj_variables = [0] * n_cpus
        obj_variables[0] = 1

        variables_sizes = {
            inputs[0]: all_sizes,
            outputs[0]: all_sizes,
            "r": obj_variables,
        }

        partition = create_user_partition(self.execution_context, variables_sizes)

        # Init the discipline only if the ranks belong to the discipline
        # subcommunicator
        partition.initialize_default_inputs(1.0)

        if self.execution_context.is_alive():
            if comm.rank == 0:
                self.default_inputs["r"] = ones(1) * sphere_radius

        if self.execution_context.is_alive():
            inputs_encode = [inp.encode() for inp in inputs]
            outputs_encode = [out.encode() for out in outputs]
            self.mesh_coupling = coupling.Py_MeshCoupling(
                inputs_encode, outputs_encode, name.encode(), comm
            )

            self.mesh_coupling.initialize(
                mesh_file.encode(),
                neutral_mesh_file.encode(),
                sphere_radius,
                sphere_neutral_radius,
                sphere_thickness,
                is_inner_sphere,
                source_intensity,
                source_direction,
                thermalDiffusivityCoefficient,
                emissivity,
                infinityTemperature,
                cell_indices_per_rank,
                self.kernel,
            )

    def _run(self):
        input_vector = self.local_data[self.inputs[0]]
        r = self.local_data["r"]

        comm = self.execution_context.comm
        r = comm.bcast(r, root=0)

        # self.mesh_coupling.set_radius(r)
        self.mesh_coupling.compute(input_vector, r, r)
        output = {self.outputs[0]: input_vector}
        self.store_local_data(**output)

    def close(self):
        self.mesh_coupling.close()

    def _compute_jacobian(self, inputs=None, outputs=None):

        inputs = list(self.get_input_data_names())
        outputs = list(self.get_output_data_names())

        # Compute global number of neutral cell
        nofNeutralLocalCells = self.mesh_coupling.getNeutralMeshSize()
        comm = self.execution_context.comm
        nofNeutralGlobalCells = comm.allreduce(nofNeutralLocalCells)

        # Initialize Jacobian matrices
        self._init_jacobian(with_zeros=True)
        self.jac[outputs[0]] = {}
        matInputs = PETSc.Mat().create(comm=comm)
        matInputs.setSizes(
            (
                (nofNeutralLocalCells, PETSC_DETERMINE),
                (nofNeutralLocalCells, PETSC_DETERMINE),
            )
        )
        matInputs.setType("dense")
        matInputs.setUp()

        for i in range(nofNeutralLocalCells):
            rowId, colIds, values = self.mesh_coupling.extractOutputInputJacobianRow(i)
            rowIds = [rowId]
            # CAVEAT: petsc4py accepts only int32 by default. bitpit indices are long integers. Cast is possible but very large meshes are not feasible
            matInputs.setValues(rowIds, colIds, values, addv=1)

        matInputs.assemblyBegin()
        matInputs.assemblyEnd()
        self.jac[outputs[0]][inputs[0]] = matInputs
        viewer = PETSc.Viewer().createASCII("CouplingJac.dat", mode=1, comm=comm)
        viewer.view(self.jac[outputs[0]][inputs[0]])

        matControl = PETSc.Mat().create(comm=comm)
        matControl.setSizes(((nofNeutralLocalCells, PETSC_DETERMINE), (1, 1)))
        matControl.setType("dense")
        matControl.setUp()

        for i in range(nofNeutralLocalCells):
            rowId, colIds, values = self.mesh_coupling.extractOutputControlJacobianRow(
                i
            )
            rowIds = [rowId]
            colIds = [1]
            # CAVEAT: petsc4py accepts only int32 by default. bitpit indices are long integers. Cast is possible but very large meshes are not feasible
            print("Jac r", i, values)
            matControl.setValues(rowIds, colIds, values, addv=1)

        matControl.assemblyBegin()
        matControl.assemblyEnd()
        self.jac[outputs[0]]["r"] = matControl

        viewer = PETSc.Viewer().createASCII("RadiusJac.dat", mode=1, comm=comm)
        viewer.view(self.jac[outputs[0]]["r"])
