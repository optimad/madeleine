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

from gems.core.discipline import MDODiscipline
from gems.parallel.api import create_execution_context, create_user_partition
import coupling 


class ToySphereDiscipline(MDODiscipline):

    def __init__(self, name, inputs, outputs, mesh_file,
                 neutral_mesh_file, sphere_radius=1.0, n_cpus=1):

        self.inputs = inputs
        self.outputs = outputs
        self.mesh_file = mesh_file
        self.neutral_mesh_file = neutral_mesh_file
        self.sphere_radius = sphere_radius

        super(ToySphereDiscipline, self).__init__(name)

        self.input_grammar.initialize_from_data_names([inputs[0], 'r'])
        self.output_grammar.initialize_from_data_names([outputs[0]])

        # Creation of the execution context and register
        self.execution_context = create_execution_context(self, n_cpus)

        comm = self.execution_context.comm
        cell_indices_per_rank = coupling.Py_computeMeshFilePartitioning(neutral_mesh_file,
                                                                     comm)
        local_size = cell_indices_per_rank.shape[0]
        all_sizes = comm.allgather(local_size)

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

        if comm.rank == 0:
            self.default_inputs['r'] = ones(1)*sphere_radius

        self.mesh_coupling = coupling.Py_MeshCoupling(inputs, outputs, name,
                                                      comm)

        self.mesh_coupling.initialize(mesh_file, neutral_mesh_file,
                                      sphere_radius, cell_indices_per_rank)

    def _run(self):
        input_vector = self.local_data[self.inputs[0]]
        r = self.local_data['r']

        comm = self.execution_context.comm
        r = comm.bcast(r, root=0)

        # self.mesh_coupling.set_radius(r)
        self.mesh_coupling.compute(input_vector)

        self.local_data[self.outputs[0]] = input_vector

    def close(self):
        self.mesh_coupling.close()
