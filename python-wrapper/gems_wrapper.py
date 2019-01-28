# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
# Copyright (c) 2019 IRT-AESE
# All rights reserved.
#
# Contributors:
#    INITIAL AUTHORS - initial API and implementation and/or
#                      initial documentation
#        :author:  Francois Gallard
#    OTHER AUTHORS   - MACROSCOPIC CHANGES


from gems.core.discipline import MDODiscipline
import coupling as coupling_interface
from numpy import array


class ToySphereDiscipline(MDODiscipline):

    def __init__(self, name, inputs, outputs, mesh_file,
                 neutral_mesh_file, sphere_radius=1.0):

        self.inputs = inputs
        self.outputs = outputs
        self.mesh_file = mesh_file
        self.neutral_mesh_file = neutral_mesh_file
        self.sphere_radius = sphere_radius

        super(ToySphereDiscipline, self).__init__(name)
        self.input_grammar.initialize_from_data_names(inputs)
        self.output_grammar.initialize_from_data_names(outputs)

        self.mesh_coupling = coupling_interface.Py_MeshCoupling(inputs,
                                                                outputs)
        self.mesh_coupling.initialize(mesh_file, neutral_mesh_file,
                                      sphere_radius)

        mesh = self.mesh_coupling.getNeutralMesh()
        self.neutral_input = coupling_interface.Py_PiercedVector()
        self.neutral_output = coupling_interface.Py_PiercedVector()

        coupling_interface.Py_initDoubleDataOnMesh(mesh, self.neutral_input)

    def _run(self):
        input_vector = self.local_data[self.inputs[0]]
        # TOdo : update self.neutral_input from input_data
        self.neutral_input[:] = input_vector
        self.sphere_coupling.compute(self.neutral_input, self.neutral_output)

        self.local_data[self.outputs[0]] = array(self.neutral_output)

    def close(self):
        self.sphere_coupling.close()
