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
from numpy import zeros

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

        self.neutral_input = coupling_interface.Py_PiercedVector()
        
    def _run(self):
        input_vector = self.local_data[self.inputs[0]]
        
        mesh = self.mesh_coupling.getNeutralMesh()
        coupling_interface.Py_initDataOnMeshFromArray(mesh, self.neutral_input,input_vector)
        self.mesh_coupling.compute(self.neutral_input)

        output_vector = zeros(input_vector.shape[0])

        coupling_interface.Py_moveDataOnMeshToArray(mesh, self.neutral_input,output_vector)

        self.local_data[self.outputs[0]] = output_vector

    def close(self):
        self.mesh_coupling.close()
