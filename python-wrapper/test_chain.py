# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
# Copyright (c) 2019 IRT-AESE
# All rights reserved.
#
# Contributors:
#    INITIAL AUTHORS - initial API and implementation and/or
#                      initial documentation
#        :author:  Francois Gallard
#    OTHER AUTHORS   - MACROSCOPIC CHANGES

import unittest
from numpy import ones

from gems.parallel.core.mpi_manager import MPIManager, get_world_comm
from gems.parallel.core.parallel_chain import MDOParallelMPIChain
from gems_wrapper import ToySphereDiscipline

COMM = get_world_comm()

class TestGEMSWrapper(unittest.TestCase):

    def test_basic(self):
        MPIManager().clear(2)

        mesh_file = "../examples/data/unitSphere5.stl"
        neutral_mesh_file = "../examples/data/unitSphere4.stl"

        toy1 = ToySphereDiscipline("Sphere1", ["Forces"], ["Pressure"], mesh_file,
                                   neutral_mesh_file, sphere_radius=1.0,
                                   n_cpus=2)

        toy2 = ToySphereDiscipline("Sphere1", ["Pressure"], ["Forces"], mesh_file,
                                   neutral_mesh_file, sphere_radius=1.0,
                                   n_cpus=2)

        chain = MDOParallelMPIChain([toy1, toy2])
        default_inputs = chain.default_inputs
        inputs = None
        if chain.execution_context.is_rank_on_mpi_group():
            r = ones(1)
            size_forces = default_inputs['Forces'].shape[0]
            forces = ones(size_forces)*10.
            inputs = {'r': r, 'Forces': forces}

        out = chain.execute(inputs)

        if chain.execution_context.is_rank_on_mpi_group():
            print out

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
