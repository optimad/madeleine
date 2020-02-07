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
from gems.parallel.mda.parallel_jacobi import ParallelMDAJacobi
from gems_wrapper import ToySphereDiscipline

COMM = get_world_comm()

class TestGEMSWrapper(unittest.TestCase):

    def test_basic(self):
        MPIManager().clear(1)

        mesh_file = "../examples/data/unitSphere5.stl"
        neutral_mesh_file = "../examples/data/unitSphere4.stl"

        toy1 = ToySphereDiscipline("Sphere1", ["Forces"], ["Pressure"], mesh_file,
                                   neutral_mesh_file, sphere_radius=1.0,
                                   n_cpus=2, kernel=1)

        toy2 = ToySphereDiscipline("Sphere2", ["Pressure"], ["Forces"], mesh_file,
                                   neutral_mesh_file, sphere_radius=1.0,
                                   n_cpus=3, kernel=2)

        disciplines = [toy1, toy2]
        acceleration = None
        # acceleration = ParallelMDAJacobi.SECANT_ACCELERATION
        # acceleration = ParallelMDAJacobi.M2D_ACCELERATION
        mda = ParallelMDAJacobi(
            disciplines=disciplines,
            name="mda",
            max_mda_iter=200,
            tolerance=1e-9,
            acceleration=acceleration,
        )

        default_inputs = mda.default_inputs
        inputs = None
        if mda.execution_context.is_rank_on_mpi_group():
            r = ones(1)
            size_forces = default_inputs['Forces'].shape[0]
            forces = ones(size_forces)*0.01
            pressure = ones(size_forces)*0.01
            inputs = {'r': r, 'Forces': forces, 'Pressure': pressure}

        out = mda.execute(inputs)

        if mda.execution_context.is_rank_on_mpi_group():
            print out

if __name__ == "__main__":
    unittest.main()
