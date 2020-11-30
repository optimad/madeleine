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
from numpy import ones, full
from copy import deepcopy

from gems_mpi_plugins.core.mpi_manager import MPIManager
from gems_wrapper import ToySphereDiscipline

COMM = MPIManager().main_comm
SIZE = COMM.size


class TestGEMSWrapper(unittest.TestCase):
    def test_outer(self):
        MPIManager().clear(0)
        mesh_file = "../examples/data/unitSphere5.stl"
        neutral_mesh_file = "../examples/data/unitSphere4.stl"
        toy1 = ToySphereDiscipline(
            "Sphere1",
            ["T_out"],
            ["T_in"],
            mesh_file,
            neutral_mesh_file,
            sphere_radius=1.0,
            n_cpus=SIZE,
            is_inner_sphere=False,
            source_intensity=1.0,
            source_direction=None,
            thermalDiffusivityCoefficient=0.25,
            emissivity=0.25,
            infinityTemperature=350.0,
        )

        # neutral_mesh_size = toy1.mesh_coupling.getNeutralMeshSize()

        # t_out = {"T_out": full(neutral_mesh_size, 300.0)}
        # res = toy1.execute(deepcopy(t_out))
        # res.update(t_out)
        # print(res)

        # toy1._compute_jacobian()
        # print(toy1.jac)
        # toy1.close()


if __name__ == "__main__":
    unittest.main()
