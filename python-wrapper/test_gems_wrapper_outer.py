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

from gems_mpi_plugins.core.mpi_manager import MPIManager
from gems_wrapper import ToySphereDiscipline

COMM = MPIManager().main_comm
SIZE = COMM.size


class TestGEMSWrapper(unittest.TestCase):
    # def test_inner(self):
    #     MPIManager().clear(0)
    #     mesh_file = "../examples/data/unitSphere5.stl"
    #     neutral_mesh_file = "../examples/data/unitSphere4.stl"
    #     toy1 = ToySphereDiscipline(
    #         "Sphere1",
    #         ["T_in"],
    #         ["T_out"],
    #         mesh_file,
    #         neutral_mesh_file,
    #         sphere_radius=1.0,
    #         n_cpus=SIZE,
    #         kernel=0,
    #     )

    #     neutral_mesh_size = toy1.mesh_coupling.getNeutralMeshSize()

    #     res = toy1.execute({"T_in": ones(neutral_mesh_size)})
    #     print(res)

    #     toy1._compute_jacobian()

    #     toy1.close()

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
            kernel=1,
        )

        neutral_mesh_size = toy1.mesh_coupling.getNeutralMeshSize()

        res = toy1.execute({"T_out": ones(neutral_mesh_size)})
        print(res)

        toy1._compute_jacobian()

        toy1.close()


if __name__ == "__main__":
    unittest.main()
