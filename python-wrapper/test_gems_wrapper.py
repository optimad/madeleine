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
from gems_wrapper import ToySphereDiscipline

COMM = get_world_comm()
SIZE = COMM.size


class TestGEMSWrapper(unittest.TestCase):

    def test_basic(self):
        MPIManager().clear(0)
        mesh_file = "../examples/data/unitSphere5.stl"
        neutral_mesh_file = "../examples/data/unitSphere4.stl"
        toy1 = ToySphereDiscipline("Sphere1", ["Forces"], ["Pressure"],
                                   mesh_file, neutral_mesh_file,
                                   sphere_radius=1.0, n_cpus=SIZE)

        neutral_mesh_size = toy1.mesh_coupling.getNeutralMeshSize()
        toy1.execute({"Forces": ones(neutral_mesh_size)})
        toy1._compute_jacobian()
        toy1.close()


if __name__ == "__main__":
    unittest.main()
