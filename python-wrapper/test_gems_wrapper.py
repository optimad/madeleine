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
from gems_wrapper import ToySphereDiscipline
from numpy import ones


class TestGEMSWrapper(unittest.TestCase):

    def test_basic(self):
        mesh_file = "../examples/data/unitSphere1.stl"
        neutral_mesh_file = "../examples/data/unitSphere2.stl"
        toy1 = ToySphereDiscipline("Sphere1", ["Forces"], ["Pressure"], mesh_file,
                                   neutral_mesh_file, sphere_radius=1.0)

        neutral_mesh_size = toy1.mesh_coupling.getNeutralMeshSize()
        toy1.execute({"Forces": ones(neutral_mesh_size)})

        toy1.close()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
