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
from copy import deepcopy
from numpy import ones, full
from numpy.random import random

from gems_mpi_plugins.core.mpi_manager import MPIManager
from gems_wrapper import ToySphereDiscipline

COMM = MPIManager().main_comm
SIZE = COMM.size


def test_inner():
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
        n_cpus=1,
        is_inner_sphere=True,
        source_intensity=1.0,
        source_direction=None,
        thermalDiffusivityCoefficient=0.0,
        emissivity=1e-6,
        infinityTemperature=350.0,
    )

    toy1.add_differentiated_inputs(["T_out"])
    toy1.add_differentiated_outputs(["T_in"])

    t_out = {}
    if toy1.execution_context.is_rank_on_mpi_group():
        neutral_mesh_size = toy1.mesh_coupling.getNeutralMeshSize()
        t_array = full(neutral_mesh_size, 300.0) + random(neutral_mesh_size) * 10.0
        t_out = {"T_out": t_array}

    res = toy1.execute(deepcopy(t_out))
    t_out.update(res)
    toy1.linearize()


if __name__ == "__main__":
    test_inner()
