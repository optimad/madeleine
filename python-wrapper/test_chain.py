# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
# Copyright (c) 2019 IRT-AESE
# All rights reserved.
#
# Contributors:
#    INITIAL AUTHORS - initial API and implementation and/or
#                      initial documentation
#        :author:  Francois Gallard
#    OTHER AUTHORS   - MACROSCOPIC CHANGES
from __future__ import print_function, unicode_literals
import sys, petsc4py

petsc4py.init(sys.argv)
import unittest
from math import pi
from numpy import ones, full, mean
from numpy.random import random

from gems_mpi_plugins.core.mpi_manager import MPIManager
from gems_mpi_plugins.core.parallel_chain import MDOParallelMPIChain
from gems_wrapper import ToySphereDiscipline

COMM = MPIManager().main_comm
SIZE = MPIManager().main_comm.size


def test_basic():
    MPIManager().clear(1)

    mesh_file = "../examples/data/unitSphere5.stl"
    neutral_mesh_file = "../examples/data/unitSphere4.stl"

    toy_inner = ToySphereDiscipline(
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

    # neutral_mesh_size = toy_inner.mesh_coupling.getNeutralMeshSize()

    toy_outer = ToySphereDiscipline(
        "Sphere1",
        ["T_in"],
        ["T_out"],
        mesh_file,
        neutral_mesh_file,
        sphere_radius=1.0,
        n_cpus=1,
        is_inner_sphere=False,
        source_intensity=30.0,
        source_direction=None,
        thermalDiffusivityCoefficient=0.25,
        emissivity=0.1,
        infinityTemperature=400.0,
    )

    mesh_file = "../examples/data/unitSphere5.stl"
    neutral_mesh_file = "../examples/data/unitSphere4.stl"

    chain = MDOParallelMPIChain([toy_outer, toy_inner])
    default_inputs = chain.default_inputs
    inputs = None

    if chain.execution_context.is_rank_on_mpi_group():
        r = full(1, 10.0)
        neutral_size = default_inputs["T_in"].shape[0]
        t_array = full(neutral_size, 300.0) + random(neutral_size) * 10.0
        t_out = {"T_out": t_array}
        inputs = {"r": r}
        inputs.update(t_out)

    out = chain.execute(inputs)

    for i in range(50):
        out = chain.execute(out)

    if chain.execution_context.is_rank_on_mpi_group():
        target = 500.0
        obj = 0.5 * mean((out["T_in"] - target) ** 2)
        obj2 = 0.5 * mean((out["T_in"] - target) ** 2) * 4 * pi * r[0] ** 2
        print("Obj is", obj, obj2)

    # COMM.Barrier()
    # if toy_inner.execution_context.is_rank_on_mpi_group():
    #     toy_inner.close()
    # if toy_outer.execution_context.is_rank_on_mpi_group():
    #     toy_outer.close()


if __name__ == "__main__":
    test_basic()
