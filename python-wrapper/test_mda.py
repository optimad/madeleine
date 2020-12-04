# -*-mode: python; py-indent-offset: 4; tab-width: 8; coding:utf-8 -*-
# Copyright (c) 2019 IRT-AESE
# All rights reserved.
#
# Contributors:
#    INITIAL AUTHORS - initial API and implementation and/or
#                      initial documentation
#        :author:  Francois Gallard
#    OTHER AUTHORS   - MACROSCOPIC CHANGES
import gc

gc.disable()
import unittest
from numpy import ones, full, mean, array

from gems_mpi_plugins.api import create_execution_context, create_user_partition
from gems_mpi_plugins.core.mpi_manager import MPIManager
from gems_mpi_plugins.mda.parallel_jacobi import ParallelMDAJacobi
from gems_wrapper import ToySphereDiscipline
from gems.core.discipline import MDODiscipline
from gems_mpi_plugins.core.partitionning import UserPartioning
from gems_mpi_plugins.core.petsc import Mat
from gems_mpi_plugins.core.utils import even_divide

COMM = MPIManager().main_comm


class ObjectiveDiscipline(MDODiscipline):
    """SellarSystem discipline."""

    def __init__(
        self, n_cpus=1, variables_sizes=None, mesh_size=48, t_ref=350.0, **kwargs
    ):
        """Constructor.

        :param n_cpus: number of cpus to use for the discipline
        :param variables_sizes: dictionnary containing the partitionning of the
        variables
        :param coupling_size: size of the coupling variable vector
        """
        super(ObjectiveDiscipline, self).__init__(**kwargs)

        self.input_grammar.initialize_from_data_names(["r", "T_in"])
        self.output_grammar.initialize_from_data_names(["obj"])

        # Creation of the execution context and register
        self.execution_context = create_execution_context(self, n_cpus)
        partition = UserPartioning(self.execution_context)
        obj_variables = [0] * n_cpus
        obj_variables[0] = 1
        if variables_sizes is None:
            variables_sizes = {
                "r": obj_variables,
                "T_in": even_divide(mesh_size, n_cpus),
                "obj": obj_variables,
            }

        partition = create_user_partition(self.execution_context, variables_sizes)

        # Init the discipline only if the ranks belong to the discipline
        # subcommunicator
        partition.initialize_default_inputs(0.0)
        self.mesh_size = mesh_size
        self.t_ref = t_ref

    def _run(self):
        """Run the discipline."""
        r = self.get_inputs_by_name("r")
        t_in = self.get_inputs_by_name("T_in")

        if self.execution_context.comm.rank == 0:
            obj = array([0.5 * mean((t_in - self.t_ref) ** 2)])

        self.store_local_data(obj=obj)

    def _compute_jacobian(self, inputs=None, outputs=None):
        """Computes the jacobian.

        :param inputs: linearization should be performed with respect
            to inputs list. If None, linearization should
            be performed wrt all inputs (Default value = None)
        :param outputs: linearization should be performed on outputs list.
            If None, linearization should be performed
            on all outputs (Default value = None)
        """
        comm = self.execution_context.comm
        self._init_jacobian(with_zeros=True)

        r = self.get_inputs_by_name("r")
        t_in = self.get_inputs_by_name("T_in")

        r_0 = None

        if comm.rank == 0:
            r_0 = r[0]

        r_0 = comm.bcast(r_0, 0)

        mat = Mat(self.jac["obj"]["T_in"])
        for i in mat.mat.owner_range:
            dobj_dt_in = t_in[i] - self.t_ref
            mat.set_value(i, i, dobj_dt_in)
        mat.assemble()


class TestGEMSWrapper(unittest.TestCase):
    def test_basic(self):

        manager = MPIManager()
        manager.clear(1)
        manager.overlap_comms = True

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

        obj_discipline = ObjectiveDiscipline()

        disciplines = [toy_inner, toy_outer, obj_discipline]

        acceleration = None
        # acceleration = ParallelMDAJacobi.SECANT_ACCELERATION
        # acceleration = ParallelMDAJacobi.M2D_ACCELERATION

        COMM.Barrier()

        mda = ParallelMDAJacobi(
            disciplines=disciplines,
            name="mda",
            max_mda_iter=500,
            tolerance=1e-6,
            acceleration=acceleration,
        )

        mda.add_differentiated_inputs(["r"])
        mda.add_differentiated_outputs(["obj"])

        default_inputs = mda.default_inputs
        inputs = None

        if mda.execution_context.is_rank_on_mpi_group():
            r = ones(1)

            size_tin = default_inputs["T_in"].shape[0]
            t_in = full(size_tin, 300.0)
            t_out = full(size_tin, 400.0)
            inputs = {"r": r, "T_in": t_in, "T_out": t_out}

        out = mda.execute(inputs)

        out_jac = mda.linearize(inputs)

        if mda.execution_context.is_rank_on_mpi_group():
            print("Execution", out)

        COMM.Barrier()
        print("Its finished")


if __name__ == "__main__":
    unittest.main()
