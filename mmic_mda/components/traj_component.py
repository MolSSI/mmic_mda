from mmelemental.models import Molecule, Trajectory
from typing import List, Tuple, Optional
import numpy

from mmic_translator import TransComponent
from mmic_translator.models import (
    TransInput,
    TransOutput,
    schema_input_default,
    schema_output_default,
)

from ..mmic_mda import units

__all__ = ["TrajToMdaComponent", "MdaToTrajComponent"]


class TrajToMdaComponent(TransComponent):
    """A component for converting Trajectory to MDAnalysis Universe object."""

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    def execute(
        self,
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        mm_traj = inputs.schema_object
        if hasattr(mm_traj, "mol"):
            uni = mm_traj.mol.to_data("MDAnalysis")

        success = True
        return success, TransOutput(
            proc_input=inputs,
            data_object=uni,
            data_units=units,
            schema_version=inputs.schema_version,
            schema_name=schema_output_default,
            success=success,
        )


class MdaToTrajComponent(TransComponent):
    """A component for converting MDAnalysis Universe to Molecule object."""

    def execute(
        self,
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        uni = inputs.data_object

        if hasattr(uni.atoms, "names"):
            from mmic_mda.components.mol_component import MdaToMolComponent

            inputs_mol = {
                "data_object": uni,
                "schema_version": inputs.schema_version,
                "schema_name": schema_input_default,
                "keywords": inputs.keywords,
            }
            out = MdaToMolComponent.compute(inputs_mol)
            if TransComponent.get(out, "schema_version"):
                assert inputs.schema_version == out.schema_version

        geometry, velocities, forces = [], [], []

        for frame in uni.trajectory:
            geometry = (
                numpy.concatenate((geometry, frame.positions.flatten()))
                if frame.has_positions
                else geometry
            )
            velocities = (
                numpy.concatenate((velocities, frame.velocities.flatten()))
                if frame.has_velocities
                else velocities
            )
            forces = (
                numpy.concatenate((forces, frame.forces.flatten()))
                if frame.has_forces
                else forces
            )

        timestep = uni.trajectory[0].dt
        timestep_units = "ps"

        # By using frames we are assuming the topology is constant. Is this always true in MDAnalysis?
        success = True
        return success, TransOutput(
            proc_input=inputs,
            schema_object=Trajectory(
                natoms=uni.trajectory.n_atoms,  # mda assumes constant natoms?
                nframes=uni.trajectory.n_frames,
                timestep=timestep,
                timestep_units=timestep_units,
                geometry=geometry if len(geometry) > 0 else None,
                geometry_units="A",
                velocities=velocities if len(velocities) > 0 else None,
                velocities_units="A/ps",
                forces=forces if len(forces) > 0 else None,
                forces_units="kJ/(mol*A)",
            ),
            schema_version=inputs.schema_version,
            schema_name=schema_output_default,
            success=success,
        )

    def get_version(self) -> str:
        """Finds program, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
        raise NotImplementedError
