from mmelemental.models.collect import Trajectory, Frame
from typing import List, Tuple, Optional
import MDAnalysis

from mmic_translator import TransComponent
from mmic_translator.models.io import (
    TransInput,
    TransOutput,
)

from ..mmic_mda import units

__all__ = ["TrajToMdaComponent", "MdaToTrajComponent"]


class TrajToMdaComponent(TransComponent):
    """ A component for converting Molecule to MDAnalysis molecule object. """

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

        return True, TransOutput(tk_object=uni, tk_units=units)


class MdaToTrajComponent(TransComponent):
    """ A component for converting MDAnalysis molecule to Molecule object. """

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

        mol = None
        uni = inputs.tk_object

        if hasattr(uni.atoms, "names"):
            from mmic_mda.components.mol_component import MdaToMolComponent

            inputs = {
                "tk_object": uni,
                "schema_version": inputs.tk_version,
                "kwargs": inputs.kwargs,
            }
            out = MdaToMolComponent.compute(inputs)
            if TransComponent.get(out, "schema_version"):
                assert inputs.schema_version == out.schema_version

        frames = [
            Frame(
                geometry=frame.positions if frame.has_positions else None,
                geometry_units="A",
                velocities=frame.velocities if frame.has_velocities else None,
                velocities_units="A/ps",
                forces=frame.forces if frame.has_forces else None,
                forces_units="kJ/(mol*A)",
                timestep=frame.dt,
                timestep_units="ps",
            )
            for frame in uni.trajectory
        ]
        # By using frames we are assuming the topology is constant. Is this always true in MDAnalysis?
        return True, TransOutput(schema_object=Trajectory(mol=mol, frames=frames))
