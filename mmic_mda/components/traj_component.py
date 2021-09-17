from mmelemental.models import Molecule, Trajectory
from cmselemental.util.decorators import classproperty
from mmic.components import TacticComponent
from typing import List, Tuple, Optional, Set
import numpy

from mmic_translator.models import (
    TransInput,
    TransOutput,
)
import MDAnalysis as mda
from ..mmic_mda import units, _supported_versions

__all__ = ["TrajToMdaComponent", "MdaToTrajComponent"]


class TrajToMdaComponent(TacticComponent):
    """A component for converting Trajectory to MDAnalysis Universe object."""

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    @classmethod
    def get_version(cls) -> str:
        """Returns distutils-style version string.

        Examples
        --------
        The string ">1.0, !=1.5.1, <2.0" implies any version after 1.0 and before 2.0
        is compatible, except 1.5.1

        Returns
        -------
        str
            Return a dist-utils valid version string.

        """
        return _supported_versions

    @classproperty
    def strategy_comps(cls) -> Set[str]:
        """Returns the strategy component(s) this (tactic) component belongs to.
        Returns
        -------
        Set[str]
        """
        return {"mmic_translator"}

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
        nframes = inputs.schema_object.nframes
        natoms = inputs.schema_object.natoms
        ndim = inputs.schema_object.ndim  # does MDA support 2d or 1D trajs?
        kwargs = inputs.keywords

        if inputs.schema_object.geometry is not None:
            coords = inputs.schema_object.geometry
        if inputs.schema_object.velocities is not None:
            velocities = inputs.schema_object.velocities
            kwargs.setdefault("velocities", velocities)
        if inputs.schema_object.forces is not None:
            forces = inputs.schema_object.forces
            kwargs.setdefault("forces", forces)

        if mm_traj.top is not None:
            top = mm_traj.top

            if isinstance(top, (list, tuple)):
                raise NotImplementedError("MDAnalysis supports only static topologies!")

            data_object = mda.Universe.empty(
                natoms,
                # n_residues=nres,
                # atom_resindex=resids,
                # residue_segindex=segindices,
                trajectory=True,
                velocities="velocities" in kwargs,
                forces="forces" in kwargs,
            )

            data_object.add_TopologyAttr("type", top.symbols)
            data_object.add_TopologyAttr("mass", top.masses)
            data_object.atoms.masses = convert(
                data_object.atoms.masses, top.masses_units, units["mass"]
            )

            if top.atom_labels is not None:
                data_object.add_TopologyAttr("name", top.atom_labels)
            if top.connectivity is not None:
                bonds = [
                    (bond[0], bond[1]) for bond in top.connectivity
                ]  # How to load bond order?
                data_object.add_TopologyAttr("bonds", bonds)

            for frame in range(nframes):
                data_object.trajectory.ts.dt = inputs.schema_object.timestep
                data_object.trajectory.ts.positions = coords[:, :, frame]
                if "velocities" in kwargs:
                    data_object.trajectory.ts.velocities = velocities[:, :, frame]
                if "forces" in kwargs:
                    data_object.trajectory.ts.forces = forces[:, :, frame]
        else:
            coords = coords.reshape(nframes, natoms, ndim)
            kwargs.update(velocities=velocities.reshape(nframes, natoms, ndim))
            kwargs.update(forces=forces.reshape(nframes, natoms, ndim))
            kwargs.setdefault("dt", inputs.schema_object.timestep)

            data_object = mda.Universe(coords, **kwargs)

        success = True
        return success, TransOutput(
            proc_input=inputs,
            data_object=data_object,
            data_units=units,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
            success=success,
        )


class MdaToTrajComponent(TacticComponent):
    """A component for converting MDAnalysis Universe to Molecule object."""

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    @classmethod
    def get_version(cls) -> str:
        """Finds program, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
        raise NotImplementedError

    @classproperty
    def strategy_comps(cls) -> Set[str]:
        """Returns the strategy component(s) this (tactic) component belongs to.
        Returns
        -------
        Set[str]
        """
        return {"mmic_translator"}

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
                "schema_name": inputs.schema_name,
                "keywords": inputs.keywords,
            }
            out = MdaToMolComponent.compute(inputs_mol)
            if getattr(out, "schema_version", None):
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
            schema_name=inputs.schema_name,
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
