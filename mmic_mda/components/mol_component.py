from mmelemental.models.molecule import Molecule
from typing import List, Tuple, Optional
from mmelemental.util.units import convert
import MDAnalysis

from mmic_translator import TransComponent
from mmic_translator.models.io import (
    TransInput,
    TransOutput,
)

from ..mmic_mda import units

__all__ = ["MolToMdaComponent", "MdaToMolComponent"]


class MolToMdaComponent(TransComponent):
    """ A component for converting Molecule to MDAnalysis molecule object. """

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

        mmol = inputs.schema_object
        natoms = len(mmol.symbols)

        if mmol.residues:
            residues = list(fast_set(mmol.residues))
            nres = len(residues)
            resnames, _ = zip(*residues)
            _, resids = zip(*mmol.residues)
            resids = [i - 1 for i in resids]
        else:
            nres = 1
            resnames, resids = "UNK", [1]

        # Must account for segments as well
        segindices = None

        mda_mol = MDAnalysis.Universe.empty(
            natoms,
            n_residues=nres,
            atom_resindex=resids,
            residue_segindex=segindices,
            trajectory=True,
            velocities=True,
            forces=True,
        )

        mda_mol.add_TopologyAttr("type", mmol.symbols)
        mda_mol.add_TopologyAttr("mass", mmol.masses)
        mda_mol.atoms.masses = convert(
            mda_mol.atoms.masses, mmol.masses_units, units["mass"]
        )

        if mmol.atom_labels is not None:
            mda_mol.add_TopologyAttr("name", mmol.atom_labels)

        if mmol.residues is not None:
            mda_mol.add_TopologyAttr("resname", resnames)

        # mda_mol.add_TopologyAttr('segid', ['SOL'])

        if mmol.geometry is not None:
            mda_mol.atoms.positions = mmol.geometry.reshape(natoms, 3)
            mda_mol.atoms.positions = convert(
                mda_mol.atoms.positions, mmol.geometry_units, units["length"]
            )

        if mmol.velocities is not None:
            mda_mol.atoms.velocities = mmol.velocities.reshape(natoms, 3)
            mda_mol.atoms.velocities = convert(
                mda_mol.atoms.velocities,
                mmol.velocities_units,
                units["speed"],
            )

        if mmol.forces is not None:
            mda_mol.atoms.forces = mmol.forces.reshape(natoms, 3)
            mda_mol.atoms.forces = convert(
                mda_mol.atoms.forces, mmol.forces_units, units["force"]
            )

        if mmol.connectivity:
            bonds = [(bond[0], bond[1]) for bond in mmol.connectivity]
            mda_mol.add_TopologyAttr("bonds", bonds)
            # How to load bond order?

        return True, TransOutput(
            proc_input=inputs, data_object=mda_mol, data_units=units
        )


class MdaToMolComponent(TransComponent):
    """ A component for converting MDAnalysis molecule to Molecule object. """

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

    def execute(
        self,
        inputs: "MdaMol",
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        # get all properties + more from Universe?
        uni = inputs.data_object
        geo = TransComponent.get(uni.atoms, "positions")
        vel = TransComponent.get(uni.atoms, "velocities")
        forces = TransComponent.get(uni.atoms, "forces")
        symbs = TransComponent.get(uni.atoms, "types")
        names = TransComponent.get(uni.atoms, "names")
        if names is not None:
            names = names.tolist()  # must be list for MMSchema
        masses = TransComponent.get(uni.atoms, "masses")

        # If bond order is none, set it to 1.
        if hasattr(uni.atoms, "bonds"):
            connectivity = [
                (bond.indices[0], bond.indices[1], bond.order or 1)
                for bond in uni.atoms.bonds
            ]
        else:
            connectivity = None

        residues = [(atom.resname, atom.resnum) for atom in uni.atoms]

        input_dict = {
            "symbols": symbs,
            "geometry": geo,
            "geometry_units": units["length"],
            "velocities": vel,
            "velocities_units": units["speed"],
            "forces": forces,
            "forces_units": units["force"],
            "residues": residues,
            "connectivity": connectivity,
            "masses": masses,
            "masses_units": units["mass"],
            "atom_labels": names,
        }

        return True, TransOutput(
            proc_input=inputs, schema_object=Molecule(**input_dict)
        )

    def get_version(self) -> str:
        """Finds program, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
        raise NotImplementedError


def fast_set(seq: List) -> List:
    """ Removes duplicate entries in a list while preserving the order. """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
