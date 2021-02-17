from mmelemental.components.trans.template_component import TransComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.molecule import Molecule
from mmic_mda.models import MdaMol
from typing import Dict, Any, List, Tuple, Optional
from mmelemental.util.decorators import require
from mmelemental.util.units import convert

__all__ = ["MolToMdaComponent", "MdaToMolComponent"]


class MolToMdaComponent(TransComponent):
    """ A component for converting Molecule to MDAnalysis molecule object. """

    @classmethod
    def input(cls):
        return Molecule

    @classmethod
    def output(cls):
        return MdaMol

    @require("MDAnalysis")
    def execute(
        self,
        inputs: Molecule,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, MdaMol]:

        import MDAnalysis
        import mmic_mda

        natoms = len(inputs.masses)

        if inputs.residues:
            residues = list(fast_set(inputs.residues))
            nres = len(residues)
            resnames, _ = zip(*residues)
            _, resids = zip(*inputs.residues)
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

        mda_mol.add_TopologyAttr("type", inputs.symbols)
        mda_mol.add_TopologyAttr("mass", inputs.masses)
        mda_mol.atoms.masses = convert(
            mda_mol.atoms.masses, inputs.masses_units, mmic_mda.units["mass"]
        )

        if inputs.names is not None:
            mda_mol.add_TopologyAttr("name", inputs.names)

        if inputs.residues is not None:
            mda_mol.add_TopologyAttr("resname", resnames)

        # mda_mol.add_TopologyAttr('segid', ['SOL'])

        if inputs.geometry is not None:
            mda_mol.atoms.positions = inputs.geometry.reshape(natoms, 3)
            mda_mol.atoms.positions = convert(
                mda_mol.atoms.positions, inputs.geometry_units, mmic_mda.units["length"]
            )

        if inputs.velocities is not None:
            mda_mol.atoms.velocities = inputs.velocities.reshape(natoms, 3)
            mda_mol.atoms.velocities = convert(
                mda_mol.atoms.velocities,
                inputs.velocities_units,
                mmic_mda.units["speed"],
            )

        if inputs.forces is not None:
            mda_mol.atoms.forces = inputs.forces.reshape(natoms, 3)
            mda_mol.atoms.forces = convert(
                mda_mol.atoms.forces, inputs.forces_units, mmic_mda.units["force"]
            )

        if inputs.connectivity:
            bonds = [(bond[0], bond[1]) for bond in inputs.connectivity]
            mda_mol.add_TopologyAttr("bonds", bonds)
            # How to load bond order?

        return True, MdaMol(data=mda_mol, units=mmic_mda.units)


class MdaToMolComponent(TransComponent):
    """ A component for converting MDAnalysis molecule to Molecule object. """

    @classmethod
    def input(cls):
        return MdaMol

    @classmethod
    def output(cls):
        return Molecule

    def execute(
        self,
        inputs: MdaMol,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Molecule]:

        import mmic_mda

        # get all properties + more from Universe?
        uni = inputs.data
        geo = TransComponent.get(uni.atoms, "positions")
        vel = TransComponent.get(uni.atoms, "velocities")
        forces = TransComponent.get(uni.atoms, "forces")
        symbs = TransComponent.get(uni.atoms, "types")
        names = TransComponent.get(uni.atoms, "names")
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
            "geometry_units": mmic_mda.units["length"],
            "velocities": vel,
            "velocities_units": mmic_mda.units["speed"],
            "forces": forces,
            "forces_units": mmic_mda.units["force"],
            "residues": residues,
            "connectivity": connectivity,
            "masses": masses,
            "masses_units": mmic_mda.units["mass"],
            "names": names,
        }

        return True, Molecule(**input_dict)


def fast_set(seq: List) -> List:
    """ Removes duplicate entries in a list while preserving the order. """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
