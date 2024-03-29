from mmic.components import TacticComponent
from mmelemental.models import Molecule
from mmelemental.util.units import convert
from cmselemental.util.decorators import classproperty
from mmic_translator.models import InputTrans, OutputTrans
from ..mmic_mda import units, __version__, __package__, _supported_versions
from typing import List, Tuple, Optional, Set
import MDAnalysis

__all__ = ["MolToMdaComponent", "MdaToMolComponent"]


provenance_stamp = {
    "creator": __package__,
    "version": __version__,
    "routine": __name__,
}


class MolToMdaComponent(TacticComponent):
    """A component for converting Molecule to MDAnalysis molecule object."""

    @classproperty
    def input(cls):
        return InputTrans

    @classproperty
    def output(cls):
        return OutputTrans

    @classproperty
    def version(cls) -> str:
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
        inputs: InputTrans,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, OutputTrans]:

        if isinstance(inputs, dict):
            inputs = self.input(**inputs)

        mmol = inputs.schema_object
        natoms = len(mmol.symbols)

        if hasattr(mmol, "substructs"):
            nres = len(mmol.substructs)
            resnames, resids = zip(*mmol.substructs)
            resids = [resid + 1 for resid in resids]
            # resnames = [resname.decode() for resname in resnames]
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
            forces=False,
        )

        mda_mol.add_TopologyAttr("type", mmol.symbols)
        mda_mol.add_TopologyAttr("mass", mmol.masses)
        mda_mol.atoms.masses = convert(
            mda_mol.atoms.masses, mmol.masses_units, units["mass"]
        )

        if mmol.atom_labels is not None:
            mda_mol.add_TopologyAttr("name", mmol.atom_labels)

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

        if mmol.connectivity is not None:
            bonds = [(bond[0], bond[1]) for bond in mmol.connectivity]
            mda_mol.add_TopologyAttr("bonds", bonds)
            # How to load bond order?

        success = True
        return success, OutputTrans(
            proc_input=inputs,
            data_object=mda_mol,
            data_units=units,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
            success=success,
            provenance=provenance_stamp,
        )


class MdaToMolComponent(TacticComponent):
    """A component for converting MDAnalysis molecule to Molecule object."""

    @classproperty
    def input(cls):
        return InputTrans

    @classproperty
    def output(cls):
        return OutputTrans

    @classproperty
    def version(cls) -> str:
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
        inputs: "MdaMol",
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, OutputTrans]:

        if isinstance(inputs, dict):
            inputs = self.input(**inputs)

        # get all properties + more from Universe?
        uni = inputs.data_object
        geo = getattr(uni.atoms, "positions", None)
        vel = getattr(uni.atoms, "velocities", None)
        symbs = getattr(uni.atoms, "types", None)
        names = getattr(uni.atoms, "names", None)
        masses = getattr(uni.atoms, "masses", None)

        # If bond order is none, set it to 1.
        if hasattr(uni.atoms, "bonds"):
            connectivity = [
                (bond.indices[0], bond.indices[1], bond.order or 1)
                for bond in uni.atoms.bonds
            ]
        else:
            connectivity = None

        residues = [(atom.resname, atom.resnum - 1) for atom in uni.atoms]

        input_dict = {
            "symbols": symbs,
            "geometry": geo,
            "geometry_units": units["length"],
            "velocities": vel,
            "velocities_units": units["speed"],
            "substructs": residues,
            "masses": masses,
            "masses_units": units["mass"],
            "atom_labels": names,
            "extras": inputs.keywords.get("extras"),
        }

        if connectivity:
            input_dict["connectivity"] = connectivity

        success = True
        return success, OutputTrans(
            proc_input=inputs,
            schema_object=Molecule(**input_dict),
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
            success=success,
            provenance=provenance_stamp,
        )
