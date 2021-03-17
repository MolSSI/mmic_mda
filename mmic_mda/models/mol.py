from typing import Dict, Any, Optional
from mmic_translator.models import ToolkitModel
from mmelemental.models.molecule import Molecule
import MDAnalysis

# MDAnalysis converter components
from mmic_mda.components.mol_component import MolToMdaComponent
from mmic_mda.components.mol_component import MdaToMolComponent

__all__ = ["MdaMol"]


class MdaMol(ToolkitModel):
    """ A model for MDAnalysis.Universe storing an MM Molecule. """

    @property
    def dtype(self):
        """ Returns the fundamental molecule object type. """
        return MDAnalysis.Universe

    @classmethod
    def isvalid(cls, data):
        """ Makes sure the Universe object stores atoms. """
        if hasattr(data, "atoms"):
            if len(data.atoms):
                return data

        raise ValueError("MDAnalysis molecule object does not contain any atoms!")

    @classmethod
    def from_file(
        cls, filename: str = None, top_filename: str = None, dtype: str = None, **kwargs
    ) -> "MdaMol":
        """
        Constructs an MdaMol object from file(s).

        Parameters
        ----------
        filename : str, optional
            The atomic positions filename to read
        top_filename: str, optional
            The topology filename to read
        dtype: str, optional
            The type of file to interpret. If unset, MDAnalysis attempts to discover dtype from the file extension.
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        MdaMol
            A constructed MdaMol class.
        """
        if dtype:
            kwargs["format"] = dtype

        if filename and top_filename:
            mol = MDAnalysis.Universe(top_filename, filename, **kwargs)
        elif filename:
            mol = MDAnalysis.Universe(filename, **kwargs)
        elif top_filename:
            mol = MDAnalysis.Universe(top_filename, **kwargs)
        else:
            raise TypeError(
                "You must supply at least one of the following: filename or top_filename."
            )

        return cls(data=mol)

    @classmethod
    def from_schema(
        cls, data: Molecule, version: Optional[str] = None, **kwargs: Dict[str, Any]
    ) -> "MdaMol":
        """
        Constructs an MdaMol object from an MMSchema Molecule object.
        Parameters
        ----------
        data: Molecule
            Data to construct Molecule from.
        version: str, optional
            Schema version e.g. 1.0.1
        **kwargs
            Additional kwargs to pass to the constructors. kwargs take precedence over data.
        Returns
        -------
        MdaMol
            A constructed MdaMol class.
        """
        inputs = {"schema_object": data, "schema_version": version}
        out = MolToMdaComponent.compute(inputs)
        return cls(data=out.tk_object, units=out.tk_units)

    def to_file(self, filename: str, dtype: str = None, **kwargs):
        """Writes the molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        """
        if dtype:
            kwargs["file_format"] = dtype
        self.data.atoms.write(filename, **kwargs)

    def to_schema(self, version: Optional[str] = None, **kwargs) -> Molecule:
        """Converts the molecule to MMSchema Molecule.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with e.g. 1.0.1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = {"tk_object": self.data, "schema_version": version, "kwargs": kwargs}
        out = MdaToMolComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
