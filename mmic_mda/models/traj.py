from pydantic import Field, validator
from typing import Dict, Any, Optional
from mmelemental.models.base import ToolkitModel
from mmelemental.models.molecule import Molecule
from mmelemental.models.trajectory import Trajectory
from mmelemental.util.decorators import require


__all__ = ["MdaTraj"]


class MdaTraj(ToolkitModel):
    """ A model for MDAnalysis.Universe storing an MM trajectory. """

    @property
    @require("MDAnalysis")
    def dtype(self):
        """ Returns the fundamental molecule object type. """
        from MDAnalysis import Universe

        return Universe

    @validator("data")
    def valid_traj(cls, data):
        """ Makes sure the Universe object stores atoms and a trajectory. """
        if not hasattr(data, "atoms"):
            raise ValueError("MDAnalysis object does not contain any atoms!")
        elif hasattr(data, "trajectory"):
            return data

        raise ValueError("MDAnalysis object does not contain any trajectories!")

    @classmethod
    @require("MDAnalysis")
    def from_file(
        cls, filename: str, top_filename: str = None, dtype: str = None, **kwargs
    ) -> "MdaTraj":
        """
        Constructs an instance of MdaTraj object from file(s).

        Parameters
        ----------
        filename : str, optional
            The trajectory filename to read
        top_filename: str, optional
            The topology filename to read
        dtype: str, optional
            The type of file to interpret. If unset, MDAnalysis attempts to discover dtype from the file extension.
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        Trajectory
            A constructed Trajectory class.
        """
        import MDAnalysis

        if dtype:
            kwargs["format"] = dtype

        if filename and top_filename:
            mol = MDAnalysis.Universe(top_filename, filename, **kwargs)
        else:
            mol = MDAnalysis.Universe(filename, **kwargs)

        return cls(data=mol)

    @classmethod
    def from_schema(
        cls,
        data: Trajectory,
        version: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> "MdaTraj":
        """
        Constructs a MdaTraj object from an MMSchema Trajectory object.
        Parameters
        ----------
        data: Trajectory
            Data to construct Molecule from.
        version: str, optional
            Schema version e.g. 1.0.1
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        MdaTraj
            A constructed MDAnalysis.Universe object.
        """
        from mmic_mda.components.traj_component import TrajToMdaComponent

        return TrajToMdaComponent.compute(data)

    def to_file(self, filename: str, **kwargs):
        """Writes the trajectory to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
        """
        self.data.atoms.write(filename, **kwargs)

    def to_schema(self, version: Optional[str] = None, **kwargs) -> Trajectory:
        """Converts the MdaTraj to MMSchema Trajectory.
        Parameters
        ----------
        version: str, optional
            Schema version e.g. 1.0.1
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        from mmic_mda.components.traj_component import MdaToTrajComponent

        return MdaToTrajComponent.compute(self)
