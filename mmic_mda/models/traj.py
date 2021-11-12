from typing import Dict, Any, Optional
from mmic_translator.models import ToolkitModel
from mmelemental.models.collect import Trajectory
from cmselemental.util.decorators import classproperty
import MDAnalysis

# MDAnalysis converter components
from mmic_mda.components.traj_component import TrajToMdaComponent
from mmic_mda.components.traj_component import MdaToTrajComponent

__all__ = ["MdaTraj"]


class MdaTraj(ToolkitModel):
    """A model for MDAnalysis.Universe storing an MM trajectory."""

    @classproperty
    def engine(cls):
        return "MDAnalysis", MDAnalysis.__version__

    @classproperty
    def dtype(cls):
        """Returns the fundamental trajectory object type."""
        return MDAnalysis.Universe

    @classmethod
    def isvalid(cls, data):
        """Makes sure the Universe object stores atoms and a trajectory."""
        if not hasattr(data, "atoms"):
            raise ValueError("MDAnalysis object does not contain any atoms!")
        elif hasattr(data, "trajectory"):
            return data

        raise ValueError("MDAnalysis object does not contain any trajectories!")

    @classmethod
    def from_file(
        cls,
        filename: str,
        top_filename: str = None,
        dtype: str = None,
        *,
        all_frames=False,
        **kwargs: Optional[Dict[str, Any]],
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
        all_frames: bool, optional
            Reads all frames in memory.
        **kwargs: Dict[str, Any], optional
            Any additional keywords to pass to the constructor.
        Returns
        -------
        Trajectory
            A constructed Trajectory class.
        """
        if dtype:
            kwargs["format"] = dtype

        if filename and top_filename:
            mol = MDAnalysis.Universe(
                top_filename, filename, in_memory=all_frames, **kwargs
            )
        else:
            mol = MDAnalysis.Universe(filename, in_memory=all_frames, **kwargs)

        return cls(data=mol)

    @classmethod
    def from_schema(
        cls,
        data: Trajectory,
        version: Optional[int] = 0,
        **kwargs: Optional[Dict[str, Any]],
    ) -> "MdaTraj":
        """
        Constructs a MdaTraj object from an MMSchema Trajectory object.
        Parameters
        ----------
        data: Trajectory
            Data to construct Trajectory from.
        version: int, optional
            Schema version e.g. 1. Overrides data.schema_version.
        **kwargs: Dict[str, Any], optional
            Additional kwargs to pass to the constructors.
        Returns
        -------
        MdaTraj
            A constructed MDAnalysis.Universe object.
        """
        inputs = {
            "schema_object": data,
            "schema_version": version or data.schema_version,
            "schema_name": data.schema_name,
        }
        return TrajToMdaComponent.compute(inputs)

    def to_file(self, filename: str, **kwargs):
        """Writes the trajectory to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
        """
        self.data.atoms.write(filename, **kwargs)

    def to_schema(self, version: Optional[int] = 0, **kwargs) -> Trajectory:
        """Converts the MdaTraj to MMSchema Trajectory.
        Parameters
        ----------
        version: int, optional
            Schema version e.g. 1
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = {
            "data_object": self.data,
            "schema_version": version,
            "schema_name": kwargs.pop("schema_name", Trajectory.default_schema_name),
            "keywords": kwargs,
        }
        out = MdaToTrajComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
