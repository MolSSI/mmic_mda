"""
mmic_mda
A short description of the project.
"""

# Add imports here
from . import components, models

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions


# Molecule (xyz) read/write file extensions
molread_ext_maps = {".gro": "gro", ".pdb": "pdb"}
molwrite_ext_maps = {".gro": "gro", ".pdb": "pdb"}

# Topology (connectivity and/or forcefield params) read/write extensions
ffread_ext_maps = {".psf": "psf"}
ffwrite_ext_maps = {".psf": "psf"}


units = {
    "length": "angstrom",
    "time": "ps",
    "energy": "kJ/mol",
    "charge": "e",
    "speed": "angstrom/ps",
    "force": "kJ/(mol*angstrom)",
    "mass": "amu",
    "angle": "degrees",
}

_classes_map = {"Molecule": models.MdaMol, "Trajectory": models.MdaTraj}
