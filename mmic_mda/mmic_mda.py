"""
mmic_mda.py
A short description of the project.

Handles the primary functions
"""

__all__ = [
    "molread_ext_maps",
    "molwrite_ext_maps",
    "trajread_ext_maps",
    "trajwrite_ext_maps",
    "ffread_ext_maps",
    "ffwrite_ext_maps",
    "units",
    "__version__",
]

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

# Trajectory read/write file extensions
trajread_ext_maps = {".dcd": "dcd", ".trr": "trr"}
trajwrite_ext_maps = {".dcd": "dcd", ".trr": "trr"}

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
