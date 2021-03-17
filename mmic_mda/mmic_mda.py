"""
mmic_mda.py
A short description of the project.

Handles the primary functions
"""

__all__ = [
    "molread_ext_maps",
    "molwrite_ext_maps",
    "ffread_ext_maps",
    "ffwrite_ext_maps",
    "units",
]

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
