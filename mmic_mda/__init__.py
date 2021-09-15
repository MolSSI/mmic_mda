"""
mmic_mda
A short description of the project.
"""

# Add imports here
from . import components, models
from .mmic_mda import *

_classes_map = {"Molecule": models.MdaMol, "Trajectory": models.MdaTraj}
