"""
Unit and regression test for the mmic_mda package.
"""

# Import package, test suite, and other packages as needed
import mmic_mda
import pytest
import sys
import os
import MDAnalysis as mda
import mmelemental as mm
import mm_data


top_file = (
    mm_data.mols["alanine.gro"],
    mm_data.mols["dialanine.pdb"],
    mm_data.mols["dialanine.gro"],
)


def pytest_generate_tests(metafunc):
    if "guess_bonds" in metafunc.fixturenames:
        metafunc.parametrize("guess_bonds", ["False", "False", "True"])


def test_mmic_mda_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_mda" in sys.modules


def test_mda_to_mol(guess_bonds, **kwargs):
    for file in top_file:
        uni = mda.Universe(file, guess_bonds=guess_bonds)
        inputs = {
            "data_object": uni,
            "keywords": kwargs,
            "schema_version": 1,
            "schema_name": "mmel_input",
        }
        mm_mol = mmic_mda.components.MdaToMolComponent.compute(inputs)


def test_mol_to_mda(guess_bonds):
    for file in top_file:
        mm_mol = mm.models.Molecule.from_file(file)
        inputs = {
            "schema_object": mm_mol,
            "schema_version": 1,
            "schema_name": "mmel_input",
        }
        mmic_mda.components.MolToMdaComponent.compute(inputs)


def test_io_methods(guess_bonds):
    for file in top_file:
        mda_mol = mmic_mda.models.MdaMol.from_file(file, guess_bonds=guess_bonds)
        assert isinstance(mda_mol.data, mda_mol.dtype)

        mm_mol = mda_mol.to_schema()
        assert isinstance(mm_mol, mm.models.Molecule)
