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


def pytest_generate_tests(metafunc):
    if "guess_bonds" in metafunc.fixturenames:
        metafunc.parametrize("guess_bonds", ["False", "True"])
    if "md_traj_file" in metafunc.fixturenames:
        metafunc.parametrize("md_traj_file", [mm_data.trajs["traj.trr"]])


def test_mmic_mda_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_mda" in sys.modules


def test_mda_to_traj(md_traj_file, **kwargs):
    traj = mda.Universe(md_traj_file)
    inputs = {
        "data_object": traj,
        "schema_version": 1,
        "schema_name": "mmel_input",
        "keywords": kwargs,
    }
    mm_traj = mmic_mda.components.MdaToTrajComponent.compute(inputs)

    return mm_traj


@pytest.mark.parametrize(
    "top_file,traj_file",
    [
        (mm_data.mols["1dzl_fixed.gro"], mm_data.trajs["traj.trr"]),
    ],
)
def test_mda_to_top_traj(guess_bonds, top_file, traj_file, **kwargs):
    traj = mda.Universe(top_file, traj_file, guess_bonds=guess_bonds)
    inputs = {
        "data_object": traj,
        "schema_version": 1,
        "schema_name": "mmel_input",
        "keywords": kwargs,
    }
    mm_traj = mmic_mda.components.MdaToTrajComponent.compute(inputs)

    return mm_traj


@pytest.mark.parametrize(
    "traj_file",
    mm_data.trajs.values(),
)
def test_traj_to_mda(traj_file, **kwargs):
    mm_traj = mm.models.Trajectory.from_file(traj_file, all_frames=True)

    if mm_traj.ndim != 3:
        return None

    inputs = {
        "schema_object": mm_traj,
        "schema_version": 1,
        "schema_name": "mmel_input",
        "keywords": kwargs,
    }
    mda_traj = mmic_mda.components.TrajToMdaComponent.compute(inputs)

    return mda_traj


@pytest.mark.parametrize(
    "top_file",
    [
        mm_data.mols["1dzl_fixed.gro"],
    ],
)
def test_io_methods(guess_bonds, top_file):
    mda_traj = mmic_mda.models.MdaTraj.from_file(top_file, guess_bonds=guess_bonds)
    assert isinstance(mda_traj.data, mda_traj.dtype)

    mm_traj = mda_traj.to_schema()
    assert isinstance(mm_traj, mm.models.collect.Trajectory)
