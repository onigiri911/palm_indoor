"""Run the Berlin test cases"""
import os

import pytest

from palm_csd.create_driver import create_driver
from palm_csd.csd_config import _reset_all_config_counters
from palm_csd.netcdf_data import remove_existing_file
from tests.tools import ncdf_equal


@pytest.fixture
def config_counters():
    """Reset config class counters to allow independent test runs."""
    yield
    _reset_all_config_counters()


@pytest.fixture
def generate_configuration_wrong_range_tree_trunk_diameter():
    """Generate a configuration file with the tree trunk diameter wrong range input file."""

    with open("tests/99_full_application/berlin_tiergarten.yml") as f:
        configuration = f.read()
    configuration = configuration.replace(
        "Berlin_trees_trunk_clean_15m_DLR.nc", "Berlin_trees_trunk_clean_15m_DLR_wrong_range.nc"
    )
    configuration = configuration.replace(
        "Berlin_trees_trunk_clean_3m_DLR.nc", "Berlin_trees_trunk_clean_3m_DLR_wrong_range.nc"
    )
    configuration = configuration.replace(
        "file_out: berlin_tiergarten", "file_out: berlin_tiergarten_wrong_range_tree_trunk_diameter"
    )
    with open(
        "tests/99_full_application/berlin_tiergarten_wrong_range_tree_trunk_diameter.yml", "w"
    ) as f:
        f.write(configuration)
    yield
    os.remove("tests/99_full_application/berlin_tiergarten_wrong_range_tree_trunk_diameter.yml")


@pytest.mark.usefixtures("config_counters")
def test_complete_run():
    """Run the Berlin test case and compare with correct output"""

    create_driver("tests/99_full_application/berlin_tiergarten.yml")

    assert ncdf_equal(
        "tests/99_full_application/output/berlin_tiergarten_root",
        "tests/99_full_application/berlin_tiergarten_root",
    ), "Root driver does not comply with reference"

    assert ncdf_equal(
        "tests/99_full_application/output/berlin_tiergarten_N02",
        "tests/99_full_application/berlin_tiergarten_N02",
    ), "Child driver does not comply with reference"

    os.remove("tests/99_full_application/berlin_tiergarten_root")
    os.remove("tests/99_full_application/berlin_tiergarten_N02")


@pytest.mark.usefixtures(
    "config_counters", "generate_configuration_wrong_range_tree_trunk_diameter"
)
def test_wrong_range_tree_trunk_diamter():
    """Run the Berlin test case with a wrong range of the tree trunk diameter.
    This should raise an error.
    """
    with pytest.raises(ValueError):
        create_driver(
            "tests/99_full_application/berlin_tiergarten_wrong_range_tree_trunk_diameter.yml"
        )

    remove_existing_file(
        "tests/99_full_application/berlin_tiergarten_wrong_range_tree_trunk_diameter_root"
    )
    remove_existing_file(
        "tests/99_full_application/berlin_tiergarten_wrong_range_tree_trunk_diameter_N02"
    )
