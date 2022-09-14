import pytest
from click.testing import CliRunner
from coastlines.raster import generate_rasters_cli
from coastlines.continental import continental_cli
from coastlines.vector import generate_vectors_cli


def test_generate_rasters_cli():
    runner = CliRunner()
    result = runner.invoke(
        generate_rasters_cli,
        [
            "--config_path",
            "configs/dea_coastlines_config_testing.yaml",
            "--study_area",
            "1098",
            "--raster_version",
            "testing",
            "--start_year",
            "2015",
            "--end_year",
            "2020",
        ]
    )
    assert result.exit_code == 0
    # assert result.output == ''


@pytest.mark.depends(on=['test_generate_rasters_cli'])
def test_generate_vector_cli():
    runner = CliRunner()
    result = runner.invoke(
        generate_vectors_cli,
        [
            "--config_path",
            "configs/dea_coastlines_config_testing.yaml",
            "--study_area",
            "1098",
            "--raster_version",
            "testing",
            "--start_year",
            "2015",
            "--end_year",
            "2020",
            "--baseline_year",
            "2020",
        ]
    )
    assert result.exit_code == 0
    # assert result.output == ''



@pytest.mark.depends(on=['test_generate_vector_cli'])
def test_generate_continental_cli():
    runner = CliRunner()
    result = runner.invoke(
        continental_cli,
        [
            "--vector_version",
            "testing",
            "--shorelines",
            "True",
            "--ratesofchange",
            "True",
            "--hotspots",
            "True",
            "--hotspots_radius",
            "100",
            "--baseline_year",
            "2020",
        ]
    )
    # assert result.output == '' # for debugging
    assert result.exit_code == 0
