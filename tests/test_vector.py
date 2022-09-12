import pytest
from click.testing import CliRunner
from coastlines.vector import generate_vectors_cli

@pytest.mark.depends(on=['test_generate_rasters_cli'])
def test_generate_vector_cli():
    runner = CliRunner()
    result = runner.invoke(
        generate_vectors_cli,
        [
            "--config_path",
            "configs/dea_coastlines_config.yaml",
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
    assert result.output == ''