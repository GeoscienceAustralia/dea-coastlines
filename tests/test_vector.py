import pytest
from click.testing import CliRunner
from coastlines.vector import generate_vectors_cli

def test_generate_rasters_cli():
    runner = CliRunner()
    result = runner.invoke(
        generate_vectors_cli,
        [
            "--config_path",
            "{{inputs.parameters.config}}",
            "--study_area",
            "{{inputs.parameters.id}}",
            "--raster_version",
            "{{inputs.parameters.result-version}}",
            "--start_year",
            "{(inputs.parameters.start-year}}",
            "--end_year",
            "{{inputs.parameters.end-year}}",
            "{{inputs.parameters.overwrite}}"
        ]
    )
    assert result.exit_code == 0
    assert result.output == ''