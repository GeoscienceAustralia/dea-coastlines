import pytest
from click.testing import CliRunner
from coastlines.continental import continental_cli

def test_generate_rasters_cli():
    runner = CliRunner()
    result = runner.invoke(
        continental_cli,
        [
            "--vector_version",
            "{{inputs.parameters.result-version}}",
            "--shorelines",
            "True",
            "--hotspots",
            "True",
            "--baseline_year",
            " {{inputs.parameters.baseline-year}}"
        ]
    )
    assert result.exit_code == 0
    assert result.output == ''