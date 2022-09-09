import pytest
from click.testing import CliRunner
from coastlines.continental import continental_cli

def test_generate_rasters_cli():
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
            "[10000, 5000, 1000]",
            "--baseline_year",
            "2020",
        ]
    )
    assert result.exit_code == 0
    assert result.output == ''