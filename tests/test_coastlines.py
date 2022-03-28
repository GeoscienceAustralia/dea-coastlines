# Hardly any tests, but a start!
import geopandas as gpd

from coastlines.utils import STYLES_FILE


def test_styles_file():
    assert STYLES_FILE.exists()

    styles = gpd.read_file(STYLES_FILE)

    assert styles.shape[0] == 3
