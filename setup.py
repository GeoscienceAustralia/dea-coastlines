#!/usr/bin/env python3

import os
from setuptools import find_packages, setup

# Where are we?
IS_SANDBOX = "sandbox" in os.getenv("JUPYTER_IMAGE", default="")

# What packages are required for this module to be executed?
REQUIRED = [
    "aiohttp",
    "affine",
    "botocore",
    "click",
    "datacube",
    "dea_tools",
    "Fiona",
    "geopandas",
    "matplotlib",
    "numpy",
    "odc-geo",
    "odc_ui",
    "pandas",
    "pygeos",
    "pyproj",
    "pytest",
    "pyTMD",
    "python_geohash",
    "pytz",
    "PyYAML",
    "rasterio",
    "scikit_image",
    "scikit_learn",
    "scipy",
    "setuptools",
    "Shapely",
    "tqdm",
    "xarray",
]

# Package metadata
NAME = "dea_coastlines"
DESCRIPTION = "Tools for running Digital Earth Australia Coastlines"
URL = "https://github.com/GeoscienceAustralia/dea-coastlines"
EMAIL = "Robbi.BishopTaylor@ga.gov.au"
AUTHOR = "Robbi Bishop-Taylor"
REQUIRES_PYTHON = ">=3.8.0"

# Setup kwargs
setup_kwargs = {
    "name": NAME,
    "description": DESCRIPTION,
    "long_description": DESCRIPTION,
    "long_description_content_type": "text/markdown",
    "author": AUTHOR,
    "author_email": EMAIL,
    "python_requires": REQUIRES_PYTHON,
    "url": URL,
    "install_requires": REQUIRED if not IS_SANDBOX else [],
    "packages": find_packages(),
    "include_package_data": True,
    "license": "Apache License 2.0",
    "entry_points": {
        "console_scripts": [
            "deacoastlines-raster = coastlines.raster:generate_rasters_cli",
            "deacoastlines-vector = coastlines.vector:generate_vectors_cli",
            "deacoastlines-continental = coastlines.continental:continental_cli",
        ]
    },
}

setup(**setup_kwargs)
