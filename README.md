> **Note: This branch contains code that is currently being refactored from the [Digital Earth Africa Coastlines](https://github.com/digitalearthafrica/deafrica-coastlines/) repository. It is unlikely to work correctly.**

![Digital Earth Australia Coastlines](visualisation/images/DEACoastlines_header.gif)

# Digital Earth Australia Coastlines

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.rse.2021.112734-0e7fbf.svg)](https://doi.org/10.1016/j.rse.2021.112734)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Access:** [Explore DEA Coastlines on the interactive Digital Earth Australia Maps platform](https://maps.dea.ga.gov.au/#share=s-DEACoastlines&playStory=1)

**License:** The code in this repository is licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0). Digital Earth Australia data is licensed under the [Creative Commons by Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).

**Product description and metadata:** For the most up-to-date information about this product, visit the official [Geoscience Australia DEA Coastlines product description](https://cmi.ga.gov.au/data-products/dea/581/dea-coastlines-landsat)

**Contact:** For assistance with any of the Python code or Jupyter Notebooks in this repository, please post a [Github issue](https://github.com/GeoscienceAustralia/DEACoastLines/issues/new). For questions or more information about this product, email dea@ga.gov.au or sign up to the [Open Data Cube Slack](https://join.slack.com/t/opendatacube/shared_invite/zt-d6hu7l35-CGDhSxiSmTwacKNuXWFUkg) and post on the [`#dea-coastlines`](https://app.slack.com/client/T0L4V0TFT/C018X6J9HLY/details/) channel.

**To cite:** 
> Bishop-Taylor, R., Nanson, R., Sagar, S., Lymburner, L. (2021). Mapping Australia's dynamic coastline at mean sea level using three decades of Landsat imagery. _Remote Sensing of Environment_, 267, 112734. Available: https://doi.org/10.1016/j.rse.2021.112734

> Nanson, R., Bishop-Taylor, R., Sagar, S., Lymburner, L., (2022). Geomorphic insights into Australia's coastal change using a national dataset derived from the multi-decadal Landsat archive. _Estuarine, Coastal and Shelf Science_, 265, p.107712. Available: https://doi.org/10.1016/j.ecss.2021.107712

> Bishop-Taylor, R., Sagar, S., Lymburner, L., Alam, I., Sixsmith, J. (2019). Sub-pixel waterline extraction: characterising accuracy and sensitivity to indices and spectra. _Remote Sensing_, 11 (24):2984. Available: https://doi.org/10.3390/rs11242984

---

[**Digital Earth Australia Coastlines**](https://maps.dea.ga.gov.au/#share=s-DEACoastlines&playStory=1) is a continental dataset that includes annual shorelines and rates of coastal change along the entire Australian coastline from 1988 to the present. 

The product combines satellite data from Geoscience Australia's [Digital Earth Australia program](https://www.ga.gov.au/dea) with tidal modelling to map the typical location of the coastline at mean sea level for each year. The product enables trends of coastal erosion and growth to be examined annually at both a local and continental scale, and for patterns of coastal change to be mapped historically and updated regularly as data continues to be acquired. This allows current rates of coastal change to be compared with that observed in previous years or decades. 

The ability to map shoreline positions for each year provides valuable insights into whether changes to our coastline are the result of particular events or actions, or a process of more gradual change over time. This information can enable scientists, managers and policy makers to assess impacts from the range of drivers impacting our coastlines and potentially assist planning and forecasting for future scenarios. 

#### Applications
* Monitoring and mapping rates of coastal erosion along the Australian coastline 
* Prioritise and evaluate the impacts of local and regional coastal management based on historical coastline change 
* Modelling how coastlines respond to drivers of change, including extreme weather events, sea level rise or human development 
* Supporting geomorphological studies of how and why coastlines have changed across time 

---

## Table of contents
* [Repository code](#repository-code)
    * [Python modules](#python-modules)
    * [Jupyter notebooks](#jupyter-notebooks)
    * [Command-line interface](#command-line-interface)
* [DEA Coastlines dataset](#dea-coastlines-dataset)
    * [Data download](#data-download)
    * [Interactive map](#interactive-map)
    * [Loading DEA Coastlines data from the Web Feature Service (WFS) using Python](#loading-dea-coastlines-data-from-the-web-feature-service-wfs-using-python)
    * [Loading DEA Coastlines data from the Web Feature Service (WFS) using R](#loading-dea-coastlines-data-from-the-web-feature-service-wfs-using-r)
* [References](#references)

---

## Repository code
The code in this repository is built on the [Digital Earth Australia](https://docs.dea.ga.gov.au/) implementation of the [Open Data Cube](https://www.opendatacube.org/) software for accessing, managing, and analyzing large quantities of Earth observation (EO) data. 
The code currently runs on the [Digital Earth Australia Sandbox](https://docs.dea.ga.gov.au/setup/Sandbox/sandbox.html) infrastructure.

#### Python modules

Code in this repository is included in the `dea_coastlines` Python package which contains three main modules. These are intended to be run in the following order:

1. [`dea_coastlines.raster`](dea_coastlines/raster.py): This module conducts raster generation for DEA Coastlines. This analysis is processed on individual study area tiles to minimise peak memory usage.

    * Load stack of all available Landsat 5, 7 and 8 satellite imagery for a location using [ODC Virtual Products](https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Virtual_products.html)
    * Convert each satellite image into a remote sensing water index (e.g. MNDWI)
    * For each satellite image, model ocean tides into a 2 x 2 km grid based on exact time of image acquisition
    * Interpolate tide heights into spatial extent of image stack
    * Mask out high and low tide pixels by removing all observations acquired outside of 50 percent of the observed tidal range centered over mean sea level
    * Combine tidally-masked data into annual median composites from 1988 to the present representing the shoreline at approximately mean sea level

2. [`dea_coastlines.vector`](dea_coastlines/vector.py): This module conducts vector subpixel coastline extraction and rates of change statistics calculation. This analysis is processed on individual study area tiles to minimise peak memory usage.

    * Apply morphological extraction algorithms to mask annual median composite rasters to a valid coastal region
    * Extract shoreline vectors using subpixel waterline extraction ([Bishop-Taylor et al. 2019b](https://doi.org/10.3390/rs11242984))
    * Compute rates of coastal change at every 30 m along Australia's non-rocky coastline using linear regression
  
3. [`dea_coastlines.continental`](dea_coastlines/continental.py): This module combines tiled layers into seamless continental-scale vector files:

    * Combines multiple output shoreline and rates of change statistics point vectors into single continental datasets
    * Aggregates this data to produce a moving window coastal change hotspot dataset that summarises coastal change at regional and continental scale.

#### Command-line interface

These three modules have a command-line interface that can be used to automate each stage of the analysis. An example of these tools is provided in the following Jupyter Notebook:
* [DEA Coastlines generation using command line tools](notebooks/DEACoastlines_generation_CLI.ipynb)


For help using these command line tools, run:
```
python -m dea_coastlines.raster --help
python -m dea_coastlines.vector --help
python -m dea_coastlines.continental --help
```

#### Jupyter notebooks
An interactive walk-through of each step of the tiled raster and vector DEA Coastlines workflow is provided in the following Jupyter Notebooks. These notebooks can be used to assist in prototyping or troubleshooting:
* [DEA Coastlines raster generation](notebooks/DEACoastlines_generation_raster.ipynb)
* [DEA Coastlines vector generation](notebooks/DEACoastlines_generation_vector.ipynb)

---

## DEA Coastlines dataset

> _For the most up-to-date product metadata, visit the official [Geoscience Australia DEA Coastlines product description](https://cmi.ga.gov.au/data-products/dea/581/dea-coastlines-landsat)_

### Data access

#### Data download

To download DEA Coastlines data for the entire Australian coastline, visit the "Access" tab of the [Geoscience Australia DEA Coastlines product description](https://cmi.ga.gov.au/data-products/dea/581/dea-coastlines#access) and follow the instructions under "Access notes". Data is available in two formats:

* GeoPackage (recommended): suitable for QGIS; includes built-in symbology for easier interpretation
* ESRI Shapefiles: suitable for ArcMap and QGIS

#### Interactive map

To explore DEA Coastlines on an interactive map, visit the [Digital Earth Australia Maps platform](https://maps.dea.ga.gov.au/#share=s-DEACoastlines&playStory=1).

![Zooming to annual rates of change and plotting chart in DEA Maps](https://data.dea.ga.gov.au/projects/coastlines/DEACoastLines_DEAMaps_1.gif)


#### Loading DEA Coastlines data from the Web Feature Service (WFS) using Python

DEA Coastlines data can be loaded directly in a Python script or Jupyter Notebook using the DEA Coastlines Web Feature Service (WFS) and `geopandas`:

```
import geopandas as gpd

# Specify bounding box
ymax, xmin = -33.65, 115.28
ymin, xmax = -33.66, 115.30

# Set up WFS requests for annual coastlines & rates of change statistics
deacl_coastlines_wfs = f'https://geoserver.dea.ga.gov.au/geoserver/wfs?' \
                       f'service=WFS&version=1.1.0&request=GetFeature' \
                       f'&typeName=dea:coastlines&maxFeatures=1000' \
                       f'&bbox={ymin},{xmin},{ymax},{xmax},' \
                       f'urn:ogc:def:crs:EPSG:4326'
deacl_statistics_wfs = f'https://geoserver.dea.ga.gov.au/geoserver/wfs?' \
                       f'service=WFS&version=1.1.0&request=GetFeature' \
                       f'&typeName=dea:coastlines_statistics&maxFeatures=1000' \
                       f'&bbox={ymin},{xmin},{ymax},{xmax},' \
                       f'urn:ogc:def:crs:EPSG:4326'

# Load DEA Coastlines data from WFS using geopandas
deacl_coastlines_gdf = gpd.read_file(deacl_coastlines_wfs)
deacl_statistics_gdf = gpd.read_file(deacl_statistics_wfs)

# Ensure CRSs are set correctly
deacl_coastlines_gdf.crs = 'EPSG:3577'
deacl_statistics_gdf.crs = 'EPSG:3577'
```

### Loading DEA Coastlines data from the Web Feature Service (WFS) using R
DEA Coastlines data can be loaded directly into `R` using the DEA Coastlines Web Feature Service (WFS) and `sf`:

```
library(magrittr)
library(glue)
library(sf)

# Specify bounding box
xmin = 115.28
xmax = 115.30
ymin = -33.66
ymax = -33.65

# Read in DEA Coastlines annual coastline data, using `glue` to insert our bounding
# box into the string, and `sf` to  load the spatial data from the Web Feature
# Service and set the Coordinate Reference System to Australian Albers (EPSG:3577)
deacl_coastlines = "https://geoserver.dea.ga.gov.au/geoserver/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=dea:coastlines&maxFeatures=1000&bbox={ymin},{xmin},{ymax},{xmax},urn:ogc:def:crs:EPSG:4326" %>% 
  glue::glue() %>%
  sf::read_sf() %>% 
  sf::st_set_crs(3577)

# Read in DEA Coastlines rates of change statistics data
deacl_statistics = "https://geoserver.dea.ga.gov.au/geoserver/wfs?service=WFS&version=1.1.0&request=GetFeature&typeName=dea:coastlines_statistics&maxFeatures=1000&bbox={ymin},{xmin},{ymax},{xmax},urn:ogc:def:crs:EPSG:4326" %>% 
  glue::glue() %>%
  sf::read_sf() %>% 
  sf::st_set_crs(3577)
```

#### Jupyter Notebook
An [Introduction to DEA Coastlines](https://docs.dea.ga.gov.au/notebooks/DEA_datasets/DEA_Coastlines.html) Jupyter notebook providing additional useful tools for loading and analysing DEA Coastlines data can be found on the [DEA Notebooks repository](https://github.com/GeoscienceAustralia/dea-notebooks). This notebook is available on the interactive [DEA Sandbox learning and analysis environment](https://docs.dea.ga.gov.au/setup/Sandbox/sandbox.html) for easy access via a web browser.

---

## References
Bishop-Taylor, R., Nanson, R., Sagar, S., Lymburner, L. (2021). Mapping Australia's dynamic coastline at mean sea level using three decades of Landsat imagery. _Remote Sensing of Environment_, 267, 112734. Available: https://doi.org/10.1016/j.rse.2021.112734

Bishop-Taylor, R., Sagar, S., Lymburner, L., & Beaman, R. J. (2019a). Between the tides: Modelling the elevation of Australia's exposed intertidal zone at continental scale. _Estuarine, Coastal and Shelf Science_, 223, 115-128. Available: https://doi.org/10.1016/j.ecss.2019.03.006

Bishop-Taylor, R., Sagar, S., Lymburner, L., Alam, I., & Sixsmith, J. (2019b). Sub-pixel waterline extraction: Characterising accuracy and sensitivity to indices and spectra. _Remote Sensing_, 11(24), 2984. Available: https://doi.org/10.3390/rs11242984

Sagar, S., Roberts, D., Bala, B., & Lymburner, L. (2017). Extracting the intertidal extent and topography of the Australian coastline from a 28 year time series of Landsat observations. _Remote Sensing of Environment_, 195, 153-169. Available: https://doi.org/10.1016/j.rse.2017.04.009
