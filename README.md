![Digital Earth Africa Coastlines](https://github.com/digitalearthafrica/deafrica-sandbox-notebooks/raw/main/Supplementary_data/Github_banner.jpg)

# Digital Earth Africa Coastlines

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.rse.2021.112734-0e7fbf.svg)](https://doi.org/10.1016/j.rse.2021.112734)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**License:** The code in this repository is licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0). Digital Earth Africa data is licensed under the [Creative Commons by Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).

**Contact:** For assistance with any of the Python code or Jupyter Notebooks in this repository, please post a [Github issue](https://github.com/GeoscienceAustralia/DEACoastLines/issues/new). For questions or more information about this workflow, email Robbi.BishopTaylor@ga.gov.au.

**To cite:** 
> Bishop-Taylor, R., Nanson, R., Sagar, S., Lymburner, L. (2021). Mapping Australia's dynamic coastline at mean sea level using three decades of Landsat imagery. _Remote Sensing of Environment_, 267, 112734. Available: https://doi.org/10.1016/j.rse.2021.112734

> Bishop-Taylor, R., Sagar, S., Lymburner, L., Alam, I., Sixsmith, J. (2019). Sub-pixel waterline extraction: characterising accuracy and sensitivity to indices and spectra. _Remote Sensing_, 11 (24):2984. Available: https://doi.org/10.3390/rs11242984

> Nanson, R., Bishop-Taylor, R., Sagar, S., Lymburner, L. (2021). In press. Geomorphic insights into Australia’s coastal change using a national dataset derived from the multi-decadal Landsat archive. _Estuarine, Coastal and Shelf Science_.

---

**Digital Earth Africa Coastlines** is a modified implementation of Digital Earth Australia (DEA) Coastlines, a continental remote sensing method designed to generate annual shorelines and rates of coastal change along the entire Australian coastline from 1988 to the present. 

The DEA Coastlines method combines satellite data with tidal modelling to map the typical location of the coastline at mean sea level for each year. The product enables trends of coastal erosion and growth to be examined annually at both a local and continental scale, and for patterns of coastal change to be mapped historically and updated regularly as data continues to be acquired. This allows current rates of coastal change to be compared with that observed in previous years or decades. 

The ability to map shoreline positions for each year provides valuable insights into whether changes to our coastline are the result of particular events or actions, or a process of more gradual change over time. This information can enable scientists, managers and policy makers to assess impacts from the range of drivers impacting our coastlines and potentially assist planning and forecasting for future scenarios. 

---

## Table of contents
* [Repository code](#repository-code)
    * [Python modules](#python-modules)
        * [Jupyter notebooks](#jupyter-notebooks)
    * [Running a DE Africa Coastlines analysis using the command-line interface (CLI)](#running-a-de-africa-coastlines-analysis-using-the-command-line-interface-cli)
        * [Analysis outputs](#analysis-outputs)
* [References](#references)

---

## Repository code
The code in this repository is built on the Digital Earth Africa implementation of the [Open Data Cube](https://www.opendatacube.org/) software for accessing, managing, and analyzing large quantities of Earth observation (EO) data. 
The code currently runs on the [Digital Earth Africa Sandbox](https://sandbox.digitalearth.africa/) infrastructure.

#### Python modules

Code in this repository is included in the `dea_coastlines` Python package which contains three main modules. These are intended to be run in the following order:

1. [`dea_coastlines.raster`](dea_coastlines/raster.py): This module conducts raster generation for DEA Coastlines. This analysis is processed on individual study area tiles to minimise peak memory usage.

    * Load stack of all available Landsat 5, 7 and 8 satellite imagery for a location using [ODC Virtual Products](https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Virtual_products.html)
    * Convert each satellite image into a remote sensing water index (e.g. MNDWI)
    * For each satellite image, model ocean tides into a tidal modelling grid based on exact time of image acquisition
    * Interpolate tide heights into spatial extent of image stack
    * Mask out high and low tide pixels by removing all observations acquired outside of 50 percent of the observed tidal range centered over mean sea level
    * Combine tidally-masked data into annual median composites representing the most representative position of the shoreline at approximately mean sea level each year

2. [`dea_coastlines.vector`](dea_coastlines/vector.py): This module conducts vector subpixel coastline extraction and rates of change statistics calculation. This analysis is processed on individual study area tiles to minimise peak memory usage.

    * Apply morphological extraction algorithms to mask annual median composite rasters to a valid coastal region
    * Extract shoreline vectors using subpixel waterline extraction ([Bishop-Taylor et al. 2019b](https://doi.org/10.3390/rs11242984))
    * Compute rates of coastal change at every 30 m along the coastline using linear regression
  
3. [`dea_coastlines.continental`](dea_coastlines/continental.py): This module combines tiled layers into seamless continental-scale vector files:

    * Combines multiple output shoreline and rates of change statistics point vectors into single continental datasets
    * Aggregates this data to produce a moving window coastal change hotspot dataset that summarises coastal change at regional and continental scale.
    
    
#### Jupyter notebooks
An interactive walk-through of each step of the tiled raster and vector DEA Coastlines workflow is provided in the following Jupyter Notebooks. These notebooks can be used to assist in prototyping or troubleshooting:
* [DE Africa Coastlines raster generation](notebooks/DEAfricaCoastlines_generation_raster.ipynb)
* [DE Africa Coastlines vector generation](notebooks/DEAfricaCoastlines_generation_vector.ipynb)

### Running a DE Africa Coastlines analysis using the command-line interface (CLI)

These three modules have a command-line interface that can be used to automate each stage of the analysis. An example of using these tools is provided in the following Jupyter Notebook:
* [DE Africa Coastlines generation using command line tools](notebooks/DEAfricaCoastlines_generation_CLI.ipynb)

For help using these command line tools, run:
```
python -m dea_coastlines.raster --help
python -m dea_coastlines.vector --help
python -m dea_coastlines.continental --help
```

#### Analysis outputs
Files generated by DE Africa Coastlines are exported to the `data` directory. 

Temporary raster and vector outputs produced by [`dea_coastlines.raster`](dea_coastlines/raster.py) and [`dea_coastlines.vector`](dea_coastlines/vector.py) for each study area grid cell are exported to:
```
data/interim/raster/{unique_analysis_name}/{unique_analysis_name}_{study_area_name}
data/interim/vector/{unique_analysis_name}/{unique_analysis_name}_{study_area_name}
```

Once all study area grid cells have been processed, these are combined into continental-scale output vector files using [`dea_coastlines.continental`](dea_coastlines/continental.py). These final outputs are exported to:
```
data/processed/{unique_analysis_name}/DEAfricaCoastlines_annualshorelines_{unique_analysis_name}.shp
data/processed/{unique_analysis_name}/DEAfricaCoastlines_ratesofchange_{unique_analysis_name}.shp
data/processed/{unique_analysis_name}/DEAfricaCoastlines_hotspots_{unique_analysis_name}_{hotspot_radius}.shp
```

---

## References
Bishop-Taylor, R., Nanson, R., Sagar, S., Lymburner, L. (2021). Mapping Australia's dynamic coastline at mean sea level using three decades of Landsat imagery. _Remote Sensing of Environment_, 267, 112734. Available: https://doi.org/10.1016/j.rse.2021.112734

Bishop-Taylor, R., Sagar, S., Lymburner, L., & Beaman, R. J. (2019a). Between the tides: Modelling the elevation of Australia's exposed intertidal zone at continental scale. _Estuarine, Coastal and Shelf Science_, 223, 115-128. Available: https://doi.org/10.1016/j.ecss.2019.03.006

Bishop-Taylor, R., Sagar, S., Lymburner, L., Alam, I., & Sixsmith, J. (2019b). Sub-pixel waterline extraction: Characterising accuracy and sensitivity to indices and spectra. _Remote Sensing_, 11(24), 2984. Available: https://doi.org/10.3390/rs11242984

Sagar, S., Roberts, D., Bala, B., & Lymburner, L. (2017). Extracting the intertidal extent and topography of the Australian coastline from a 28 year time series of Landsat observations. _Remote Sensing of Environment_, 195, 153-169. Available: https://doi.org/10.1016/j.rse.2017.04.009

Nanson, R., Bishop-Taylor, R., Sagar, S., Lymburner, L. (2021). In press. Geomorphic insights into Australia’s coastal change using a national dataset derived from the multi-decadal Landsat archive. _Estuarine, Coastal and Shelf Science_.
