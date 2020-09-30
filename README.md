![Digital Earth Australia Coastlines](visualisation/deacoastlines_example_header.gif)

# Digital Earth Australia Coastlines

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**Access:** [Explore DEA Coastlines on the interactive Digital Earth Australia Maps platform](https://maps.dea.ga.gov.au/#share=s-DEACoastlines&playStory=1)

**License:** The code in this repository is licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0). Digital Earth Australia data is licensed under the [Creative Commons by Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).

**Product description and metadata:** For the most up-to-date information about this product, visit the official [Geoscience Australia DEA Coastlines product description](https://cmi.ga.gov.au/data-products/dea/581/dea-coastlines-landsat)

**Contact:** For assistance with any of the Python code or Jupyter Notebooks in this repository, please post a [Github issue](https://github.com/GeoscienceAustralia/DEACoastLines/issues/new). For questions or more information about this product, email dea@ga.gov.au or sign up to the [Open Data Cube Slack](https://join.slack.com/t/opendatacube/shared_invite/zt-d6hu7l35-CGDhSxiSmTwacKNuXWFUkg) and post on the [`#dea-coastlines`](https://app.slack.com/client/T0L4V0TFT/C018X6J9HLY/details/) channel.


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
* [Data access](#data-access)
    * [Digital Earth Australia Maps](#digital-earth-australia-maps)
    * [Loading DEA Coastlines data using Web Feature Service (WFS)](#loading-dea-coastlines-data-using-web-feature-service-wfs)
* [DEA Coastlines dataset](#dea-coastlines-dataset)
    * [Annual coastlines](#annual-coastlines)
    * [Rates of change statistics](#rates-of-change-statistics)
    * [Coastal change hotspots](#coastal-change-hotspots)
* [Caveats and limitations](#caveats-and-limitations)
* [References](#references)

---

## Repository code
The code in this repository is built on the [Digital Earth Australia](https://docs.dea.ga.gov.au/) implementation of the [Open Data Cube](https://www.opendatacube.org/) software for accessing, managing, and analyzing large quantities of Earth observation (EO) data. 
The software currently runs on [Australia's National Computational Infrastructure (NCI)](https://nci.org.au/). 
Instructions for setting up an account on the NCI's Virtual Desktop Infrastructure or Gadi supercomputer [can be found here](https://docs.dea.ga.gov.au/setup/NCI/README.html).

This repository contains three main scripts (and corresponding Jupyter notebooks) that are intended to be run in the following order:

1. [`deacoastlines_generation.py`](deacoastlines_generation.py)/[`DEACoastlines_generation.ipynb`](DEACoastlines_generation.ipynb): This code conducts raster generation for DEA Coastlines:

    * Load stack of all available Landsat 5, 7 and 8 satellite imagery for a location using [ODC Virtual Products](https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Virtual_products.html)
    * Convert each satellite image into a remote sensing water index (MNDWI)
    * For each satellite image, model ocean tides into a 2 x 2 km grid based on exact time of image acquisition
    * Interpolate tide heights into spatial extent of image stack
    * Mask out high and low tide pixels by removing all observations acquired outside of 50 percent of the observed tidal range centered over mean sea level
    * Combine tidally-masked data into annual median composites from 1988 to the present representing the coastline at approximately mean sea level

2. [`deacoastlines_statistics.py`](deacoastlines_statistics.py)/[`DEACoastlines_statistics.ipynb`](DEACoastlines_statistics.ipynb): This code conducts vector subpixel coastline extraction:

    * Apply morphological extraction algorithms to mask annual median composite rasters to a valid coastal region
    * Extract waterline vectors using subpixel waterline extraction (Bishop-Taylor et al. 2019b)
    * Compute rates of coastal change at every 30 m along Australia's non-rocky coastlines using linear regression
  
3. [`deacoastlines_summary.py`](deacoastlines_summary.py): This script combines individual datasets into continental DEA Coastlines layers:

    * Combines output coastline and rates of change statistics point vectors into single continental datasets
    * Aggregates this data to produce a moving window coastal change hotspot dataset that summarises coastal change at regional and continental scale.

An additional Jupyter notebook provides useful tools for analysing DEA Coastlines data:

* [`DEACoastlines_tools.ipynb`](DEACoastlines_tools.ipynb): 

    * Selecting and loading DEA Coastlines data using an interactive map
    * Interactively drawing a transect across DEA Coastlines annual coastlines and generating a plot of coastal change through time
    * Interactively plotting the distribution of retreating and growing coastlines within a selected region

---

## Data access

### Digital Earth Australia Maps

To explore DEA Coastlines on the interactive Digital Earth Australia Maps platform, visit the link below:

[https://maps.dea.ga.gov.au/#share=s-DEACoastlines](https://maps.dea.ga.gov.au/#share=s-DEACoastlines&playStory=1)

To add DEA Coastlines to DEA Maps manually:

1. Open DEA Maps: http://maps.dea.ga.gov.au/ 
2. Select `Add data` on the top-left. 
3. Select `Coastal > Digital Earth Australia Coastlines > Digital Earth Australia Coastlines`
4. Click blue 'Add to the map' button on top-right. 

By default, the map will show hotspots of coastal change at continental scale. 
Red dots represent retreating coastlines (e.g. erosion), while blue dots indicate seaward growth. The larger the dots and the brighter the colour, the more coastal change that is occurring at the location. 

More detailed rates of change will be displayed as you zoom in. 
To view a time series chart of how an area of coastline has changed over time, click on any labelled point (press "Expand" on the pop-up for more detail):

![Zooming to annual rates of change and plotting chart in DEA Maps](https://data.dea.ga.gov.au/projects/coastlines/DEACoastLines_DEAMaps_1.gif)

Zoom in further to view individual annual coastlines: 

![Zooming to annual coastlines in DEA Maps](https://data.dea.ga.gov.au/projects/coastlines/DEACoastLines_DEAMaps_2.gif)
 
> Note: To view a DEA Coastlines layer that is not currently visible (e.g. rates of change statistics at full zoom), each layer can be added to the map individually from the `Coastal > Digital Earth Australia Coastlines > Supplementary data` directory.

### Loading DEA Coastlines data using Web Feature Service (WFS)

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

---

## DEA Coastlines dataset

> _For the most up-to-date product metadata, visit the official [Geoscience Australia DEA Coastlines product description](https://cmi.ga.gov.au/data-products/dea/581/dea-coastlines-landsat)_

The **DEA Coastlines** product contains three layers:

### Annual coastlines
Annual coastline vectors from 1988 to 2019 that represent the median or ‘typical’ position of the coastline at approximately mean sea level tide (0 m AHD) for each year.
   * Semi-transparent coastlines have low certainty due to either few non-cloudy satellite observations, poor tidal modelling performance, or aerosol issues (see [Caveats and limitations](#caveats-and-limitations))

![DEA Coastlines coastlines layer](visualisation/deacoastlines_example_coastlines.JPG)

### Rates of change statistics
A point dataset providing robust rates of coastal change statistics for every 30 m along Australia’s non-rocky (clastic) coastlines.
The most recent 2019 coastline is used as a baseline for measuring rates of change.

![DEA Coastlines statistics layer](visualisation/deacoastlines_example_statistics.JPG)

On the [interactive DEA Coastlines web map](https://maps.dea.ga.gov.au/#share=s-DEACoastlines&playStory=1), points are shown for locations with statistically significant rates of change only (p-value < 0.01, see `sig_time` below). 
Each point shows annual rates of change (in metres per year; see `rate_time` below), and an estimate of uncertainty in this rate in brackets (95% confidence interval; see `se_time` below). 
For example, there is a 95% chance that the rate of erosion at a point with the following label is between -9.0 and -11.0 metres per year: 

> -10.0 m (± 1.0 m)

The rates of change statistics dataset contains the following attribute columns that can be accessed by clicking on labelled points in the web map: 

##### Annual coastline distances
   * `dist_1990`, `dist_1991` etc: Annual coastline positions/distances (in metres) relative to the 2019 baseline coastline. Negative values indicate that an annual coastline was located inland of the 2019 baseline coastline.    
   
##### Rates of change statistics
   * `rate_time`: Annual rates of change (in metres per year) calculated by linearly regressing all annual coastline distances against time. Negative values indicate retreat, while positive values indicate growth. 
   * `sig_time`: Significance (p-value) of the linear relationship between annual coastline distances and time. Small values (e.g. p-value < 0.01 or 0.05) may indicate a coastline is undergoing consistent coastal change through time. 
   * `se_time`: Standard error (in metres) of the linear relationship between annual coastline distances and time. This can be used to generate confidence intervals around the rate of change given by rate_time (e.g. 95% confidence interval = rate_time * 1.96)
   * `outl_time`: Individual annual coastlines are noisy estimators of coastline position that can be influenced by environmental conditions (e.g. clouds, breaking waves, sea spray) or modelling issues (e.g. poor tidal modelling results or limited clear satellite observations). To obtain robust rates of change, outlying years are excluded using a robust outlier detection algorithm, and recorded in this column.
   
##### Climate driver statistics
   * `rate_soi`: Annual rates of change (in metres per year) calculated by linearly regressing all annual coastline distances against the Southern Oscillation Index (SOI). Negative values indicate retreat during La Nina years.
   * `sig_soi`: Significance (p-value) of the linear relationship between annual coastline distances and SOI. 
   * `se_soi`: Standard error (in metres) of the linear relationship between annual coastline distances and SOI.
   * `outl_soi`: A list of any years excluded from the SOI regression by the robust outlier detection algorithm.

##### Other derived statistics
   * `retreat`, `growth`: True/False columns indicating whether a shoreline was retreating (i.e. moving inland) or growing (i.e. moving seaward) based on the rate_time column.
   * `sce`: Shoreline Change Envelope (SCE). A measure of the maximum change or variability across all annual coastlines, calculated by computing the maximum distance between any two annual coastlines (excluding outliers).
   * `nsm`: Net Shoreline Movement (NSM). The distance between the oldest (1988) and most recent (2019) annual coastlines (excluding outliers). Negative values indicate the shoreline retreated between the oldest and most recent coastline; positive values indicate growth.
   * `max_year`, `min_year`: The year that annual coastlines were at their maximum (i.e. located furthest towards the ocean) and their minimum (i.e. located furthest inland) respectively (excluding outliers).
   * `breaks`: An experimental list of any years identified as non-linear breakpoints in the time series. This can be useful for verifying that a significant trend is indeed linear, or identifying areas of rapid non-linear change (e.g. associated with coastal development or management).

### Coastal change hotspots
A point layer giving the average rate of change (in metres per year) for significant statistics points within a moving 5 km window along the coastline. This is useful for visualising regional or continental-scale patterns of coastal change. 

![DEA Coastlines coastal change hotspots](visualisation/deacoastlines_example_summary.JPG)

---

## Caveats and limitations
> _For the most up-to-date information about caveats and limitations, visit the official [Geoscience Australia DEA Coastlines product description](https://cmi.ga.gov.au/data-products/dea/581/dea-coastlines-landsat)_

* Rates of change statistics may be inaccurate or invalid for some complex mouthbars, or other coastal environments undergoing rapid non-linear change through time. In these regions, it is advisable to visually assess the underlying annual coastline data when interpreting rates of change to ensure these values are fit-for-purpose. Regions significantly affected by this issue include:
    * Cambridge Gulf, Western Australia
    * Joseph Bonaparte Gulf, Western Australia/Northern Territory
* Annual coastlines may be less accurate in regions with complex tidal dynamics or large tidal ranges, and low-lying intertidal flats where small tidal modelling errors can lead to large horizontal offsets in coastline positions. Annual coastline accuracy in intertidal environments may also be reduced by the influence of wet muddy substrate or intertidal vegetation, which can make it difficult to extract a single unambiguous coastline (Bishop-Taylor et al. 2019a, 2019b). It is anticipated that future versions of this product will show improved results due to integrating more advanced methods for waterline detection in intertidal regions, and through improvements in tidal modelling methods. Regions significantly affected by intertidal issues include:
    * The Pilbara coast, Western Australia from Onslow to Pardoo
    * The Mackay region, Queensland from Proserpine to Broad Sound
    * The upper Spencer Gulf, South Australia from Port Broughton to Port Augusta
    * Western Port Bay, Victoria from Tooradin to Pioneer Bay
    * Hunter Island Group, Tasmania from Woolnorth to Perkins Island
    * Moreton Bay, Queensland from Sandstone Bay to Wellington Point
* Coastlines may be noisier and more difficult to interpret in regions with low availability of satellite observations caused by persistent cloud cover. In these regions it can be difficult to obtain the minimum number of clear satellite observations required to generate clean, noise-free annual coastlines. Regions significantly affected by cloud cover issues include:
    * South-western Tasmania from Macquarie Heads to Southport
* In some urban locations, the spectra of bright white buildings located near the coastline may be inadvertently confused with water, causing a land-ward offset from true coastline positions. 
* Some areas of extremely dark and persistent shadows (e.g. steep coastal cliffs across southern Australia) may be inadvertently mapped as water, resulting in a landward offset from true coastline positions. 
* 1991 and 1992 coastlines are currently affected by aerosol-related issues caused by the 1991 Mount Pinatubo eruption. These coastlines should be interpreted with care, particularly across northern Australia. 

---

## References
Bishop-Taylor, R., Sagar, S., Lymburner, L., & Beaman, R. J. (2019a). Between the tides: Modelling the elevation of Australia's exposed intertidal zone at continental scale. _Estuarine, Coastal and Shelf Science_, 223, 115-128. Available: https://doi.org/10.1016/j.ecss.2019.03.006

Bishop-Taylor, R., Sagar, S., Lymburner, L., Alam, I., & Sixsmith, J. (2019b). Sub-pixel waterline extraction: Characterising accuracy and sensitivity to indices and spectra. _Remote Sensing_, 11(24), 2984. Available: https://doi.org/10.3390/rs11242984

Sagar, S., Roberts, D., Bala, B., & Lymburner, L. (2017). Extracting the intertidal extent and topography of the Australian coastline from a 28 year time series of Landsat observations. _Remote Sensing of Environment_, 195, 153-169. Available: https://doi.org/10.1016/j.rse.2017.04.009
