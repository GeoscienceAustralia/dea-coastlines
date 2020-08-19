![](https://github.com/GeoscienceAustralia/dea-notebooks/blob/develop/Supplementary_data/dea_logo_wide.jpg)

# Digital Earth Australia CoastLines

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

**License:** The code in this repository is licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0). Digital Earth Australia data is licensed under the [Creative Commons by Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/).

**Contact:** For assistance with any of the Python code or Jupyter Notebooks in this repository, please post a [Github issue](https://github.com/GeoscienceAustralia/DEACoastLines/issues/new). For questions or more information about this product, sign up to the [Open Data Cube Slack](https://join.slack.com/t/opendatacube/shared_invite/zt-d6hu7l35-CGDhSxiSmTwacKNuXWFUkg) and post on the [`#dea-coastlines`](https://app.slack.com/client/T0L4V0TFT/C018X6J9HLY/details/) channel.

---

**Digital Earth Australia CoastLines** is a continental dataset that includes annual shorelines and rates of coastal change along the entire Australian coastline from 1988 to the present. 

The product combines satellite data from Geoscience Australia's [Digital Earth Australia program](https://www.ga.gov.au/dea) with tidal modelling to map the typical location of the coastline at mean sea level for each year. The product enables trends of coastal erosion and growth to be examined annually at both a local and continental scale, and for patterns of coastal change to be mapped historically and updated regularly as data continues to be acquired. This allows current rates of coastal change to be compared with that observed in previous years or decades. 

The ability to map shoreline positions for each year provides valuable insights into whether changes to our coastline are the result of particular events or actions, or a process of more gradual change over time. This information can enable scientists, managers and policy makers to assess impacts from the range of drivers impacting our coastlines and potentially assist planning and forecasting for future scenarios. 

#### Applications
* Monitoring and mapping rates of coastal erosion along the Australian coastline 
* Prioritise and evaluate the impacts of local and regional coastal management based on historical coastline change 
* Modelling how coastlines respond to drivers of change, including extreme weather events, sea level rise or human development 
* Supporting geomorphological studies of how and why coastlines have changed across time 


## Data structure and features
The **DEA CoastLines** product contains three layers:

### Annual coastlines
Annual coastline vectors from 1988 to 2019 that represent the median or ‘typical’ position of the coastline at approximately mean sea level tide (0 m AHD).
   * Semi-transparent coastlines have low certainty due to either few non-cloudy satellite observations, or poor tidal modelling performance. 

![DEA CoastLines coastlines layer](visualisation/deacl_coastlines.JPG)

### Rates of change statistics
A point dataset providing robust rates of coastal change statistics for every 30 m along Australia’s non-rocky (clastic) coastlines. The most recent 2019 coastline is used as a baseline for measuring rates of change. By default, points are shown for significant rates of change only (p-value < 0.01, see sig_time below). The dataset contains the following attribute columns: 

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
   
![DEA CoastLines statistics layer](visualisation/deacl_statistics.JPG)

### Summary
A point layer giving the average rate of change (in metres per year) for significant statistics points within a moving 5 km window along the coastline. This is useful for visualising regional or continental-scale patterns of coastal change. 

![DEA CoastLines summary layer](visualisation/deacl_summary.JPG)

## Key limitations and caveats
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
