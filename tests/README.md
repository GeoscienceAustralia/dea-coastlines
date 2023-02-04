
Integration tests
=================


This directory contains integration tests that are run to verify that the entire Coastlines code runs correctly. The ``test_coastline.py`` file runs a simple Coastlines analysis over a single beach (Narrabeen Beach in northern Sydney), using the DEA Coastlines [Command Line Interface (CLI) tools](../notebooks/DEACoastlines_generation_CLI.ipynb) to run the raster, vector, continental layers and validation analysis steps.
# Validation


In addition to testing whether the code runs, we also run a small-scale validation of the results of the integration tests by comparing them to validation data from the [Narrabeen-Collaroy Beach Survey Program](https://doi.org/10.1038/sdata.2016.24). The ensures that the code both works, and generates sensible results. Note that this integration test validation is for a single site and a limited number of years - it is not intended to be representative of DEA Coastline's accuracy overall.
## Latest integration test validation results


The latest integration test completed at **2023-02-04 19:01**. Compared to the previous run, it had an RMSE accuracy of **4.57 m (no change)**, an MAE accuracy of **3.68 m (no change)**, and a Pearson correlation of **0.996 (no change)**.

<img src="stats_tests.png" width="850"/>