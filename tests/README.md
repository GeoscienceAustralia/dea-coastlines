
Integration tests
=================


This directory contains integration tests that are run to verify that the entire Coastlines code runs correctly. The ``test_coastline.py`` file runs a simple Coastlines analysis over a single beach (Narrabeen Beach in northern Sydney), using the DEA Coastlines [Command Line Interface (CLI) tools](../notebooks/DEACoastlines_generation_CLI.ipynb) to run the raster, vector, continental layers and validation analysis steps.
# Validation


In addition to testing whether the code runs without errors, we also run a small-scale validation of the results of the integration tests by comparing them to validation data from the [Narrabeen-Collaroy Beach Survey Program](https://doi.org/10.1038/sdata.2016.24). The ensures that the code both works, and generates sensible results. Note that this integration test validation is for a single site and a limited number of years - it is not intended to be representative of DEA Coastline's accuracy overall.
## Latest integration test validation results


The latest integration test completed at **2023-02-04 21:28**. Compared to the previous run, it had an RMSE accuracy of **6.19 m (:small_red_triangle_down: decreased by -0.3)**, an MAE accuracy of **4.38 m (:small_red_triangle_down: decreased by -0.1)**, and a Pearson correlation of **0.990 (heavy_minus_sign: no change)**.

<img src="stats_tests.png" width="950"/>