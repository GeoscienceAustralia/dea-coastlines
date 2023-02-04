
Integration tests
=================


This directory contains integration tests that are run to verify that the entire Coastlines code runs correctly. The ``test_coastline.py`` file runs a simple Coastlines analysis over a single beach (Narrabeen Beach in northern Sydney), using the DEA Coastlines [Command Line Interface (CLI) tools](../notebooks/DEACoastlines_generation_CLI.ipynb) to run the raster, vector, continental layers and validation analysis steps.
# Validation


In addition to testing whether the code runs without errors, we also run a small-scale validation of the results of the integration tests by comparing them to validation data from the [Narrabeen-Collaroy Beach Survey Program](https://doi.org/10.1038/sdata.2016.24). The ensures that the code both works, and generates sensible results. Note that this integration test validation is for a single site and a limited number of years - it is not intended to be representative of DEA Coastline's accuracy overall.
## Latest integration test validation results


The latest integration test completed at **2023-02-04 19:23**. Compared to the previous run, it had an RMSE accuracy of **6.52 m (improved by 1.9499999999999993)**, an MAE accuracy of **4.55 m (improved by 0.8699999999999997)**, and a Pearson correlation of **0.99 (deteriorated by -0.006000000000000005)**.

<img src="stats_tests.png" width="950"/>