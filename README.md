# EARS 202 Final Project
**Constraining a sub-annual chronology for the Denali ice core**
*Jacob Chalif*

Our general method is to leverage a known relationship between ice core water isotope variability and precipitation-weighted temperature at the cloud condensation level in order to develop a sub-annual chronology for the Denali ice core. These programs align ice core water isotope data to reanalysis precipitation weighted temperature (PWT) via 2 methods: **1: Dynamic time warping (DTW)** and **2: Markov Chain Monte Carlo (MCMC)**.

Most of this analysis was done in Julia. The DTW approach is found in `dtw_timewarping.jl` while the MCMC approach is found in `mc_timewarping.jl`. Some plotting was done in Matlab, code for which can be found in `Plots.m`.