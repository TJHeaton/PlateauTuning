# PlateauTuning

This folder will contain the R code used to test plateau tuning in the Climate of the Past paper:

Bard, E. and Heaton, T. J.: On the tuning of plateaus in atmospheric and oceanic 14C records to derive calendar chronologies of deep-sea cores and records of 14C marine reservoir age changes, Clim. Past Discuss. [preprint], https://doi.org/10.5194/cp-2020-164, in review, 2021.

Plateau tuning is a method proposed by Sarnthein et al. (2020) to create chronologies for ocean sediment cores by matching suites of hypothesized atmospheric 14C plateaus to hpothesized 14C plateaus in the various sediment cores. Our paper aims to discuss our strong reservations about the reliability of the plateau tuning approach.

The code is provided in R and consists of three files:

1) PlateauMatchingFunctions.R - code to estimate the local gradient of the 14C data  
2) SimulatedAtmosphericCores.R - code to create 3 simulated pseudo-atmospheric 14C cores and recreate Figs 5 and 6  
3) SimulatedMarineCores - code to create 2 simulated pseudo-marine 14C cores and recreate Fig 7 

You will also need to download the three files:

1) IntCal20TreeRing.csv - pointwise estimate of the tree-ring based section of IntCal20 
2) Suigetsu2013.csv - Lake Suigetsu 14C data 
3) Cariaco2013Raw.csv - Cariaco 14C data

We will also provide a DOI for this code with the paper
