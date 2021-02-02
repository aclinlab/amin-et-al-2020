# Amin et al., 2020
 Code used in [Amin et al. 2020](https://elifesciences.org/articles/56954) to measure distances between synapses along neurite skeletons
 
### Fetching data from hemibrain server

1. Install neuprint-python from:
https://github.com/connectome-neuprint/neuprint-python
(On a Mac, you may need to install Command Line Tools.)
2. Run getAPLskelAndSynapses_v1.1.py

`python3 getAPLskelAndSynapses_v1.1.py`

This will create the files APLskelv1.1.csv (APL's skeleton), APLtoKCv1.1.csv (list of all APL-KC synapses), KCtoAPLv1.1.csv (list of all KC-APL synapses). These are already provided in this repository under the folder 'data'.

### Analyzing APL's skeleton and synapses

See the file `runAPL.m` to see what commands to run to reproduce Figure 8 and related figure supplements from Amin et al. 

### General explanation of code structure

The class `neurSkel` is intended to represent the skeleton of any neuron from the hemibrain connectome. The class `APLskel` extends this class with some APL-specific functions. Most core functions are contained in these classes but much analysis for the paper was done with standalone scripts as detailed in runAPL.m.

### Code dependencies

Mapping the connectome APL onto the standard APL of our data requires the class `activityMap` in [https://github.com/aclinlab/calcium-imaging](https://github.com/aclinlab/calcium-imaging)
