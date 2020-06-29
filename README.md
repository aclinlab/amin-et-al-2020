# hemibrain-synapse-spacing
 Code used in Amin et al. 2020 to measure distances between synapses along neurite skeletons

### Fetching data from hemibrain server

1. Install neuprint-python from:
https://github.com/connectome-neuprint/neuprint-python
(On a Mac, you may need to install Command Line Tools.)
2. Run getAPLskelAndSynapses.py

`python3 getAPLskelAndSynapses.py`

This will create the files APLskel.csv (APL's skeleton), APLtoKC.csv (list of all APL-KC synapses), KCtoAPL.csv (list of all KC-APL synapses). These are already provided in this repository under the folder 'data'.

### Analyzing APL's skeleton and synapses

See the file `runAPL.m` to see what commands to run to reproduce Figure 8 and related figure supplements from Amin et al. 

### General explanation of code structure

The class `neurSkel` is intended to represent the skeleton of any neuron from the hemibrain connectome. The class `APLskel` extends this class with some APL-specific functions. Most core functions are contained in these classes but much analysis for the paper was done with standalone scripts as detailed in runAPL.m.

### Code dependencies

Mapping the connectome APL onto the standard APL of our data requires the class `activityMap` in [https://github.com/aclinlab/calcium-imaging](https://github.com/aclinlab/calcium-imaging)