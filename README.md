## These scripts reproduce the analysis and figures from Trautmann et al. 2023, ERC
 [![DOI](https://zenodo.org/badge/658690923.svg)](https://zenodo.org/badge/latestdoi/658690923)

- they were written and run with MATLAB R2019a

- open MATLAB from the folder, e.g. by opening `getStarted.mat`
- `getStarted.mat` sets all paths and runs all other scripts
- if paths are set correctly, all files in the main folder can be run independently

- to execute the scripts, please also download the _input_ and _model_output_ data from the linked zenodo repository, and save them in a _data/_ subfolder
- the folder further contains:
	* _data/input_ 		- forcing data and observational constraints 
 						- ancillary data such as the regional classification
 						- all data is provided in *.mat files and adjusted to the study area and considered time period
 						- all original data is publical available, please see the data availability of Trautmann et al. 2022, HESS & Trautmann et al. 2023, ERC
	* _data/model_output_ 	- output from different model experiments, including simulated fluxes & states as .netcdf files; information on model settings and model code

