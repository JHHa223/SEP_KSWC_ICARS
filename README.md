This repository includes the model for predicting solar energetic particles at Korea Space Weather Center (KSWC).

The original contents were primarily presented at International Conference on Advanced Remote Sensing (ICARS).

The detailed description is provided as follows:

Sample data
1. 0209_proton_1day.txt : A sample observed data measured by GOES geostationary satellite downloaded from the SWPC NOAA.
2. spectrum_LB3.npy : A sample particle distribution function produced by the shock in quasi-stationary state.
3. goes_proton_sequence.py : A sample code for reading the observed data (i.e., 0209_proton_1day.txt).

Model
4. distribution_plots_ensemble_final.py : A model code for generating the time-sequence of proton flux at L1-point. A sample file of observed data is used for comparing the model prediction and ground truth.
