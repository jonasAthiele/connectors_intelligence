## EEG analysis

### Data

`rest` - rest raw data can be obtained from [Ociepka et al. 2022](https://osf.io/kv2sx/) (Rest_1.zip)

`raven` - Raven (RPM) data (mostly preprocessed) can be obtained from [Ociepka et al. 2023](https://osf.io/htrsg/)

`1_outcome.set` & `1_outcome_fdt` - arbitrary RPM eeg set from `raven` used as template for topoplots

`descr_behav_data` - descriptive and behavioral data of subjects (age, sex, RPM sum scores)

`channel_locations_64.ced` - electrode locations (coordinates)

	

### Main scripts

`eeg_analysis.mat`

- preprocessing EEG data
- get MSE of rest and RPM trials for theta band
- save measures in "brain_measures.mat" files for each subject 

`get_results_eeg.mat`

- get association between MSE measures and RPM sum scores
- make plots of results
- outputs "info_subjects" & "info_subjects_all"

`compute_partial_spearman.mat`
- function for computing partial Spearman correlation

`cluster_permutation_test_signed.mat` & `permute_corr_func.mat`
- functions for cluster based permutation test

### External functions

`multiscaleSampleEntropy.mat` - function for computing MSE (https://www.mathworks.com/matlabcentral/fileexchange/62706-multiscale-sample-entropy)

`slanCM.mat` & `slanCM_Data.mat` - function & data for creating colormaps (https://de.mathworks.com/matlabcentral/fileexchange/120088-200-colormap)
