## EEG analysis

### Data:

`rest` - rest raw data can be obtained from [Ociepka et al. 2022](https://osf.io/kv2sx/) (Rest_1.zip)

`raven` - Raven (RPM) data (mostly preprocessed) can be obtained from [Ociepka et al. 2023](https://osf.io/htrsg/)

`descr_behav_data` - descriptive and behavioral data of subjects (age, sex, RPM sum scores)

`channel_locations_64.ced` - electrode locations (coordinates)

`1_outcome.set` & `1_outcome_fdt` - arbitrary RPM eeg set used as template for topoplots
	

### Main scripts:

`eeg_analysis.mat`

	- preprocessing EEG data
	- get FCs of rest and RPM trials for theta band
	- save measures in "brain_measures.mat" files for each subject 


`get_results_eeg.mat`

	- compute centrality measures
	- get association between centrality measures and RPM sum scores
	- make plots of results
	- outputs "info_subjects" & "info_subjects_all" & "results_degree_uncorrected"
		& "results_participation_uncorrected"

`centrality_eeg_50.mat` - interim result: centrality measures


### External functions

`threshold_proportional.mat` - function from the [Brain Conenctivity Toolbox (BCT)](https://sites.google.com/site/bctnet/home) to threshold FC

`community_louvain.mat` & `consensus_und.mat` & `agreement.mat` - functions from the BCT to compute partitions

`participation_coef.mat` - function of BCT for computing participation coefficient

`fdr_bh` - function for FDR correction

`slanCM.mat` & `slanCM_Data.mat` - external function & data for creating colormaps (https://de.mathworks.com/matlabcentral/fileexchange/120088-200-colormap)
