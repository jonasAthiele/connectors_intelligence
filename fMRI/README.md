
## fMRI analysis

### Data

`aCC24Treg_ts_20p2p7_pkg` - folder with preprocessed fMRI data (3 runs Raven (RPM), 1 run rest)

`onsets` - folder with task onsets of each participant
	
`node_network_assignment_200_17.mat` - assignment of the 200 nodes to 17 networks

`table_descriptives.mat` - age, sex of participants

`table_accuracy.mat` - RPM performance of participants

`Schaefer200_Pfit_CoreRegions_Annotated.csv` - Schaefer nodes assigned to P-FIT-regions

`table_accuracy.mat` - RPM performance of participants

`table_handedness.mat` - handedness of participants

### Main scripts	

`get_fcs.mat`

  - filter time-ranges of problem solving 
  - estimate fcs
  - get parameters of framewise displacement for later exclusions due to head motion	

`get_results.mat`

  - exclude due to high motion and too short raven runs
  - compute centrality measures 
  - get associations centrality measures and intelligence
  - get info of subjects
  - get main figures and figures for supplement
	

### External functions	

`parcplot` - [parc_plotter](https://github.com/faskowit/parc_plotter) (needed for plotting)
		
`threshold_proportional` - from [Brain Connectivity Toolbox, (BCT)](https://sites.google.com/site/bctnet/home), thresholds fc

`participation_coef` - from the BCT, computes participation coefficient

`fisherZTransform` - from the BCT, Fisher's z transformation

`apply_CanonicalHRF.mat` - concolves onsets with HRF to get timings of expected response

`slanCM.mat` & `slanCM_Data.mat` - function & data for creating colormaps (https://de.mathworks.com/matlabcentral/fileexchange/120088-200-colormap)
