
function [ taskdesignmat_hrf ] = apply_CanonicalHRF(taskdesignmat, tr_seconds, upsampleFactor) 
% adapted from Michael Cole lab:
% Preproc_HCPData_PostMinPreproc_TaskRest_BasisHRFModel,
% https://github.com/ColeLab/TaskFCRemoveMeanActivity/blob/master/empiricalfMRIAnalyses/Preproc_HCPData_PostMinPreproc_TaskRest_CanonicalHRFModel_v2.m
% https://github.com/ColeLab/TaskFCRemoveMeanActivity/
% see this paper: https://doi.org/10.1016/j.neuroimage.2018.12.054

if nargin < 3
   upsampleFactor = 1000 ;  
end

%% define basis functions 

%disp('Preparing HRF basis functions for the GLM')
%Define HRF basis functions (obtained from FSL's FLOBS); sampled at 50ms resolution (20 Hz)

samplingRateOfHRF=20;   %In Hz
HRFCanonical = spm_hrf(1/samplingRateOfHRF);


%% interpolate

%disp(['Interpolating design matrix to be at ' num2str(upsampleFactor) 'Hz sampling rate'])
samplingRateInHz=1/tr_seconds;
time = (1:size(taskdesignmat,1))*1/samplingRateInHz;
newtime = time(1):1/upsampleFactor:time(end); % new sampling
new_designmat = interp1(time,taskdesignmat,newtime,'linear');

if size(taskdesignmat,2) == 1
    new_designmat = new_designmat'; 
end

% disp(['Interpolating HRFs to be at ' num2str(upsampleFactor) 'Hz sampling rate'])
samplingRateInHz=samplingRateOfHRF;
time = (1:size(HRFCanonical,1))*1/samplingRateInHz;
newtime = time(1):1/upsampleFactor:time(end); % new sampling
new_HRFCanonical = interp1(time,HRFCanonical,newtime,'linear');

if size(taskdesignmat,2) == 1
    new_HRFCanonical = new_HRFCanonical'; 
end


%% convolve

%disp('Convolving with hemodynamic response function (HRF) basis set')
taskdesignmat_hrf1=zeros(size(new_designmat,1),size(new_designmat,2));

regressorCount=1;
for regressorNum_ForThisBasis=1:size(new_designmat,2)
    convData=conv(new_designmat(:,regressorNum_ForThisBasis),new_HRFCanonical);
    taskdesignmat_hrf1(:,regressorCount)=convData(1:size(new_designmat,1),:);
    regressorCount=regressorCount+1;
end


%% downsamp

% disp('Downsampling designmatrix back to data sampling rate')
taskdesignmat_hrf=downsample(taskdesignmat_hrf1,tr_seconds*upsampleFactor);

end % main func
