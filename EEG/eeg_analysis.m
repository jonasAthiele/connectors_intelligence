
% Script for eeg preprocessing and mse calculations

% subjects to be excluded
BAD_SUBJECTS = [6,24,86,210]; % excluded
% 6, 24 eeg not completed
% 210 - bad signal in Raven (not in the data set)
% 86 - bad signal in rest file

% file paths with resting state and raven (rpm) eeg data
path_data = pwd;
FOLDER_DIRECTORY_RAVEN = fullfile(path_data,"data","raven");
FOLDER_DIRECTORY_REST = fullfile(path_data,"data", "rest");
FOLDER_DIRECTORY_REST_PREPROC = fullfile(path_data,"data", "rest_preproc");
FOLDER_DIRECTORY_RAVEN_PREPROC = fullfile(path_data,"data", "raven_preproc");
FOLDER_DIRECTORY_MSE = fullfile(path_data,"data", "eeg_mse_2");

% behavioral data
beh_file = fullfile(path_data,"descr_behav_data.csv");
BEHAVIORAL_DATA = readtable(beh_file);

SUBJECT_IDS = BEHAVIORAL_DATA.ID;
n_subjects = length(SUBJECT_IDS);

eeglab;
close

for s = 1:n_subjects

    id = SUBJECT_IDS(s);

    if ~ismember(id, BAD_SUBJECTS)

        brain_measures = [];

        % raven (rpm) data is already preprocessed as described in https://doi.org/10.1016/j.intell.2023.101780
        % raven (rpm) data available at https://osf.io/htrsg/
        % rest data is available in raw form at: https://osf.io/kv2sx/
        % data belongs to Adam Chuderski, Michał Ociepka
        
        %% Preprocessing of rest data
        % preprocess rest data as described in https://doi.org/10.1016/j.intell.2023.101780
        
        file_name = append("rest1_", string(id), ".bdf");
        rest_file_path = fullfile(FOLDER_DIRECTORY_REST, file_name);
        [EEG, stats] = preprocess_rest(rest_file_path);
        
        save_name = append("rest1_preprocessed_", string(id), ".set");
        pop_saveset(EEG, 'filename',char(save_name),'filepath',char(FOLDER_DIRECTORY_REST_PREPROC));
         

        %% MSE
        % load rest file for data analysis 
        file_name = append("rest1_preprocessed_", string(id), ".set");
        rest_path = char(fullfile(FOLDER_DIRECTORY_REST_PREPROC, file_name));
        brain_measures.rest = get_brain_measures(rest_path);

        % load rpm file for data analysis 
        file_name = append(string(id), "_outcome.set");
        if isfile(char(fullfile(FOLDER_DIRECTORY_RAVEN_PREPROC, file_name)))
            file_path = char(fullfile(FOLDER_DIRECTORY_RAVEN_PREPROC, file_name));
        else
            file_name = append(string(id), "_a_outcome.set");
            file_path = char(fullfile(FOLDER_DIRECTORY_RAVEN_PREPROC, file_name));
        end
        
        if isfile(file_path)
            brain_measures.raven = get_brain_measures(file_path);
            save_name = append(string(id), '_brain_measures.mat');
            file_path = fullfile(FOLDER_DIRECTORY_MSE, save_name);
            save(file_path,'brain_measures');
        end
      

    end

end
  

% function for preprocessing rest data
function [EEG_preprocessed, stats] = preprocess_rest(eeg_file_path)
    
    % load file
    EEG = pop_biosig(char(eeg_file_path)); 
    set_name = split(eeg_file_path,'\');
    set_name = set_name(end);
    
    
    % check if channels are ordered correctly (needs to be adapted for other nets)
    label_template = ["Fp1"	"AF7"	"AF3"	"F1"	"F3"	"F5"	"F7"...
    	"FT7"	"FC5"	"FC3"	"FC1"	"C1"	"C3"	"C5"	"T7"...
    	"TP7"	"CP5"	"CP3"	"CP1"	"P1"	"P3"	"P5"	"P7"...
    	"P9"	"PO7"	"PO3"	"O1"	"Iz"	"Oz"	"POz"	"Pz"...
    	"CPz"	"Fpz"	"Fp2"	"AF8"	"AF4"	"AFz"	"Fz"	"F2"...
    	"F4"	"F6"	"F8"	"FT8"	"FC6"	"FC4"	"FC2"	"FCz"...
        "Cz"	"C2"	"C4"	"C6"	"T8"	"TP8"	"CP6"	"CP4"...
    	"CP2"	"P2"	"P4"	"P6"	"P8"	"P10"	"PO8"	"PO4"	"O2"];
    
    labels = string({EEG.chanlocs.labels});
    if ~all(labels(1:64) == label_template(1:64))
        error('wrong channel order!')
    end
   
    % set channel locations
    EEG = pop_chanedit(EEG, 'settype',{'65:68','EOG'},'settype',{'69:72','MISC'},'settype',{'1:64','EEG'});
    EEG = pop_chanedit(EEG, 'lookup','C:\\Program Files\\MATLAB\\R2021b\\toolbox\\eeglab_current\\eeglab2023.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
    EEG = pop_chanedit(EEG, 'load',[],'lookup','.\channel_locations_64.ced');
    EEG = pop_reref(EEG, [71 72]); %linked mastiods
    EEG = pop_chanedit(EEG, 'settype',{'1:64','EEG'});
    EEG = pop_chanedit(EEG, 'settype',{'65:68','EOG'});
    EEG = pop_chanedit(EEG, 'settype',{'69:72','MISC'});


    %EEG = pop_rmbase( EEG, [],[]); %baseline correction
    %pop_eegplot( EEG, 1, 1, 1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',1,'plotfreqz',1);
    %pop_eegplot( EEG, 1, 1, 1);
    EEG = pop_eegfiltnew(EEG, 'locutoff',49.5,'hicutoff',50.5,'revfilt',1,'plotfreqz',1);
    close all
    %figure; pop_spectopo(EEG, 1, [0      379996.0938], 'EEG' , 'percent', 50, 'freq', [6 10 22], 'freqrange',[0 80],'electrodes','off');
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    EEG = pop_crls_regression(EEG, [65  66  67  68], 3, 0.9999, 0.01);
    %pop_eegplot( EEG, 1, 1, 1);
    
    
    t_events = [EEG.event.latency];
    t_start = t_events(end-1);
    t_end = t_events(end);

    t_range = floor((t_end-t_start)/EEG.srate);
    if t_range ~= 360
        error('t start strange')
    end
    
    EEG = pop_select(EEG, 'point',[t_start+1 t_end-1]);
    
    % epoching
    EpochLength = 10;
    EpochOverlap = 0;

    % create spatially distributed markers
    EEG = eeg_regepochs(EEG,'recurrence',EpochLength * (1-EpochOverlap),'eventtype','epoch_start','extractepochs','off');
    
    % epoch
    EEG = pop_epoch(EEG,{'epoch_start'},[0 EpochLength],'epochinfo','yes');
    
    % drop non EEG channels
    EEG = pop_select(EEG, 'rmchantype',{'MISC','EOG'});
    
    % check epochs for peak to peak amplitudes
    % identify which epochs of each channel have amplitudes above threshold
    amplitudeThreshold = 160;
    n_epochs = EEG.trials;
    bad_channels = false(64, n_epochs);
    for e = 1:n_epochs
        
        e_data = EEG.data(:,:,e);
        amplitudeMax = max(abs(e_data),[],2);
        bad_channels(:,e) = (amplitudeMax > amplitudeThreshold);
            
    end

    % remove or interpolate bad epochs
    [EEG, ~, badlist] = pop_TBT(EEG, bad_channels, 10, 1, 0); %https://doi.org/10.5281/zenodo.1241518
    
    stats.bad_epochs = badlist.nbadtrial;
    EEG.setname=set_name;
    EEG_preprocessed = EEG;
    
    
end

% function for computing MSE
function brain_measures = get_brain_measures(file_path)


    cfg = []; 
    cfg.dataset = file_path;  
    cfg.channel = {'EEG'};
    data = ft_preprocessing(cfg);
    n_trial = size(data.trial,2);
    n_elec = length(data.label);
    
    mse = []; % initializing
    n_epochs = nan(n_trial,1); % initializing
    for trial = 1:n_trial
        
        %disp(trial)
        cfg = []; 
        cfg.trials = trial;
        data_trial = ft_preprocessing(cfg, data);
       
        % epoch in 10 second epochs
        cfg = [];
        cfg.length = 10; % 10 second blocks
        cfg.overlap = 0;
        data_trial_epoched = ft_redefinetrial(cfg, data_trial);

        % filter data in freq range
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [1,40];
        data_trial_filt = ft_preprocessing(cfg, data_trial_epoched);


        n_epochs_trial = size(data_trial_filt.trial,2);  
        
        if n_epochs_trial < 1
            error("no epochs, check epoching")
        end
        
        % compute MSE
        mse_trial = zeros(n_epochs_trial, n_elec, 20);
        for epoch = 1:n_epochs_trial
            
            cfg = []; 
            cfg.trials = epoch;
            data_epoch = ft_preprocessing(cfg, data_trial_filt);
            ts = cell2mat(data_epoch.trial);
            m = 2;
            r = 0.2;
            for tau = 1:20
                for elec = 1:n_elec
                    x = ts(elec,:);
                    [mse_trial(epoch,elec,tau), ~, ~] = multiscaleSampleEntropy(x, m, r, tau);
                end
    
            end
        end

        mse(trial,:,:,:) = mse_trial;
        n_epochs(trial) = n_epochs_trial;
           
    end
    
    brain_measures.mse = mse;
    brain_measures.n_epochs = n_epochs;


end
