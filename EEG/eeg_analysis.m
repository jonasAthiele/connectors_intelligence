
% Script for preprocessing and connectivity analysis 

% subjects to be excluded
BAD_SUBJECTS = [6,24,86,210]; % excluded
% 6, 24 eeg not completed
% 210 - bad signal in Raven (not in the data set)
% 86 - Bad signal in rest file

% file paths with resting state and raven (rpm) eeg data
path_data = pwd;
FOLDER_DIRECTORY_RAVEN = fullfile(path_data,"data","raven");
FOLDER_DIRECTORY_REST = fullfile(path_data,"data", "rest");
FOLDER_DIRECTORY_REST_PREPROC = fullfile(path_data,"data", "rest_preproc");
FOLDER_DIRECTORY_RAVEN_PREPROC = fullfile(path_data,"data", "raven_preproc");
FOLDER_DIRECTORY_FC = fullfile(path_data,"data", "fc_eeg");

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
        
        %% Preprocessing
        % preprocess rest data as described in https://doi.org/10.1016/j.intell.2023.101780
        
        file_name = append("rest1_", string(id), ".bdf");
        rest_file_path = fullfile(FOLDER_DIRECTORY_REST, file_name);
        [EEG, stats] = preprocess_rest(rest_file_path);
        
        save_name = append("rest1_preprocessed_", string(id), ".set");
        pop_saveset(EEG, 'filename',char(save_name),'filepath',char(FOLDER_DIRECTORY_REST_PREPROC));
         

        file_name = append("rest1_preprocessed_", string(id), ".set");
        rest_path = char(fullfile(FOLDER_DIRECTORY_REST_PREPROC, file_name));
        EEG_rest = pop_loadset(rest_path);
        
        % read rpm data
        try
            file_name = append(string(id), "_outcome.set");
            raven_file_path = char(fullfile(FOLDER_DIRECTORY_RAVEN, file_name));
            EEG_raven = pop_loadset(raven_file_path);
        catch
            file_name = append(string(id), "_a_outcome.set");
            raven_file_path = char(fullfile(FOLDER_DIRECTORY_RAVEN, file_name));
            EEG_raven = pop_loadset(raven_file_path);
        end

        n_epochs = EEG_raven.trials;
        epoch_info = string({EEG_raven.event.type});

        % rereference to average
        EEG_raven = pop_reref(EEG_raven); % average reference


        % sanity check
        if n_epochs~=length(epoch_info)
            error('trials unequal epochs')
        end
        
        % test for high amplitudes in rpm data
        amplitudeThreshold = 160;
        bad_channels = false(64, n_epochs);
        for e = 1:n_epochs 
            e_data = EEG_raven.data(:,:,e);
            amplitudeMax = max(abs(e_data),[],2);
            bad_channels(:,e) = (amplitudeMax > amplitudeThreshold);     
        end

        
        if any(bad_channels(:))
            [EEG_raven, ~, badlist] = pop_TBT(EEG_raven, bad_channels, 10, 1, 0);
            error('raven bad channels, why?')
        end

        % save eeg sets for rpm data
        eeg_raven_temp = EEG_raven;
        save_name = append(string(id), "_raven.set");
        pop_saveset(eeg_raven_temp, 'filename',char(save_name),'filepath',char(FOLDER_DIRECTORY_RAVEN_PREPROC));

        %% Connectivity
        % reload rest file for data analysis 
        file_name = append("rest1_preprocessed_", string(id), ".set");
        rest_path = char(fullfile(FOLDER_DIRECTORY_REST_PREPROC, file_name));
        brain_measures.rest = get_brain_measures(rest_path, false);

        % reload rpm file for data analysis 
        file_name = append(string(id), "_raven.set");
        file_path = char(fullfile(FOLDER_DIRECTORY_RAVEN_PREPROC, file_name));
        EEG_raven = [];
        if isfile(file_path)
            EEG_raven = pop_loadset(file_path);
            brain_measures.raven = get_brain_measures(file_path, true);
            save_name = append(string(id), '_brain_measures.mat');
            file_path = fullfile(FOLDER_DIRECTORY_FC, save_name);
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
    EEG = pop_reref(EEG, [],'exclude',[65:72]); %average reference
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

% function for computing connectivity
function brain_measures = get_brain_measures(file_path, filter_alpha)


    cfg = []; cfg.dataset = file_path;  
    cfg.channel = {'EEG'};
    data = ft_preprocessing(cfg);
    srate = data.fsample;

    % dynamic frequency analysis for filtering out high occipital alpha periods
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.foi          = 9:0.5:11;                         % analysis 9 to 11 Hz in steps of 0.5 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.25;   % length of time window = 0.25 sec
    cfg.toi = 'all';
    TFRhann = ft_freqanalysis(cfg, data);

    % get time points without NANs (introduced by dynamic freq analysis)
    test_data = squeeze(TFRhann.powspctrm(1,1,1,:));
    idx_data = ~isnan(test_data);
    
    n_trial = size(TFRhann.powspctrm,1);
    theta_range = [4,7];
    theta_freqs = theta_range(1):0.5:theta_range(2);
    connMatrix= zeros(n_trial,64,64,size(theta_freqs,2));
    fc_epochs = zeros(n_trial,1);
    
    for trial = 1:n_trial
        
        disp(trial)

        %% Identify high alpha epochs
        % power spectrum over time of the trial
        dat_temp = squeeze(TFRhann.powspctrm(trial,[26,63],:,idx_data)); % electrode 26 = P03, electrode 63 = P04 (alpha has been observed to get high at this electrodes when trial was solved)   
        dat_temp = squeeze(mean(dat_temp,2));
        dat_temp = squeeze(mean(dat_temp));
        % smooth spectrum with movemean
        dat_temp_movmean = movmean(dat_temp,3*srate);
        
        % get mean of lowest 30% of alpha power across trial as threshold
        aa_dat = sort(dat_temp_movmean);
        alpha_thresh = mean(aa_dat(1:int16(size(aa_dat,2)*0.3)));
        
        % get standard deviation of alpha power across trial
        SD_alpha = std(dat_temp_movmean);
    
        % timepoints of alpha is lower than threshold + 2 standard
        % deviations of whole trial alpha
        idx_alpha_valid_temp = dat_temp_movmean < alpha_thresh + 2*SD_alpha;
        
        frames_idx_alpha = find(idx_alpha_valid_temp==1);
        idx_alpha_start = frames_idx_alpha(1);  % start of low alpha
        T_alpha = size(idx_alpha_valid_temp,2); % end of trial

        % get 2 second epochs in which no more than 10% of alpha is above thresh 
        number_alpha_epochs = floor((T_alpha - idx_alpha_start) / (2*srate)); % 2 second blocks
        tt_start = idx_alpha_start:(2*srate):idx_alpha_start+(2*srate)*number_alpha_epochs;
        tt_end = (idx_alpha_start+(2*srate):(2*srate):T_alpha)-1;

        idx_alpha_valid = false(size(idx_alpha_valid_temp));
        for b = 1:number_alpha_epochs
            idx_alpha_block = idx_alpha_valid_temp(tt_start(b):tt_end(b));
            not_valids = sum(~idx_alpha_block);
            percentage_not_valid = not_valids/(2*srate);
            if percentage_not_valid < 0.1 % 10 percent not valids allowed
                idx_alpha_valid(tt_start(b):tt_end(b)) = true;
            end
        end
          
        cfg = []; 
        cfg.trials = trial;
        data_trial = ft_preprocessing(cfg, data);
        ts_temp = cell2mat(data_trial.trial);
        if filter_alpha
            % only use low alpha epochs
            ts_temp = ts_temp(:,idx_data);
            ts_temp = ts_temp(:,idx_alpha_valid);
        end
        data_trial.trial{1,1} = ts_temp;
        T_trial = size(ts_temp,2);
        data_trial.sampleinfo = [1, T_trial];
        data_trial.hdr.nSamples = T_trial;
        data_trial.hdr.nTrials = 1;
        time_temp = cell2mat(data_trial.time);
        data_trial.time{1,1} = time_temp(1:T_trial);
        
        % epoch in 2 second epochs
        cfg = [];
        cfg.length = 2; % 2 second blocks
        cfg.overlap = 0;
        data_trial_epoched = ft_redefinetrial(cfg, data_trial);

        % filter data in theta range
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = theta_range;
        data_trial_epoched = ft_preprocessing(cfg, data_trial_epoched);
        
        
        %% Connectivity 
        % Fourier components
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.taper = 'dpss';
        cfg.output = 'fourier';
        cfg.keeptrials = 'yes';
        cfg.pad = 'nextpow2';
        cfg.foi = theta_freqs;
        cfg.tapsmofrq = 1;
        Freq = ft_freqanalysis(cfg, data_trial_epoched);
        
        % connectivity
        cfg = [];
        cfg.method = 'wpli_debiased';
        source_conn = ft_connectivityanalysis(cfg, Freq);
        clear Freq;
        
        % average across frequency bins
        connMatrix_trial = abs(source_conn.wpli_debiasedspctrm);
        if all(all(isnan(connMatrix_trial)))
            error('Connectivity matrix only contains NaN');
        end
        %connMatrix_trial(1:size(connMatrix_trial,1)+1:end) = 0; % set diagonal elements to zero
        connMatrix(trial,:,:,:) = connMatrix_trial;
        fc_epochs(trial) = size(data_trial_epoched.trial,2);


           
    end
    
    brain_measures.freq = theta_freqs;
    brain_measures.fcdwpli = connMatrix;
    brain_measures.fc_epochs = fc_epochs;


end