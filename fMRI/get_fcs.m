


name_dir = pwd;


PATH_DATA = "aCC24Treg_ts_20p2p7_pkg";
PATH_DATA_REGRESSORS = "aCC24Treg_ts_20p2p7_pkg";

names = string({dir(fullfile(name_dir, PATH_DATA)).name});
names(strcmp(names,'.')) = [];
names(strcmp(names,'..')) = [];
% exclude subjects with mising data
names(strcmp(names,"sub-M87123913")) = []; % no Raven run 3 available
names(strcmp(names,"sub-M87192637")) = []; % no Raven run 2 available
names(strcmp(names,"sub-M87136332")) = []; % no Raven run 3 available
names(strcmp(names,"sub-M87141906")) = []; % Raven run 2 incomplete
names(strcmp(names,"sub-M87174803")) = []; % Raven run 2 incomplete
names(strcmp(names,"sub-M87192557")) = []; % Raven run 3 incomplete


NAMES_RUNS = ["rest","raven_run1","raven_run2", "raven_run3"];
NAMES_RUNS_DATA = ["task-rest_out","task-ravenr01_out","task-ravenr02_out","task-ravenr03_out"];
NAMES_RUNS_REGRESSORS = ["task-rest_desc","task-ravenr01_desc","task-ravenr02_desc","task-ravenr03_desc"];

NAME_FILE_DATA = "schaefer200.ptseries.nii.gz";
NUMNODES = 200;
TR = 2;
S = length(names);
NUMRUNS = 4;
RUNLENGTH = [158, 173, 173, 173];
NUMPERMS = 100;

NUMCONDITIONS = 4;

for s = 1:S
    
    disp(s);
    name_sub = names(s);
    name_sub_short = char(name_sub);
    name_sub_short = string(name_sub_short(end-4:end));
    folder_sub_timings_split = split(name_sub,'-');
    folder_sub_timings = folder_sub_timings_split(2);
   
    onsets_convolved = zeros(sum(RUNLENGTH),NUMCONDITIONS); % 4 conditions = raven_problem,solution,jitter,keypresses
    onsets_combined_concat = zeros(sum(RUNLENGTH),NUMCONDITIONS); % 4 = raven_problem,solution,jitter,keypresses
    data_concat = zeros(sum(RUNLENGTH),NUMNODES);
    
    %% Concatenate runs and nuisance regressors of runs 
    for run = 1:NUMRUNS
        
        time_run_start = 1 + sum(RUNLENGTH(1:run)) - RUNLENGTH(run);
        time_run_end = sum(RUNLENGTH(1:run));

        % load run's data
        % fMRI time series
        file_path_data = fullfile(PATH_DATA,name_sub,NAMES_RUNS_DATA(run),NAME_FILE_DATA);
        data = squeeze(niftiread(file_path_data));
        data_concat(time_run_start:time_run_end,:) = normalize(data);

        % task onsets
        if run == 1
            onsets_combined = zeros(RUNLENGTH(run),NUMCONDITIONS);
        else
            file_name_onsets = append(name_sub,"_",NAMES_RUNS_DATA(run),".csv");
            file_path_onsets = fullfile("onsets", name_sub, file_name_onsets);
            onsets = readtable(file_path_onsets);
            onsets_combined = [onsets.problem, onsets.solution,...
                onsets.jitter, onsets.key_press];   
        end

        onsets_combined_concat(time_run_start:time_run_end,:) = onsets_combined; 

        % get framewise displacement fd
        name_regress_file = append(name_sub, "_", NAMES_RUNS_REGRESSORS(run), "-confounds_timeseries.tsv");
        name_regress_info_file = append(name_sub, "_", NAMES_RUNS_REGRESSORS(run), "-confounds_timeseries.json");
        file_path_regressors = fullfile(PATH_DATA_REGRESSORS,name_sub,'func',name_regress_file);
        file_path_info_regressors = fullfile(PATH_DATA_REGRESSORS,name_sub,'func',name_regress_info_file);
        regressor_struct = tdfread(file_path_regressors);
        
        fd = regressor_struct.framewise_displacement;
        fd = string(fd);
        fd = str2double(fd);
        fd(isnan(fd)) = 0;
        head_motion.(folder_sub_timings).(NAMES_RUNS(run)) = fd; 

        %% Get convolved task onsets for each condition 
        onsets_convolved(time_run_start:time_run_end,:) = apply_CanonicalHRF(onsets_combined, TR); 
        
%         % try without interplotation for comparison only 
%         for i = 1:NUMCONDITIONS
%             cov_hrf = conv(onsets_combined(:,i),spm_hrf(TR));
%             onsets_convolved2(time_run_start:time_run_end,i)=cov_hrf(1:RUNLENGTH(run));
%         end
        
        X_rundesign(time_run_start:time_run_end,run) = 1; 

        
    end
    
    % total timepints of rest and raven trials
    tp_raven = sum(onsets_combined_concat(:,1));
    tp_rest = RUNLENGTH(1);  
    tp_rr = min([tp_raven, tp_rest]);
    tp.(folder_sub_timings).raven = tp_raven;
    tp.(folder_sub_timings).rest = tp_rest;
    tp.(folder_sub_timings).rr = tp_rr;

    % get fcs
    fc_rest = zeros(NUMNODES,NUMNODES,NUMPERMS);
    fc_raven = zeros(NUMNODES,NUMNODES,NUMPERMS);

    for perm = 1:NUMPERMS

        resid_bold = data_concat;
        resid_rest = resid_bold(1:RUNLENGTH(1),:);
        resid_raven = resid_bold(RUNLENGTH(1)+1:end,:);
    
        % get fc of condition specific timepoints
        onsets_convolved_raven = onsets_convolved(RUNLENGTH(1)+1:end,:);
        onsets_raven = onsets_combined_concat(RUNLENGTH(1)+1:end,:);
    
        % rest
        idx_xx = randperm(158);
        idx_xx = sort(idx_xx(1:tp_rr));
        rest_xx = resid_rest(idx_xx,:);
        [fc_rest(:,:,perm), ~] = get_fc(rest_xx, NUMNODES, "fc");
   
        % raven 
        [cut_mask_raven, ~] = get_cut_mask(onsets_raven(:,1), onsets_convolved_raven(:,1));
        idx_cm = find(cut_mask_raven);
        idx_xx = randperm(size(idx_cm,1));
        idx_xx = sort(idx_cm(idx_xx(1:tp_rr)));
        bold_cut = resid_raven(idx_xx,:);
        [fc_raven(:,:,perm),~] = get_fc(bold_cut, NUMNODES, "fc");


    end

    fc.(folder_sub_timings).rest_rr = fc_rest;
    fc.(folder_sub_timings).raven_rr = fc_raven;
 
    
end

clearvars -except fc names head_motion tp 
save('ws_fc.mat', 'fc','head_motion','names','tp', '-v7.3') % -v7.3 in case of big file

function [cut_mask, onsets_cut_mask] = get_cut_mask(onset_vec, onset_vec_conv)
    %%% create a mask (cut_mask) of the n timepoints per block with highest expected bold response
    %%% n is the lenghts of the onsets (not convolved)
    
    cut_mask = false(size(onset_vec));
  
    % lengths of onsets
    if onset_vec(end) == 1 % attach zero at the end to get size of last block 
        onset_vec = [onset_vec; 0];
    end

    if onset_vec(1) == 1 % attach zero at the start to get size of first block
        onset_vec = [0; onset_vec];
    end
    
    block_sizes = find(diff(onset_vec)==-1) - find(diff(onset_vec)==1);
    onsets_cut_mask = zeros(2,length(block_sizes));
    
    % get indexes around expected activation peaks
    [~,idx_peaks] = findpeaks(onset_vec_conv);
    if length(idx_peaks)~=length(block_sizes)
        msg = "peak detection went wrong";
        error(msg)
    end
     
    for numblock = 1:length(block_sizes)

        block_peak = idx_peaks(numblock);
        block_size = block_sizes(numblock);
        
        start_block = block_peak-block_size+1;
        end_block = block_peak+block_size-1;
        block_range = start_block:end_block;
        block_cut = onset_vec_conv(start_block:end_block);
        [~, idx_max_moving_mean_block_cut] = max(movmean(block_cut,block_size)); 
        
        if mod(block_size,2) == 0 % if block_size is even number
            

            start_highest_cluster_block_cut = idx_max_moving_mean_block_cut - (block_size/2);
            end_highest_cluster_block_cut = idx_max_moving_mean_block_cut + (block_size/2) - 1;
        
        else

            start_highest_cluster_block_cut = idx_max_moving_mean_block_cut - ((block_size-1)/2);
            end_highest_cluster_block_cut = idx_max_moving_mean_block_cut + ((block_size-1)/2);

        end
        start_highest_cluster = block_range(start_highest_cluster_block_cut);
        end_highest_cluster = block_range(end_highest_cluster_block_cut);
        cut_mask(start_highest_cluster:end_highest_cluster) = 1;
        onsets_cut_mask(1, numblock) = start_highest_cluster;
        onsets_cut_mask(2, numblock) = end_highest_cluster;
    end

    
end



function [fc, ets] = get_fc(ts, NUMNODES, option)


    if size(ts,1) == NUMNODES

        ts = ts';
        error('check if matrix has format: time x nodes')
    end
    

    if option == "fc"


        fc = corr(ts); % corellating time series of nodes
        fc = (fc+fc')./2; % symmetrize matrix
        fc(1:size(fc,1)+1:end) = 0; % set diagonal elements to zero
        fc = fisherZTransform(fc); % Fisher z-transform all correlations
        ets = 0;
    
    elseif option == "edgefc" % https://github.com/brain-networks/edge-ts/blob/master/main.m
        
        fc = zeros(NUMNODES);
        
        % calculate number of edges
        nedges = NUMNODES*(NUMNODES - 1)/2;

        % indices of unique edges (upper triangle)
        [u,v] = find(triu(ones(NUMNODES),1));
        %idx = (v - 1)*nnodes + u;

        % generate edge time series
        ets = ts(:,u).*ts(:,v);
        
        % mean of ets
        ets_mean = mean(ets);

        % transform to matrix
        fc(triu(ones(NUMNODES),1) == 1) = ets_mean;
        fc = fc + fc'; % symmetrize
        fc = fisherZTransform(fc);

    end

end



