
% Script for computing centrality measures and their relationship with rpm scores

% behavioral data
path_data = pwd;
beh_file = fullfile(path_data,"descr_behav_data.csv");
BEHAVIORAL_DATA = readtable(beh_file);

SUBJECT_IDS = BEHAVIORAL_DATA.ID;

FOLDER_DIRECTORY_FC = fullfile(path_data,"data","fc_eeg");

% set up for plots
eeglab
EEG = pop_loadset('1_outcome.set');
setup_plots();

NUMPERMS = 100; % number of permutations
threshold_fc = 0.5; % threshold for thresholding FC
rescale_meth = 'rescale'; % rescale FC to be in the range 0 to 1


% analyses
file_names={dir(FOLDER_DIRECTORY_FC).name};
file_names=file_names(~ismember(file_names,{'.','..'}));

for s = 1:length(file_names)
    
    f = file_names(s);
    fsplit = split(f,'_');
    id(s) = str2double(string(fsplit(1)));
end
id = sort(id);
% exclude subjects due to missing data
id(id == 72) = []; % no age
id(id == 140) = []; % no age


n_subjects = length(id); % number subjects


%% Get FCs of all subjects to make a group average FC
for s = 1:n_subjects

    disp(s)
    id_s = id(s);

    file_name = append(string(id_s), '_brain_measures.mat');
    file_path = fullfile(FOLDER_DIRECTORY_FC, file_name);
    
    load(file_path)
        
    rest_epochs = brain_measures.rest.fc_epochs;
    raven_epochs = brain_measures.raven.fc_epochs;

    rest_epochs_sum = sum(rest_epochs);
    raven_epochs_sum = sum(raven_epochs);

    n_trials_raven = size(raven_epochs,1);
    n_trials_rest = size(rest_epochs,1);

    if rest_epochs_sum < raven_epochs_sum
        
        perm_idx = randperm(n_trials_raven);
        cum_sum_epochs = cumsum(raven_epochs(perm_idx));
        [~,idx_end_raven] = min(abs(cum_sum_epochs - rest_epochs_sum));
        idx_raven = perm_idx(1:idx_end_raven);
        idx_rest = 1:n_trials_rest;

    elseif rest_epochs_sum > raven_epochs_sum
        
        perm_idx = randperm(n_trials_rest);
        cum_sum_epochs = cumsum(rest_epochs(perm_idx));
        [~,idx_end_rest] = min(abs(cum_sum_epochs - raven_epochs_sum));
        idx_rest = perm_idx(1:idx_end_rest);
        idx_raven = 1:n_trials_raven;

    else
        idx_rest = 1:n_trials_rest;
        idx_raven = 1:n_trials_raven;
    end

    fc_raven(:,:,:,s) = squeeze(mean(brain_measures.raven.fcdwpli(idx_raven,:,:,:)));
    fc_rest(:,:,:,s) = squeeze(mean(brain_measures.rest.fcdwpli(idx_rest,:,:,:)));
   
end

group_fc_raven = squeeze(mean(fc_raven,4));
group_fc_rest = squeeze(mean(fc_rest,4));

names_freqs = brain_measures.raven.freq;
n_freqs = size(names_freqs,2);

%% Compute partitions on group average FC of each frequency
for f = 1:n_freqs
    
    fcf = squeeze(group_fc_raven(:,:,f));
    fcf(1:size(fcf,1)+1:end) = 0; % set diagonal elements to zero 
    fcf = threshold_proportional(fcf,threshold_fc); % threshold fc
    fcf = rescale_fc(fcf,rescale_meth); % rescale in range 0 to 1
    
    [Q_raven, CI_raven, CI2_raven] = get_communities(fcf, 2); % gamma = 2
   
    C_raven_group(:,f) = CI2_raven;


    fcf = squeeze(group_fc_rest(:,:,f));
    fcf(1:size(fcf,1)+1:end) = 0; % set diagonal elements to zero
    fcf = threshold_proportional(fcf,threshold_fc);
    fcf = rescale_fc(fcf,rescale_meth); % rescale in range 0 to 1
    
    [Q_raven, CI_raven, CI2_raven] = get_communities(fcf, 2); % gamma = 2
   
    C_rest_group(:,f) = CI2_raven;


end

%% Compute centrality measures (degree and participation coefficient)
for s = 1:n_subjects
    
    disp(s)
    id_s = id(s);

    file_name = append(string(id_s), '_brain_measures.mat');
    file_path = fullfile(FOLDER_DIRECTORY_FC, file_name);
    
    load(file_path)


    for p = 1:NUMPERMS
        disp(s)
        disp(p)
        
        rest_epochs = brain_measures.rest.fc_epochs;
        raven_epochs = brain_measures.raven.fc_epochs;

        rest_epochs_sum = sum(rest_epochs);
        raven_epochs_sum = sum(raven_epochs);

        n_trials_raven = size(raven_epochs,1);
        n_trials_rest = size(rest_epochs,1);

        trials_raven(s) = n_trials_raven;
        epochs_raven(s) = raven_epochs_sum;

        if rest_epochs_sum < raven_epochs_sum
            
            perm_idx = randperm(n_trials_raven);
            cum_sum_epochs = cumsum(raven_epochs(perm_idx));
            [~,idx_end_raven] = min(abs(cum_sum_epochs - rest_epochs_sum));
            idx_raven = perm_idx(1:idx_end_raven);
            idx_rest = 1:n_trials_rest;

        elseif rest_epochs_sum > raven_epochs_sum
            
            perm_idx = randperm(n_trials_rest);
            cum_sum_epochs = cumsum(rest_epochs(perm_idx));
            [~,idx_end_rest] = min(abs(cum_sum_epochs - raven_epochs_sum));
            idx_rest = perm_idx(1:idx_end_rest);
            idx_raven = 1:n_trials_raven;

        else
            idx_rest = 1:n_trials_rest;
            idx_raven = 1:n_trials_raven;

        end


        fc = squeeze(mean(brain_measures.raven.fcdwpli(idx_raven,:,:,:)));
        for f = 1:n_freqs
            
            fcf = squeeze(fc(:,:,f));
            fcf(1:size(fcf,1)+1:end) = 0; % set diagonal elements to zero
            fcf = threshold_proportional(fcf,threshold_fc);
            fcf = rescale_fc(fcf,rescale_meth); % rescale in range 0 to 1
    
            degree_raven_p(:,f,p) = mean(fcf);
            
            W = fcf;
            C = squeeze(C_rest_group(:,f));
            partiC_raven_p(:,f,p) = participation_coef(W,C);
        
        end
    
    
        fc = squeeze(mean(brain_measures.rest.fcdwpli(idx_rest,:,:,:)));
        for f = 1:n_freqs
    
            fcf = squeeze(fc(:,:,f));
            fcf(1:size(fcf,1)+1:end) = 0; % set diagonal elements to zero
            fcf = threshold_proportional(fcf,threshold_fc);
            fcf = rescale_fc(fcf,rescale_meth); % rescale in range 0 to 1
            
            degree_rest_p(:,f,p) = mean(fcf);
            W = fcf;
            C = squeeze(C_rest_group(:,f));
            partiC_rest_p(:,f,p) = participation_coef(W,C);

    
        end



    end

    degree_rest(:,:,s) = squeeze(mean(degree_rest_p,3));
    degree_raven(:,:,s) = squeeze(mean(degree_raven_p,3));

    partiC_rest(:,:,s) = squeeze(mean(partiC_rest_p,3));
    partiC_raven(:,:,s) = squeeze(mean(partiC_raven_p,3));


    idx_subject = find(BEHAVIORAL_DATA.ID == id_s);
    
    age(s) = BEHAVIORAL_DATA.WIEK(idx_subject);
    sex_temp = BEHAVIORAL_DATA.PLEC(idx_subject);
    if sex_temp == "M" % male
        sex_translate = 0;
    elseif sex_temp == "K" % female
        sex_translate = 1;
    end
    sex(s) = sex_translate;
    raven_score(s) = BEHAVIORAL_DATA.RAVEN(idx_subject);
   


end

% test which subjetcs have enough rpm trials
idx_quest = epochs_raven < 180;
idx_valid = ~idx_quest; % 1 if subject has at least 180 2 second trials for rpm items (same number as resting state)

C = [sex(idx_valid)', age(idx_valid)', epochs_raven(idx_valid)']; % confounds
y = raven_score(idx_valid)'; % rpm sum score


electrode_names = string({EEG.chanlocs.labels}); % get electrode names
names_freqs = brain_measures.raven.freq; % get names of evaluated frequencies

%% Plot associations with intelligence (rpm sum score)

% participation coefficient of one frequency
plot_freq = 7;
idx_plot_freq = find(names_freqs == plot_freq);
pc_rest = squeeze(partiC_rest(:,idx_plot_freq,idx_valid));
pc_raven = squeeze(partiC_raven(:,idx_plot_freq,idx_valid));

X = pc_raven - pc_rest;
[rho,p,p_adj,num_comparisons] = get_relation(X',y,C, "fdr");
cMap = flipud(slanCM('RdBu'));
figure()
topoplot(rho,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-0.3, 0.3]);
%cMap = slanCM('binary');
pmask = zeros(64,1);
pmask(p<0.05)=0.6;
figure()
topoplot(pmask,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-1, 1]);

% participation coefficient for whole theta range
pc_rest = squeeze(mean(partiC_rest(:,1:7,idx_valid),2));
pc_raven = squeeze(mean(partiC_raven(:,1:7,idx_valid),2));
X = pc_raven - pc_rest;
[rho,p,p_adj,num_comparisons] = get_relation(X',y,C, "fdr");
cMap = flipud(slanCM('RdBu'));
figure()
topoplot(rho,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-0.3, 0.3]);
%cMap = slanCM('binary');
pmask = zeros(64,1);
pmask(p<0.05)=0.6;
figure()
topoplot(pmask,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-1, 1]);

% table participation coeff p < 0.05 (not corrected)
names_sig = electrode_names(p<0.05);
p_sig = p(p<0.05);
rho_sig = rho(p<0.05);
table_res_parti = array2table([names_sig', rho_sig', p_sig'],...
    'VariableNames',["region","rho","p_uncorrected"]);
writetable(table_res_parti, "results_participation_uncorrected.csv")


% degree of one frequency (frequency defined above)
deg_rest = squeeze(degree_rest(:,idx_plot_freq,idx_valid));
deg_raven = squeeze(degree_raven(:,idx_plot_freq,idx_valid));

X = deg_raven - deg_rest;
[rho,p,p_adj,num_comparisons] = get_relation(X',y,C, "fdr");
cMap = flipud(slanCM('RdBu'));
figure()
topoplot(rho,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-0.3, 0.3]);
pmask = zeros(64,1);
pmask(p<0.05)=0.6;
figure()
topoplot(pmask,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-1, 1]);

% degree for whole theta range
deg_rest = squeeze(mean(degree_rest(:,1:7,idx_valid),2));
deg_raven = squeeze(mean(degree_raven(:,1:7,idx_valid),2));
X = deg_raven - deg_rest;
[rho,p,p_adj,num_comparisons] = get_relation(X',y,C, "fdr");
cMap = flipud(slanCM('RdBu'));
figure()
topoplot(rho,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-0.3, 0.3]);
pmask = zeros(64,1);
pmask(p<0.05)=0.6;
figure()
topoplot(pmask,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-1, 1]);

% table degree p < 0.05 (not corrected)
names_sig = electrode_names(p<0.05);
p_sig = p(p<0.05);
rho_sig = rho(p<0.05);
table_res_degree = array2table([names_sig', rho_sig', p_sig'],...
    'VariableNames',["region","rho","p_uncorrected"]);
writetable(table_res_degree, "results_degree_uncorrected.csv")


%% Plots for Supplement
% plot histogram of raven scores
figure()
histogram(y,10, 'FaceColor',[0.4 0.4 0.4])
box off
set(gca,'fontsize', 16) 

% plot fcs
% raven
fci = squeeze(mean(fc_raven(:,:,1:7,idx_valid),3));
fcf = squeeze(mean(fci,3));
value_background = nan;
fcf(1:size(fcf,1)+1:end) = value_background; % set diagonal elements to zero
mask = ones(64)*value_background;
mask(triu(ones(64))==1) = fcf(triu(ones(64))==1);
f = figure;
h=imagesc(mask);
f.Position = [100 100 900 800];
ax = gca;
ax.TickLength = [0 0];
xticks(1:64);
yticks([1:64]);
yticklabels(electrode_names)
xticklabels(electrode_names)
set(h, 'AlphaData', ~isnan(mask));
cMap = slanCM('Reds');
colormap(cMap)  
bb = colorbar;

%% rest
fci = squeeze(mean(fc_rest(:,:,1:7,idx_valid),3));
fcf = squeeze(mean(fci,3));
value_background = nan;
fcf(1:size(fcf,1)+1:end) = value_background; % set diagonal elements to zero
mask = ones(64)*value_background;
mask(triu(ones(64))==1) = fcf(triu(ones(64))==1);
f = figure;
h=imagesc(mask);
f.Position = [100 100 900 800];
ax = gca;
ax.TickLength = [0 0];
xticks(1:64);
yticks([1:64]);
yticklabels(electrode_names)
xticklabels(electrode_names)
set(h, 'AlphaData', ~isnan(mask));
cMap = slanCM('Reds');
colormap(cMap)  
bb = colorbar;


%% raven - rest
fci_rest = squeeze(mean(fc_rest(:,:,1:7,idx_valid),3));
fci_raven = squeeze(mean(fc_raven(:,:,1:7,idx_valid),3));

fcf = fci_raven-fci_rest;
fcf = squeeze(mean(fcf,3));
value_background = nan;
fcf(1:size(fcf,1)+1:end) = value_background; % set diagonal elements to zero
mask = ones(64)*value_background;
mask(triu(ones(64))==1) = fcf(triu(ones(64))==1);

f = figure;
f.Position = [100 100 900 800];
cMap = flipud(slanCM('Blues'));
colormap(cMap)  
h=imagesc(mask);
ax = gca;
ax.TickLength = [0 0];
xticks(1:64);
yticks([1:64]);
yticklabels(electrode_names)
xticklabels(electrode_names)
set(h, 'AlphaData', ~isnan(mask));
bb = colorbar;


% PC rest
pc = squeeze(mean(partiC_rest(:,1:7,idx_valid),2));
pc = squeeze(mean(pc,2));
pc = normalize(pc,'range',[0,1]);
cMap = slanCM('Reds');
figure()
topoplot(pc,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[0 1]);
colorbar

% PC raven
pc = squeeze(mean(partiC_raven(:,1:7,idx_valid),2));
pc = squeeze(mean(pc,2));
pc = normalize(pc,'range',[0,1]);
cMap = slanCM('Reds');
figure()
topoplot(pc,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[0 1]);
colorbar

% PC raven - rest
pc_raven = squeeze(mean(partiC_raven(:,1:7,idx_valid),2));
pc_rest = squeeze(mean(partiC_rest(:,1:7,idx_valid),2));

for s = 1:size(pc_raven,2)
    pcs(:,s) = pc_raven(:,s) - pc_rest(:,s); 
end

pc = squeeze(mean(pcs,2));
%pc = normalize(pi,'range',[0,1]);
cMap = flipud(slanCM('RdBu'));
figure()
topoplot(pc,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-0.01 0.01]);
colorbar


% degree rest
deg = squeeze(mean(degree_rest(:,1:7,idx_valid),2));
deg = squeeze(mean(deg,2));
deg = normalize(deg,'range',[0,1]);
cMap = slanCM('Reds');
figure()
topoplot(deg,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[0 1]);
colorbar

% degree raven
deg = squeeze(mean(degree_raven(:,1:7,idx_valid),2));
deg = squeeze(mean(deg,2));
deg = normalize(deg,'range',[0,1]);
cMap = slanCM('Reds');
figure()
topoplot(deg,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[0 1]);
colorbar

% degree raven - rest
deg_raven = squeeze(mean(degree_raven(:,1:7,idx_valid),2));
deg_rest = squeeze(mean(degree_rest(:,1:7,idx_valid),2));

for s = 1:size(deg_raven,2)
    degs(:,s) = deg_raven(:,s) - deg_rest(:,s); 
end

deg = squeeze(mean(degs,2));
%deg = normalize(pi,'range',[0,1]);
cMap = flipud(slanCM('Blues'));
figure()
topoplot(deg,EEG.chanlocs,'electrodes','labels', 'colormap', cMap, 'interplimits', 'electrodes','maplimits',[-0.1 -0.06]);
colorbar


%% Some information to store

varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd"];

% total sample:
no_subjects = length(idx_valid);
age_min = min(age);
age_max = max(age);
age_mean = mean(age);
age_sd = std(age);
no_males = sum(sex==0);
no_females = sum(sex==1);

rpm = raven_score;
rpm_min = min(rpm);
rpm_max = max(rpm);
rpm_mean = mean(rpm);
rpm_sd = std(rpm);


data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd];

info_all_subs = array2table(data,'VariableNames',varNames);
writetable(info_all_subs, "info_subjects_all.csv")

% included subjects
varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd",...
    "r_rpm_age", "p_rpm_age","r_rpm_sex","p_rpm_sex","r_rpm_epochs",...
    "p_rpm_epochs"];

no_subjects = sum(idx_valid==1);
age_min = min(age(idx_valid));
age_max = max(age(idx_valid));
age_mean = mean(age(idx_valid));
age_sd = std(age(idx_valid));
no_males = sum(sex(idx_valid)==0);
no_females = sum(sex(idx_valid)==1);

rpm = raven_score(idx_valid);
rpm_min = min(rpm);
rpm_max = max(rpm);
rpm_mean = mean(rpm);
rpm_sd = std(rpm);

[r_rpm_age, p_rpm_age] = corr(raven_score(idx_valid)',age(idx_valid)','type','Spearman'); % correlation between age and RPM scores
[r_rpm_sex, p_rpm_sex] = corr(raven_score(idx_valid)',sex(idx_valid)','type','Spearman');
[r_rpm_epochs, p_rpm_epochs] = corr(raven_score(idx_valid)',epochs_raven(idx_valid)','type','Spearman');


data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd, ...
    r_rpm_age,p_rpm_age,r_rpm_sex,p_rpm_sex,-r_rpm_epochs,p_rpm_epochs]; 
% minus at correlation with number of epochs to get correlation between rpm and removed epochs to be consistent with fMRI analysis

info_subs = array2table(data,'VariableNames',varNames);

writetable(info_subs, "info_subjects.csv")



function W_nrm = rescale_fc(W, option)
        
        if option == "pos"
            
            W = W.*~eye(size(W));  % set diagonale zero
            W(W<0) = 0;
            W_nrm = W./max(abs(W(:)));

        elseif option == "abs"
            
            W = W.*~eye(size(W));  % set diagonale zero
            W = abs(W);
            W_nrm = W./max(abs(W(:)));

        elseif option == "rescale"
            
            % normalize FC
            W_min = min(min(W));
            W_max = max(max(W));
            W_nrm = (W - W_min) / (W_max - W_min);
            W_nrm = W_nrm.*~eye(size(W_nrm)); % set diagonale zero
        end

end



function [Q, CI, CI2] = get_communities(FC, gam)

    [N,~] = size(FC); % N-number nodes
    
    R = 100; % runs of Louvain
    
    % parameters for calculating consensus between runs of Louvain
    tau = 0.1; % threshold which controls the resolution of the reclustering
    reps = 10; % number of times that the clustering algorithm is reapplied
    
    G = length(gam); % gamma influences number of clusters/communities (higher gamma = more clusters)
    

    FCi = FC; % FC 
    for g=1:G % loop over range of gamma

        gamma = gam(g); % resolution parameter gamma

        
        qmax = -10; % modularity
        ciall = zeros(N,R); % community affiliation vectors
        for r=1:R % loop over runs (of Louvain)
            
            % compute community affiliation and modularity
            [ci, q] = community_louvain(FCi,gamma,[],'negative_asym');
            
            % save community affiliation CI vector
            % that maximizes modularity Q
            if(q>qmax)
                qmax = q;
                CI(:,g) = ci; % CI of max Q
                Q(g) = q; % max Q
            end
            
            ciall(:,r) = ci; % CIs of all runs

        end
        
        % consensus community affiliation over all Louvain runs
        CI2(:,g) = consensus_und(agreement(ciall),tau,reps); 
        

    end
    

        

end



function setup_plots()

    %% setup paths
    addpath(strcat(pwd,'/src/'))
    addpath(strcat(pwd,'/data/'))
    addpath(genpath(strcat(pwd,'/src/external/')))

end



function [rho,p,p_adj,num_comparisons] = get_relation(X,y,C, option)

        
    % X shape = Variables x Subjects
    num_comparisons = size(X,2);

    for c = 1:num_comparisons
        Xc = X(:,c);
        [rho(c), p(c)] = partialcorr(Xc,y,C,'type','spearman');
        
    end

    
    if option == "bonferroni"

        p_adj = p*num_comparisons;


    elseif option == "fdr"
        q = 0.05;
        
        [~, ~, ~, p_adj]=fdr_bh(p,q);

    else
    
        p_adj = p;

    end

end

   