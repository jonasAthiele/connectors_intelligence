
% Script for computing relationship between mse and rpm scores

% behavioral data
path_data = pwd;
beh_file = fullfile(path_data,"descr_behav_data.csv");
BEHAVIORAL_DATA = readtable(beh_file);
SUBJECT_IDS = BEHAVIORAL_DATA.ID;

% folder with mse values (brain_measures)
FOLDER_DIRECTORY_MSE = fullfile(path_data,"data","eeg_mse");

NUMPERMS = 100; % number of permutations
NUMEPOCHS = 3; % numbers of epochs per trial (1 = first 10 seconds, 2 = 20 seconds, 3 = 30 seconds)
idx_epochs = 1:NUMEPOCHS;

%% get ids of participants
file_names={dir(FOLDER_DIRECTORY_MSE).name};
file_names=file_names(~ismember(file_names,{'.','..'}));

for s = 1:length(file_names)
    f = file_names(s);
    fsplit = split(f,'_');
    id(s) = str2double(string(fsplit(1)));
end
id = sort(id);

%% exclude subjects due to missing data
id(id == 72) = []; % no age
id(id == 140) = []; % no age


%% exclude subjects with short data length
n_subjects = length(id); % number subjects
trials_raven = nan(n_subjects,1); % initializing
epochs_rest = nan(n_subjects,1); % initializing
age_all = nan(n_subjects,1); % initializing
sex_all = nan(n_subjects,1); % initializing
rpm_all = nan(n_subjects,1); % initializing
for s = 1:n_subjects
 
    id_s = id(s);

    % descriptive data
    idx_subject = find(BEHAVIORAL_DATA.ID == id_s);

    age_all(s) = BEHAVIORAL_DATA.WIEK(idx_subject);
    sex_temp = BEHAVIORAL_DATA.PLEC(idx_subject);
    if sex_temp == "M" % male
        sex_translate = 0;
    elseif sex_temp == "K" % female
        sex_translate = 1;
    end
    sex_all(s) = sex_translate;
    rpm_all(s) = BEHAVIORAL_DATA.RAVEN(idx_subject);

    % check for trial/epoch number
    file_name = append(string(id_s), '_brain_measures.mat');
    file_path = fullfile(FOLDER_DIRECTORY_MSE, file_name);
    
    load(file_path)

    trials_raven(s) = size(brain_measures.raven.n_epochs,2);
    epochs_rest(s) = size(brain_measures.rest.n_epochs,2);

end

[~,~,n_electrodes,n_timescales] = size(brain_measures.raven.mse);

% criteria
n_trials_raven_min = 36/2; % subjects with at least half of the rpm trials
n_epochs_rest_min = 33; % each rest run consists of 360 seconds epoched into 10 sec epochs

id_too_short = (trials_raven < n_trials_raven_min | epochs_rest < n_epochs_rest_min);
id(id_too_short) = [];

%% get mse values for included subjects
n_subjects = length(id);
avg_mse_raven_perm = nan(NUMPERMS,n_subjects,n_electrodes,n_timescales); % initializing
avg_mse_rest_perm = nan(NUMPERMS,n_subjects,n_electrodes,n_timescales); % initializing
age = nan(n_subjects,1); % initializing
sex = nan(n_subjects,1); % initializing
trials = nan(n_subjects,1); % initializing
rpm = nan(n_subjects,1); % initializing

for s = 1:n_subjects
 
    id_s = id(s);
    file_name = append(string(id_s), '_brain_measures.mat');
    file_path = fullfile(FOLDER_DIRECTORY_MSE, file_name);
    load(file_path)

    n_trials_raven = size(brain_measures.raven.n_epochs,2);
    n_epochs_rest = size(brain_measures.rest.n_epochs,2);

    idx_subject = find(BEHAVIORAL_DATA.ID == id_s);

    age(s) = BEHAVIORAL_DATA.WIEK(idx_subject);
    sex_temp = BEHAVIORAL_DATA.PLEC(idx_subject);
    if sex_temp == "M" % male
        sex_translate = 0;
    elseif sex_temp == "K" % female
        sex_translate = 1;
    end
    sex(s) = sex_translate;
    rpm(s) = BEHAVIORAL_DATA.RAVEN(idx_subject);
    trials(s) = n_trials_raven;

    n_epoch_rest_select = NUMEPOCHS * n_trials_raven_min;
    n_trials_raven_select = n_trials_raven_min;
    
    % make sure rest and rpm data lengths are equal
    if n_epoch_rest_select > n_epochs_rest_min
        n_trials_raven_select = floor(n_epochs_rest_min/NUMEPOCHS);
        n_epoch_rest_select = n_trials_raven_select*NUMEPOCHS;
    end

    for p = 1:NUMPERMS
        perm_idx_trials_raven = randperm(n_trials_raven);
        perm_idx_epochs_rest = randperm(n_epochs_rest);

        idx_raven = perm_idx_trials_raven(1:n_trials_raven_select);
        idx_rest = perm_idx_epochs_rest(1:n_epoch_rest_select);
    
        mse_raven = brain_measures.raven.mse(idx_raven,idx_epochs,:,:);
        avg_mse_raven_perm(p,s,:,:) = squeeze(mean(mean(mse_raven, 1), 2));  % average across trials and epochs - result: [electrodes x timescales]
        
        mse_rest = squeeze(brain_measures.rest.mse);
        avg_mse_rest_perm(p,s,:,:) = squeeze(mean(mse_rest(idx_rest,:,:))); % average across trials and epochs - result: [electrodes x timescales]
    
    end
    
end  

avg_mse_raven = squeeze(mean(avg_mse_raven_perm,1)); % average across permutations
avg_mse_rest = squeeze(mean(avg_mse_rest_perm,1)); % average across permutations

mse = avg_mse_raven - avg_mse_rest; % difference rpm minus rest

%% Compute correlation and p for association between mse and
% rpm for each electrode and timescale
confounds = [sex, age, trials];
[rho_rs, p_rs] = compute_partial_spearman(mse, rpm, confounds);

% permutation test based on locations (to get corrected p-vals)
electrode_names = string({readlocs('channel_locations_64.ced').labels});
nPerms = 100;
alpha = 0.05;

chanlocs = readlocs('channel_locations_64.ced');
coords3D = [[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]']; 
distMat = squareform(pdist(coords3D));  % Euclidean distances in 3D

% figure out a neighboring radius
distMat_noDiag = distMat(~eye(size(distMat)));
meanDist = mean(distMat_noDiag);
disp(meanDist)
neighbor_radius = 0.5 * meanDist;  

adjacency = distMat < neighbor_radius;
adjacency(1:size(adjacency,1)+1:end) = 0;  % zero diagonal (no self-edges)
% Number of neighbors per electrode
numNeighbors = sum(adjacency, 2);

% plot histogram
histogram(numNeighbors)
xlabel('Number of neighbors per electrode')
ylabel('Count of electrodes')

% plot neighbors for electrode 1
elec = 1;
neighbors = find(adjacency(elec, :));

figure;
scatter3(coords3D(:,1), coords3D(:,2), coords3D(:,3), 'filled'); hold on
scatter3(coords3D(elec,1), coords3D(elec,2), coords3D(elec,3), 200, 'r', 'filled')
scatter3(coords3D(neighbors,1), coords3D(neighbors,2), coords3D(neighbors,3), 200, 'g', 'filled')
title('Neighbors of electrode 1 (red)')

% initialize permutation function data
permute_corr_func('set_data', mse, rpm, confounds);

% run cluster permutation test
[obsClusters, obsClusterStats, correctedPvals, maxClusterStatsPerm, corrected_p_matrix] = ...
    cluster_permutation_test_signed(rho_rs, p_rs, adjacency, @(perm) permute_corr_func(perm), nPerms, alpha);



%% figure: associations mse (whole rpm trials) - rpm sum scores 

alpha = 0.05;
timescale_ticks = 1:20;
n_electrodes = length(electrode_names);

f=figure;
% plot the correlation matrix
imagesc(rho_rs);

% customize plot
f.Position = [100 100 900 800];
ax = gca;
ax.TickLength = [0 0];
xticks(1:20);
yticks(1:n_electrodes);
yticklabels(electrode_names);
xticklabels(timescale_ticks);
set(gca, 'XTickLabelRotation', 0);  % rotate x tick labels to horizontal
set(gca, 'FontSize', 9);  % tick label size
colormap(flipud(slanCM('RdBu')));
caxis([-0.35 0.35])
colorbar;
cb = colorbar;
cb.FontSize = 9.5;
ylabel('Electrodes', 'FontSize', 12);
xlabel('Timescales', 'FontSize', 12);
hold on;
% grid for columns (timescales)
for x = 0.5 : 1 : length(timescale_ticks)+0.5
    line([x x], [0.5 n_electrodes+0.5], 'Color', [0 0 0], 'LineWidth', 0.5);
end
% grid for rows (electrodes)
for y = 0.5 : 1 : n_electrodes+0.5
    line([0.5 length(timescale_ticks)+0.5], [y y], 'Color', [0 0 0], 'LineWidth', 0.5);
end


% overlay asterisks
[row, col] = find(corrected_p_matrix < alpha);
yOffset = +0.3; 
for k = 1:length(row)
    text(col(k), row(k) + yOffset, '*', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'Color', 'k');
end

hold off;
box off;
f.Position = [10 10 450 800]; % [left bottom width height] in pixels
exportgraphics(f, 'mse_raven_30sec.png', 'Resolution', 300);


%% group mean rest, rpm, delta rpm-rest - average across electrodes
 
% average across electrodes per subject at each timescale
mean_mse_raven_subj = squeeze(mean(avg_mse_raven, 2)); % subjects x timescales
mean_mse_rest_subj  = squeeze(mean(avg_mse_rest, 2));  % subjects x timescales

nTimescales = size(avg_mse_raven, 3);
timescales = 1:nTimescales;

% paired t-tests per timescale across subjects
p_values = nan(1, nTimescales);
for t = 1:nTimescales
    [~, p] = ttest(mean_mse_rest_subj(:, t), mean_mse_raven_subj(:, t));
    p_values(t) = p;
end

% FDR correction for multiple comparisons
corrected_p = mafdr(p_values, 'BHFDR', true);

% calculate group averages (mean across subjects and electrodes)
group_avg_raven = squeeze(mean(mean(avg_mse_raven, 1), 2)); % 1 x timescales
group_avg_rest  = squeeze(mean(mean(avg_mse_rest, 1), 2));  % 1 x timescales

sig_timescales = find(corrected_p < 0.05);

fig = figure; hold on;
plot(timescales, group_avg_rest, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
plot(timescales, group_avg_raven, 'k', 'LineWidth', 1.5);

ylims = ylim;
set(gca, 'FontSize', 13);  % tick label size
set(gca, 'XTickLabelRotation', 0);  % rotate x tick labels to horizontal

xlabel('Timescale', 'FontSize',14);
ylabel('Sample entropy','FontSize',14);

% add vertical grid lines separating timescales using 'line' (behind plot)
for x = timescales
    line([x x], ylims, 'Color', [0.8 0.8 0.8], 'LineStyle', '-', 'LineWidth', 0.7, 'HandleVisibility','off');
end

ax = gca;
ax.XTick = timescales;
ax.XLim = [min(timescales)-0.5 max(timescales)+0.5];
grid on;

y_asterisk = ylims(2) - 0.05 * (ylims(2) - ylims(1)); % a bit below top
for t = sig_timescales
    plot(t, y_asterisk, '*k', 'MarkerSize', 6);  % original size
end

plot(timescales, group_avg_rest, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
plot(timescales, group_avg_raven, 'k', 'LineWidth', 1.5);

% legend adjustments
lgd = legend({'  Rest', '  Raven'}, 'Location', 'southeast');
lgd.FontSize = 14;
lgd.Box = 'off';

hold off;

fig.Position = [100 100 600 400]; % [left bottom width height] in pixels
exportgraphics(fig, 'mse_average_plot.png', 'Resolution', 300);



%% group mean rest, raven, delta raven-rest - not average across electrodes

% calculate group means (average over subjects)
mean_rest = squeeze(mean(avg_mse_rest, 1));   % [electrodes x timescales]
mean_raven = squeeze(mean(avg_mse_raven, 1)); % [electrodes x timescales]

% number of electrodes and timescales
n_electrodes = size(mean_rest, 1);

% transform rest and rpm MSE with sqrt
mean_rest_t = sqrt(mean_rest);   % [electrodes x timescales]
mean_raven_t = sqrt(mean_raven); % same
% keep the original diff (no transform)
diff_mse = mean_raven - mean_rest;

% define timescales and color limits
timescale_ticks = 1:20;
cax_rest = [0.7 1.4];   % adjust based on transformed data
cax_raven = [0.7 1.4];  % same here
cax_diff = [-0.12 0.12]; % based on original differences

fig = figure;

% plot sqrt-transformed rest MSE
ax1 = subplot(1,3,1);
imagesc(mean_rest_t(:, timescale_ticks));
%title('Rest MSE (sqrt)');
yticks(1:n_electrodes);
yticklabels(electrode_names);
xticks(timescale_ticks);
xticklabels(string(timescale_ticks));

set(gca, 'FontSize', 12);  % tick label size
set(gca, 'XTickLabelRotation', 0);  % rotate x tick labels to horizontal
set(gca, 'TickLength', [0.005 0.005]);
caxis(cax_rest);
colormap(ax1, slanCM('Reds'));
colorbar;
cb = colorbar;
cb.FontSize = 13;
ylabel('Electrodes', 'FontSize', 16);
xlabel('Timescales', 'FontSize', 16);
hold on;
% grid for columns (timescales)
for x = 0.5 : 1 : length(timescale_ticks)+0.5
    line([x x], [0.5 n_electrodes+0.5], 'Color', [0 0 0], 'LineWidth', 0.5);
end
% grid for rows (electrodes)
for y = 0.5 : 1 : n_electrodes+0.5
    line([0.5 length(timescale_ticks)+0.5], [y y], 'Color', [0 0 0], 'LineWidth', 0.5);
end
hold off;
box off;

% plot sqrt-transformed Raven MSE
ax2 = subplot(1,3,2);
imagesc(mean_raven_t(:, timescale_ticks));
%title('RPM MSE (sqrt)');
yticks(1:n_electrodes);
yticklabels(electrode_names);
xticks(timescale_ticks);
xticklabels(string(timescale_ticks));

set(gca, 'FontSize', 12);  % tick label size
set(gca, 'XTickLabelRotation', 0);  % rotate x tick labels to horizontal
set(gca, 'TickLength', [0.005 0.005]);
caxis(cax_raven);
colormap(ax2, slanCM('Reds'));
colorbar;
cb = colorbar;
cb.FontSize = 13;
ylabel('Electrodes', 'FontSize', 16);
xlabel('Timescales', 'FontSize', 16);
hold on;
% grid for columns (timescales)
for x = 0.5 : 1 : length(timescale_ticks)+0.5
    line([x x], [0.5 n_electrodes+0.5], 'Color', [0 0 0], 'LineWidth', 0.5);
end
% grid for rows (electrodes)
for y = 0.5 : 1 : n_electrodes+0.5
    line([0.5 length(timescale_ticks)+0.5], [y y], 'Color', [0 0 0], 'LineWidth', 0.5);
end
hold off;
box off;

% plot raw difference (RPM - Rest)
ax3 = subplot(1,3,3);
imagesc(diff_mse(:, timescale_ticks));
%title('Difference (RPM - Rest)');
yticks(1:n_electrodes);
yticklabels(electrode_names);
xticks(timescale_ticks);
xticklabels(string(timescale_ticks));
set(gca, 'FontSize', 12);  % tick label size
set(gca, 'XTickLabelRotation', 0);  % rotate x tick labels to horizontal
set(gca, 'TickLength', [0.005 0.005]);
caxis(cax_diff);
colormap(ax3, flipud(slanCM('RdBu')));
colorbar;
cb = colorbar;
cb.FontSize = 13;
ylabel('Electrodes', 'FontSize', 16);
xlabel('Timescales', 'FontSize', 16);
hold on
% Grid for columns (timescales)
for x = 0.5 : 1 : length(timescale_ticks)+0.5
    line([x x], [0.5 n_electrodes+0.5], 'Color', [0 0 0], 'LineWidth', 0.5);
end
% Grid for rows (electrodes)
for y = 0.5 : 1 : n_electrodes+0.5
    line([0.5 length(timescale_ticks)+0.5], [y y], 'Color', [0 0 0], 'LineWidth', 0.5);
end
hold off;
box off;

set(fig, 'Color', 'w');
fig.Position = [100 100 2200 1100]; % [left bottom width height] in pixels
exportgraphics(fig, 'mse_plot.png', 'Resolution', 300);

%% Supplementary Figures

%% distribution of rpm scores
figure()
histogram(rpm,10, 'FaceColor',[0.4 0.4 0.4])
box off
set(gca,'fontsize', 16) 



%% Some information to store

% total sample:
no_subjects = length(rpm_all);
age_min = min(age_all);
age_max = max(age_all);
age_mean = mean(age_all);
age_sd = std(age_all);
no_males = sum(sex_all==0);
no_females = sum(sex_all==1);

rpm_min = min(rpm_all);
rpm_max = max(rpm_all);
rpm_mean = mean(rpm_all);
rpm_sd = std(rpm_all);


data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd];

varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd"];

info_all_subs = array2table(data,'VariableNames',varNames);
writetable(info_all_subs, "info_subjects_all.csv")


% included subjects

no_subjects = length(rpm);
age_min = min(age);
age_max = max(age);
age_mean = mean(age);
age_sd = std(age);
no_males = sum(sex==0);
no_females = sum(sex==1);

rpm_min = min(rpm);
rpm_max = max(rpm);
rpm_mean = mean(rpm);
rpm_sd = std(rpm);

[r_rpm_age, p_rpm_age] = corr(rpm,age,'type','Spearman'); % correlation between age and RPM scores
[r_rpm_sex, p_rpm_sex] = corr(rpm,sex,'type','Spearman');
[r_rpm_epochs, p_rpm_epochs] = corr(rpm,trials,'type','Spearman');


data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd, ...
    r_rpm_age,p_rpm_age,r_rpm_sex,p_rpm_sex,-r_rpm_epochs,p_rpm_epochs]; 
% minus at correlation with number of epochs to get correlation between rpm and removed epochs

varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd",...
    "r_rpm_age", "p_rpm_age","r_rpm_sex","p_rpm_sex","r_rpm_epochs",...
    "p_rpm_epochs"];

info_subs = array2table(data,'VariableNames',varNames);

writetable(info_subs, "info_subjects.csv")