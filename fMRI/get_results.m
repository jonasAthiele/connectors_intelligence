
%% Script for computing associations between fMRI centrality measures and intelligence

%% load data and set parameters
load ws_fc.mat
load node_network_assignment_200_17.mat
load table_descriptives.mat
acc = readtable('table_accuracy.csv');
acc_names = string(acc.Subject);
names_regions = string(readtable('Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv').ROIName);


NAMES_RUNS = ["rest","raven_run1","raven_run2", "raven_run3"];


S = length(names);
no_subs_all = S; % number of subjetcs

% number nodes
fields = fieldnames(fc);
firstFieldName = fields{1};
[N,~,~] = size(fc.(firstFieldName).rest_rr);


% set up parc plotter (https://github.com/faskowit/parc_plotter)
addpath(strcat(pwd,'/parcplot/'))
run("setup_data.m")

rescale_meth = 'rescale'; % rescale approach: rescale FC to be in the range 0 to 1
thresh_prop = 0.5; % threshold for FCs

%% get info of all subjects
accs_all = nan(S,1); % initializing
age_all = nan(S,1); % initializing
sex_all = nan(S,1); % initializing
for s=1:S

    name_split = split(names(s),"-");
    name_short = name_split(2);

    accs_all(s) = acc.Accuracy(strcmp(acc_names, names(s)));
    age_all(s) = table_descriptives.Age(strcmp(table_descriptives.Subject, name_short));
    sex_temp = table_descriptives.Sex(strcmp(table_descriptives.Subject, name_short));
    
    if sex_temp == "M"
        sex_all(s) = 0;
    elseif sex_temp == "F"
        sex_all(s) = 1;
    end

end


%% exclude high motion subjects
no_high_motion = true(S,1);
thresh_mean_fd = 0.5;
thresh_spike = 5; 
fd_mean = nan(S,1); % initializing
for s=1:S
    
    name_split = split(names(s),"-");
    name_short = name_split(2);

    fd_mean_temp = ones(length(NAMES_RUNS),1)*10;
    for r = 1:length(NAMES_RUNS)
    
        fd_r = head_motion.(name_short).(NAMES_RUNS(r));
        fd_r_mean = mean(fd_r);
        fd_r_max = max(fd_r);

        fd_mean_temp(r) = fd_r_mean;

        if fd_r_mean > thresh_mean_fd
            no_high_motion(s) = 0;
        end

        if fd_r_max > thresh_spike
            no_high_motion(s) = 0;
        end

    end
    
    fd_mean(s) = mean(fd_mean_temp);


end


names = names(no_high_motion); % exclude high motion subjects
fd_mean = fd_mean(no_high_motion); % exclude high motion subjects
S = length(names);

% compute mean FCs, participation coefficient and degree
accs = nan(S,1); % initializing 
age = nan(S,1); % initializing
sex = nan(S,1); % initializing
degree_rest = nan(N,S); % initializing
degree_raven = nan(N,S); % initializing
partiC_rest = nan(N,S); % initializing
partiC_raven = nan(N,S); % initializing
for s=1:S
    
    disp(s)
    name_split = split(names(s),"-");
    name_short = name_split(2);

    accs(s) = acc.Accuracy(strcmp(acc_names, names(s)));
    age(s) = table_descriptives.Age(strcmp(table_descriptives.Subject, name_short));
    sex_temp = table_descriptives.Sex(strcmp(table_descriptives.Subject, name_short));
    
    if sex_temp == "M"
        sex(s) = 0;
    elseif sex_temp == "F"
        sex(s) = 1;
    end

    numperms = size(fc.(name_short).rest_rr,3);
  

    degree_rest_p = nan(N,numperms); % initializing
    degree_raven_p = nan(N,numperms); % initializing
    partiC_rest_p = nan(N,numperms); % initializing
    partiC_raven_p = nan(N,numperms); % initializing
    for p = 1:numperms

        fc_temp = fc.(name_short).rest_rr;
        fc_temp = squeeze(fc_temp(:,:,p));
        fc_temp = threshold_proportional(fc_temp,thresh_prop);
        fc_rest_temp = rescale_fc(fc_temp,rescale_meth);
    
        fc_temp = fc.(name_short).raven_rr;
        fc_temp = squeeze(fc_temp(:,:,p));
        fc_temp = threshold_proportional(fc_temp,thresh_prop);
        fc_raven_temp = rescale_fc(fc_temp,rescale_meth);


        degree_rest_p(:,p) = mean(fc_rest_temp);
        degree_raven_p(:,p) = mean(fc_raven_temp);
    

        C = list7; % Yeo 7 networks as partitions
        W = fc_rest_temp;
        partiC_rest_p(:,p) = participation_coef(W,C);

        W = fc_raven_temp;
        partiC_raven_p(:,p) = participation_coef(W,C);

    end

    degree_rest(:,s) = mean(degree_rest_p,2);
    degree_raven(:,s) = mean(degree_raven_p,2);

    partiC_rest(:,s) = mean(partiC_rest_p,2);
    partiC_raven(:,s) = mean(partiC_raven_p,2);

end



%% get info of all subjects
% total sample
no_subjects = no_subs_all;
age_min = min(age_all);
age_max = max(age_all);
age_mean = mean(age_all);
age_sd = std(age_all);
no_males = sum(sex_all==0);
no_females = sum(sex_all==1);

rpm = accs_all*30;
rpm_min = min(rpm);
rpm_max = max(rpm);
rpm_mean = mean(rpm);
rpm_sd = std(rpm);

rpm36 = accs_all*36;
rpm36_min = min(rpm36);
rpm36_max = max(rpm36);
rpm36_mean = mean(rpm36);
rpm36_sd = std(rpm36);

handednessTable = readtable('Handedness.xls');
handednessTable = rmmissing(handednessTable);
% get handedness data
cleanSubjectIDs = erase(names, 'sub-');  % Now just ["M87100268", "M86000432", ...]
isMatch = ismember(handednessTable.URSI, cleanSubjectIDs);
handedness_table_valid = handednessTable(isMatch, :);
handedness = handedness_table_valid.WhatIsYourDominantHand_;
num_right = sum(handedness == 1);
num_left = sum(handedness == 0);

varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd",...
    "rpm36_min","rpm36_max","rpm36_mean","rpm36_sd","no_righthanded",...
    "no_lefthanded"];

data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd,...
    rpm36_min,rpm36_max,rpm36_mean,rpm36_sd, num_right, num_left];

info_all_subs = array2table(data,'VariableNames',varNames);
writetable(info_all_subs, "info_subjects_all.csv")


%% only subjects with long enoguh raven runs included for analysis
% get indices of subjects with at least 158 time points (normal length of rest scan) in rest and raven
tp_rr = nan(S,1); % initializing 
for s=1:S
    name_split = split(names(s),"-");
    name_short = name_split(2);
    tp_rr(s) =  tp.(name_short).rr; % time points used for analysis (rpm data was shortened to match rest data length = 158 frames)
end
idx_quest = tp_rr < 158;
idx_158 = ~idx_quest;

no_subjects = sum(idx_158==1);
age = age(idx_158);
sex = sex(idx_158);
accs = accs(idx_158);
fd_mean = fd_mean(idx_158);
degree_raven = degree_raven(:,idx_158);
degree_rest = degree_rest(:,idx_158);
partiC_raven = partiC_raven(:,idx_158); 
partiC_rest = partiC_rest(:,idx_158);
names = names(idx_158)';

age_min = min(age);
age_max = max(age);
age_mean = mean(age);
age_sd = std(age);

no_males = sum(sex==0);
no_females = sum(sex==1);

rpm = accs*30;
rpm_min = min(rpm);
rpm_max = max(rpm);
rpm_mean = mean(rpm);
rpm_sd = std(rpm);

rpm36 = accs*36;
rpm36_min = min(rpm36);
rpm36_max = max(rpm36);
rpm36_mean = mean(rpm36);
rpm36_sd = std(rpm36);

% get handedness info
handednessTable = readtable('table_handedness.csv');
handednessTable = rmmissing(handednessTable);
cleanSubjectIDs = erase(names, 'sub-');  % Now just ["M87100268", "M86000432", ...]
isMatch = ismember(handednessTable.URSI, cleanSubjectIDs);
handedness_table_valid = handednessTable(isMatch, :);
handedness = handedness_table_valid.DominantHand;
num_right = sum(handedness == 1);
num_left = sum(handedness == 0);

[r_acc_age, p_acc_age] = corr(accs,age,'type','Spearman'); % correlation between age and accuracies
[r_acc_sex, p_acc_sex] = corr(accs,sex,'type','Spearman');
[r_acc_fd, p_acc_fd]  = corr(accs,fd_mean,'type','Spearman');

% some information to store: age, intelligence scores, sex of participants
varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd",...
    "rpm36_min","rpm36_max","rpm36_mean","rpm36_sd","r_acc_age",...
    "p_acc_age","r_acc_sex","p_acc_sex","r_acc_fd","p_acc_fd",...
    "no_righthanded", "no_lefthanded"];

data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd,...
    rpm36_min,rpm36_max,rpm36_mean,rpm36_sd, r_acc_age, ...
    p_acc_age,r_acc_sex,p_acc_sex,r_acc_fd, p_acc_fd, num_right, num_left];

info_subs = array2table(data,'VariableNames',varNames);
writetable(info_subs, "info_subjects.csv")

%% Plot distribution of RPM scores
figure()
histogram(accs*36,10, 'FaceColor',[0.4 0.4 0.4]) % multiplied by 36 which is the number of rpm items in the EEG sample
box off
set(gca,'fontsize', 16) 

%% Get relations of centrality measures with intelligence
C = [age,sex,fd_mean]; % confounds
y = accs; % rpm accuracy

% degree
X = degree_raven - degree_rest;
[rho_deg,p_deg,p_adj_deg,~] = get_relation(X',y,C);
crange = [-0.5,0.5];
cMap = flipud(slanCM('RdBu'));
plot_nodes(rho_deg, crange, cMap, p_adj_deg)

% regions with associations with p < 0.05 (not corrected for multiple
% comparisons)
names_sig = names_regions(p_deg<0.05);
p_unadj = p_deg(p_deg<0.05);
p_adj = p_adj_deg(p_deg<0.05);
rho_sig = rho_deg(p_deg<0.05);

table_res_degree = array2table([names_sig, rho_sig, p_unadj, p_adj],...
    'VariableNames',["region","rho","p_uncorrected","p_corrected"]);
writetable(table_res_degree, "results_degree_uncorrected.csv")

% participation coefficient
X = partiC_raven - partiC_rest;
crange = [-0.5,0.5];
[rho_pc,p_pc,p_adj_pc,comparisons] = get_relation(X',y,C);
cMap = flipud(slanCM('RdBu'));
plot_nodes(rho_pc, crange, cMap, p_adj_pc)

% regions with associations with p < 0.05 (fdr corrected)
names_sig = names_regions(p_adj_pc<0.05);
p_unadj = p_pc(p_adj_pc<0.05);
p_adj = p_adj_pc(p_adj_pc<0.05);
rho_sig = rho_pc(p_adj_pc<0.05);
table_res_parti = array2table([names_sig, rho_sig, p_unadj, p_adj],...
    'VariableNames',["region","rho","p_uncorrected","p_corrected"]);
writetable(table_res_parti, "results_participation_corrected.csv")

%% Group average degree and pc for rest and rpm
% plot degree raven-rest
cMap = flipud(slanCM('RdBu'));
crange = [-0.05,0.05];
plot_nodes(mean(degree_raven-degree_rest,2), crange, cMap)

% plot degree raven
cMap = slanCM('Reds');
crange = [0,0.2];
plot_nodes(mean(degree_raven,2), crange, cMap)

% plot degree rest
cMap = slanCM('Reds');
crange = [0,0.2];
plot_nodes(mean(degree_rest,2), crange, cMap)


% plot pc raven-rest
cMap = flipud(slanCM('RdBu'));
crange = [-0.1,0.1];
plot_nodes(mean(partiC_raven-partiC_rest,2), crange, cMap)

% pc raven
cMap = slanCM('Reds');
crange = [0.6,0.9];
plot_nodes(mean(partiC_raven,2), crange, cMap)

% pc rest
cMap = slanCM('Reds');
crange = [0.6,0.9];
plot_nodes(mean(partiC_rest,2), crange, cMap)


%% correlation between pc and intelligence for P-FIT nodes vs rest of nodes
% load data
T = readtable('Schaefer200_Pfit_CoreRegions_Annotated.csv');  % Annotated table
corrs = rho_pc;  % Replace with your actual Nx1 Spearman correlations

% check matching size
assert(height(T) == length(corrs), 'Mismatch between table and correlation vector');

% define groups
isPfit = logical(T.IncludedInPfit);
pfit_corrs = corrs(isPfit);
nonpfit_corrs = corrs(~isPfit);

% observed mean difference
obs_diff = mean(pfit_corrs) - mean(nonpfit_corrs);

% permutation test
nPerm = 10000;
perm_diffs = zeros(nPerm,1);
n = length(corrs);

for i = 1:nPerm
    idx = randperm(n);
    shuffled_isPfit = isPfit(idx);
    
    perm_diffs(i) = mean(corrs(shuffled_isPfit)) - mean(corrs(~shuffled_isPfit));
end

% one-tailed p-value (is observed > chance?)
p = mean(perm_diffs >= obs_diff);

% display results
fprintf('Observed mean difference (P-FIT - non-P-FIT): %.4f\n', obs_diff);
fprintf('Permutation p-value (one-tailed): %.4f\n', p);

figure();
hold on;

% histogram
histogram(perm_diffs, 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5);

% vertical line for observed difference
xline(obs_diff, 'k', 'LineWidth', 2);

xlabel('Mean Difference (P-FIT - non-P-FIT)');
ylabel('Frequency');
%title('Permutation Test Null Distribution');

legend({'Null distribution', 'Observed difference'}, 'Location', 'northwest');

set(gca, 'FontSize', 12);
box on;
hold off;
exportgraphics(gcf, 'PermutationTestHistogram.png', 'Resolution', 300);

%% Plots for Supplement

% group mean fcs for rest and rpm
% reshape fcs for plotting
N_triu = N * (N - 1) / 2; % elements of upper triangle NxN matrix - N = number of nodes
fc_rest_s = nan(N_triu,no_subjects); % initializing
fc_raven_s = nan(N_triu,no_subjects); % initializing
for s=1:no_subjects
    
    name_split = split(names(s),"-");
    name_short = name_split(2);
    fc_temp = squeeze(mean(fc.(name_short).rest_rr,3));
    fc_rest_s(:,s) = fc_temp(triu(true(N),1));
    
    fc_temp = squeeze(mean(fc.(name_short).raven_rr,3));
    fc_raven_s(:,s) = fc_temp(triu(true(N),1));
end

% plot fc raven-rest
figure()
cMap = flipud(slanCM('RdBu'));
plot_matrix(mean(fc_raven_s - fc_rest_s,2), ones(19900,1), 200, [-0.3,0.3], cMap, "nodes", [names7(1:7), names7(1:7)], "whole")

% plot fc raven
figure()
cMap = flipud(slanCM('RdBu'));
plot_matrix(mean(fc_raven_s,2), ones(19900,1), 200, [-0.8,0.8], cMap, "nodes", [names7(1:7), names7(1:7)], "whole")

% plot fc rest
figure()
cMap = flipud(slanCM('RdBu'));
plot_matrix(mean(fc_rest_s,2), ones(19900,1), 200, [-0.8,0.8], cMap, "nodes", [names7(1:7), names7(1:7)], "whole")


%% spatial correlations
% calculate mean across subjects for each region
mean_partiC_raven = mean(partiC_raven, 2);
mean_degree_raven = mean(degree_raven, 2);
mean_partiC_rest = mean(partiC_rest, 2);
mean_degree_rest = mean(degree_rest, 2);

% correlation between Rest and Raven for Participation Coefficient
[r_pc_rest_raven, p_pc_rest_raven] = corr(mean_partiC_rest, mean_partiC_raven);

% Correlation between Rest and Raven for Degree
[r_deg_rest_raven, p_deg_rest_raven] = corr(mean_degree_rest, mean_degree_raven);

% differences between Raven and Rest
diff_pc = mean_partiC_raven - mean_partiC_rest;
diff_deg = mean_degree_raven - mean_degree_rest;

% correlation between differences for Participation Coefficient and Degree
[r_diff_pc_deg, p_diff_pc_deg] = corr(diff_pc, diff_deg);

% correlation of correlations rpm and Participation Coefficient vs rpm and Degree
[r_spatial, p_spatial] = corr(rho_pc, rho_deg);


% display results
fprintf('Correlation PC Rest vs Raven: r=%.3f, p=%.3g\n', r_pc_rest_raven, p_pc_rest_raven);
fprintf('Correlation Degree Rest vs Raven: r=%.3f, p=%.3g\n', r_deg_rest_raven, p_deg_rest_raven);
fprintf('Correlation differences PC vs Degree: r=%.3f, p=%.3g\n', r_diff_pc_deg, p_diff_pc_deg);
fprintf('Spatial correlation between degree and participation coefficient for RPM associations: r=%.3f, p=%.3g\n', r_spatial, p_spatial);

%% Functions
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

function plot_nodes(vals,crange, cMap, p)
    

    % prepare
    annotName = 'schaefer200-yeo17';
    surf = load([pwd '/parcplot/data/fsaverage/mat/' 'fsaverage_inflated.mat']) ;
    surfStruct = surf.surfStruct ;
    annots = load([pwd '/parcplot/data/fsaverage/mat/' 'fsaverage_annots.mat']) ;
    annotMap = annots.allAnnots ;

    % plot
    dataVec = vals;
    parc_plot(surfStruct,annotMap,annotName,dataVec,'cMap',cMap,'valRange',crange)
  
    figure
    colormap(cMap); 
    caxis(crange)
    hCB = colorbar('north');
    set(gca,'Visible',false)
    hCB.Position = [0.1 0.4 0.8 0.12];
    hCB.FontSize = 20;
    hCB.Ticks = crange(1):crange(2)/2:crange(2);

    if nargin == 4
        dataVec = zeros(200,1);
        dataVec(p<0.05) = 0.6;
    
        parc_plot(surfStruct,annotMap,annotName,dataVec,...
            'cMap',cMap,'valRange',[-1,1])
    end

end

function [rho,p_values,corrected_p,num_comparisons] = get_relation(X,y,C)

        
    % X shape = Variables x Subjects
    num_comparisons = size(X,2);
    rho = nan(num_comparisons,1);
    p_values = nan(num_comparisons,1);
    for c = 1:num_comparisons
        Xc = X(:,c);
        [rho(c), p_values(c)] = partialcorr(Xc,y,C,'type','spearman');
        
    end

    corrected_p = mafdr(p_values, 'BHFDR', true);

end


function plot_matrix(vals, p, N, range, cMap, label_axes, label_ticks, option)
    
    matrix_form = false;
    if size(vals,1) == size(vals,2) %check if vals are in matrix form
        
        matrix_form = true;

    end

    plot_nodewise = false;
    vals_m = nan(N);
    p_m = nan(N);
    if option == "whole"
        
        if matrix_form
            vals_m = vals;
            p_m = p;
        else
            vals_m(triu(true(N),1)) = vals;
            p_m(triu(true(N),1)) = p;
        end
        plot_nodewise = true;
        fontSize = 6;
        line_x = [12 28 39 50 56 74 100 112 130 141 156 164 183] + 1.5;
        ticks_x = [6, 20, 33.5, 44.5, 53, 65, 87, 106, 121, 135.5, 148.5, 160, 173.5, 191.5] + 1;
    
    elseif option == "condition"
        if matrix_form
            vals_m = vals;
            p_m = p;
        else
            vals_m(triu(true(N),1)) = vals;
            p_m(triu(true(N),1)) = p;
        end
        fontSize = 18;
    else
        if matrix_form
            vals_m = vals;
            p_m = p;
        else
            vals_m(triu(true(N))) = vals;
            p_m(triu(true(N))) = p;
        end
        fontSize = 18;
    end
    

    
    h=imagesc(vals_m);
    if plot_nodewise
        line(repmat(line_x,2,1),repmat([0;N+1], 1, length(line_x)), 'Color', 'black'); % vertical
        line(repmat([0;N+1], 1, length(line_x)), repmat(line_x,2,1), 'Color', 'black'); % horizontal
    else
        line(repmat([1.5:1:N-0.5],2,1),repmat([0;N+1], 1, N-1), 'Color', 'black'); % vertical
        line(repmat([0;N+1], 1, N-1), repmat([1.5:1:N-0.5],2,1), 'Color', 'black'); % horizontal
    end
    
    ax = gca;
    ax.TickLength = [0 0];
    caxis([-max(abs(range)),max(abs(range))]); % centering white color
    set(h, 'AlphaData', ~isnan(vals_m));
%     cMap = redblue;
    colormap(cMap)  
    bb = colorbar;
    set(bb, 'ylim', range)
   
    [indp1, indp2] = find(p_m<0.05);
    for c = 1: 1:length(indp1)
       text(indp2(c),indp1(c), '*', 'HorizontalAlignment', 'center','fontsize',fontSize)
    end


    if plot_nodewise
        
        xticks(ticks_x);
        yticks(ticks_x);
    
    elseif N == 7
        ticks_x = 1:7;
        xticks(ticks_x);
        yticks(ticks_x);

    elseif N == 17
        ticks_x = 1:17;
        xticks(ticks_x);
        yticks(ticks_x);
    end



    yticklabels(label_ticks)
    xticklabels(label_ticks)

    

    xlabel(label_axes) 
    ylabel(label_axes) 


end



