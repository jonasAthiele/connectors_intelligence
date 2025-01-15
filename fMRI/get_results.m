
%% Load data and set parameters
load ws_fc.mat
load node_network_assignment_200_17.mat
load table_descriptives.mat
acc = readtable('table_accuracy.csv');
acc_names = string(acc.Subject);
names_regions = string(readtable('Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv').ROIName);


NAMES_RUNS = ["rest","raven_run1","raven_run2", "raven_run3"];


S = length(names);
no_subs_all = S; % number of subjetcs

setup_plots(); % some function to setup parc_plotter

rescale_meth = 'rescale'; % rescale FC to be in the range 0 to 1
thresh_prop = 0.5; % threshold for FCs

%% Get info of all subjects
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


%% Exclude high motion subjects
no_high_motion = true(S,1);
thresh_mean_fd = 0.5;
thresh_spike = 5; 
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

% Compute mean FCs, participation coefficient and degree
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


        fc_rest_p(:,p) = fc_rest_temp(triu(true(200),1));
        fc_raven_p(:,p) = fc_raven_temp(triu(true(200),1));

    end

    fc_rest(:,s) = squeeze(mean(fc_rest_p,2));
    fc_raven(:,s) = mean(fc_raven_p,2);

    degree_rest(:,s) = mean(degree_rest_p,2);
    degree_raven(:,s) = mean(degree_raven_p,2);

    partiC_rest(:,s) = mean(partiC_rest_p,2);
    partiC_raven(:,s) = mean(partiC_raven_p,2);

end

% get indices of subjects with at least 158 time points (normal length of rest scan) in rest and raven
for s=1:S
    disp(s)
    name_split = split(names(s),"-");
    name_short = name_split(2);
    tp_rr(s) =  tp.(name_short).rr;
    

end
idx_quest = tp_rr < 158;
idx_158 = ~idx_quest;


%% Get some info of subjects
% some information to store: age, intelligence scores, sex of participants
varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd",...
    "rpm36_min","rpm36_max","rpm36_mean","rpm36_sd","r_acc_age","p_acc_age","r_acc_sex","p_acc_sex","r_acc_fd","p_acc_fd"];


% included subjects
no_subjects = sum(idx_158==1);
age_min = min(age(idx_158));
age_max = max(age(idx_158));
age_mean = mean(age(idx_158));
age_sd = std(age(idx_158));
no_males = sum(sex(idx_158)==0);
no_females = sum(sex(idx_158)==1);

rpm = accs(idx_158)*30;
rpm_min = min(rpm);
rpm_max = max(rpm);
rpm_mean = mean(rpm);
rpm_sd = std(rpm);

rpm36 = accs(idx_158)*36;
rpm36_min = min(rpm36);
rpm36_max = max(rpm36);
rpm36_mean = mean(rpm36);
rpm36_sd = std(rpm36);

[r_acc_age, p_acc_age] = corr(accs(idx_158)',age(idx_158)','type','Spearman'); % correlation between age and accuracies
[r_acc_sex, p_acc_sex] = corr(accs(idx_158)',sex(idx_158)','type','Spearman');
[r_acc_fd, p_acc_fd]  = corr(accs(idx_158)',fd_mean(idx_158)','type','Spearman');

data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd,...
    rpm36_min,rpm36_max,rpm36_mean,rpm36_sd, r_acc_age, ...
    p_acc_age,r_acc_sex,p_acc_sex,r_acc_fd, p_acc_fd];

info_subs = array2table(data,'VariableNames',varNames);
writetable(info_subs, "info_subjects.csv")

% total sample
varNames = ["no_subjects","age_min","age_max","age_mean","age_sd",...
    "no_males","no_females","rpm_min","rpm_max","rpm_mean","rpm_sd",...
    "rpm36_min","rpm36_max","rpm36_mean","rpm36_sd"];

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

data = [no_subjects,age_min,age_max,age_mean,age_sd,...
    no_males,no_females,rpm_min,rpm_max,rpm_mean,rpm_sd,...
    rpm36_min,rpm36_max,rpm36_mean,rpm36_sd];

info_all_subs = array2table(data,'VariableNames',varNames);
writetable(info_all_subs, "info_subjects_all.csv")

%% Plot distribution of RPM scores
figure()
histogram(accs(idx_158)*36,10, 'FaceColor',[0.4 0.4 0.4])
box off
set(gca,'fontsize', 16) 

%% Get relations of centrality measures with intelligence
C = [age',sex',fd_mean']; % confounds
y = accs'; % rpm accuracy

% degree
X = degree_raven(:,idx_158) - degree_rest(:,idx_158);
[rho,p,p_adj,comparisons] = get_relation(X',y(idx_158),C(idx_158,:), "fdr");
crange = [-0.5,0.5];
cMap = flipud(slanCM('RdBu'));
plot_nodes(rho, p_adj, crange, cMap)
plot_nodes(rho, p, crange, cMap)

% table degree p < 0.05 (not corrected)
names_sig = names_regions(p<0.05);
p_sig = p(p<0.05);
rho_sig = rho(p<0.05);

table_res_degree = array2table([names_sig, rho_sig', p_sig'],...
    'VariableNames',["region","rho","p_uncorrected"]);
writetable(table_res_degree, "results_degree_uncorrected.csv")


% participation coefficient
X = partiC_raven(:,idx_158) - partiC_rest(:,idx_158);
crange = [-0.5,0.5];
[rho,p,p_adj,comparisons] = get_relation(X',y(idx_158),C(idx_158,:), "fdr");
cMap = flipud(slanCM('RdBu'));
plot_nodes(rho, p_adj, crange, cMap)
plot_nodes(rho, p, crange, cMap)


% table participation coeff p < 0.05 (not corrected)
names_sig = names_regions(p<0.05);
p_sig = p(p<0.05);
rho_sig = rho(p<0.05);
table_res_parti = array2table([names_sig, rho_sig', p_sig'],...
    'VariableNames',["region","rho","p_uncorrected"]);
writetable(table_res_parti, "results_participation_uncorrected_50_100perms.csv")

% fdr corrected 
names_sig = names_regions(p_adj<0.05);
p_sig = p_adj(p_adj<0.05);
rho_sig = rho(p_adj<0.05);
table_res_parti = array2table([names_sig, rho_sig', p_sig'],...
    'VariableNames',["region","rho","p_uncorrected"]);
writetable(table_res_parti, "results_participation_corrected_50_100perms.csv")


%% Plots for Supplement

% reshape fcs for plotting
for s=1:S
    
    name_split = split(names(s),"-");
    name_short = name_split(2);
    fc_temp = squeeze(mean(fc.(name_short).rest_rr,3));
    fc_rest_s(:,s) = fc_temp(triu(true(200),1));
    
    fc_temp = squeeze(mean(fc.(name_short).raven_rr,3));
    fc_raven_s(:,s) = fc_temp(triu(true(200),1));
end

% plot fc raven-rest
cMap = flipud(slanCM('RdBu'));
plot_matrix(mean(fc_raven_s(:,idx_158) - fc_rest_s(:,idx_158),2), ones(19900,1), 200, [-0.3,0.3], cMap, "nodes", [names7(1:7), names7(1:7)], "whole")


% plot fc raven
cMap = flipud(slanCM('RdBu'));
plot_matrix(mean(fc_raven_s(:,idx_158),2), ones(19900,1), 200, [-0.8,0.8], cMap, "nodes", [names7(1:7), names7(1:7)], "whole")

% plot fc rest
cMap = flipud(slanCM('RdBu'));
plot_matrix(mean(fc_rest_s(:,idx_158),2), ones(19900,1), 200, [-0.8,0.8], cMap, "nodes", [names7(1:7), names7(1:7)], "whole")


% plot degree raven-rest
cMap = flipud(slanCM('RdBu'));
crange = [-0.05,0.05];
plot_nodes(mean(degree_raven(:,idx_158)-degree_rest(:,idx_158),2), ones(200,1), crange, cMap)

% plot degree raven
cMap = slanCM('Reds');
crange = [0,0.2];
plot_nodes(mean(degree_raven(:,idx_158),2), ones(200,1), crange, cMap)

% plot degree rest
cMap = slanCM('Reds');
crange = [0,0.2];
plot_nodes(mean(degree_rest(:,idx_158),2), ones(200,1), crange, cMap)


% plot pc raven-rest
cMap = flipud(slanCM('RdBu'));
crange = [-0.1,0.1];
plot_nodes(mean(partiC_raven(:,idx_158)-partiC_rest(:,idx_158),2), ones(200,1), crange, cMap)

% pc raven
cMap = slanCM('Reds');
crange = [0.6,0.9];
plot_nodes(mean(partiC_raven(:,idx_158),2), ones(200,1), crange, cMap)

% pc rest
cMap = slanCM('Reds');
crange = [0.6,0.9];
plot_nodes(mean(partiC_rest(:,idx_158),2), ones(200,1), crange, cMap)

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

function setup_plots()

    %% setup paths
    addpath(strcat(pwd,'/parcplot/src/'))
    addpath(strcat(pwd,'/parcplot/data/'))
    addpath(genpath(strcat(pwd,'/parcplot/src/external/')))

end

function plot_nodes(vals, p, crange, cMap)
    

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

    dataVec = zeros(200,1);
    dataVec(p<0.05) = 0.6;

    parc_plot(surfStruct,annotMap,annotName,dataVec,...
        'cMap',cMap,'valRange',[-1,1])


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



