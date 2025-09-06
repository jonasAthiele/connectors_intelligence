function [obsClusters, obsClusterStats, correctedPvals, maxClusterStatsPerm, corrected_p_matrix] = ...
    cluster_permutation_test_signed(rho_rs, p_rs, adjacency, compute_corr_func, nPerms, alpha)
% CLUSTER_PERMUTATION_TEST_SIGNED
% Cluster-based permutation test for both positive and negative correlations, with single output matrix.
%
% INPUTS:
%   rho_rs          - observed correlation matrix (channels x conditions)
%   p_rs            - observed p-value matrix (channels x conditions)
%   adjacency       - adjacency matrix (channels x channels), binary
%   compute_corr_func - function handle returning [rho_perm, p_perm] for a given permutation
%   nPerms          - number of permutations
%   alpha           - cluster-forming p-value threshold
%
% OUTPUTS:
%   obsClusters         - struct with .pos and .neg fields: cell arrays of clusters
%   obsClusterStats     - struct with .pos and .neg fields: cluster stats
%   correctedPvals      - struct with .pos and .neg fields: corrected p-values
%   maxClusterStatsPerm - struct with .pos and .neg: null distribution
%   corrected_p_matrix  - single matrix: corrected p-values from both pos and neg clusters

[nChannels, nConds] = size(rho_rs);

% Init outputs
obsClusters.pos = cell(nConds,1); obsClusters.neg = cell(nConds,1);
obsClusterStats.pos = cell(nConds,1); obsClusterStats.neg = cell(nConds,1);
correctedPvals.pos = cell(nConds,1); correctedPvals.neg = cell(nConds,1);

maxClusterStatsPerm.pos = zeros(nPerms, nConds);
maxClusterStatsPerm.neg = zeros(nPerms, nConds);

corrected_p_matrix = ones(nChannels, nConds);  % <== Combined matrix

% === Step 1: Observed clusters (pos and neg)
for direction = ["pos", "neg"]
    switch direction
        case "pos"
            sigMaskObs = (p_rs < alpha) & (rho_rs > 0);
        case "neg"
            sigMaskObs = (p_rs < alpha) & (rho_rs < 0);
    end

    for c = 1:nConds
        sigElectrodes = find(sigMaskObs(:, c));
        if isempty(sigElectrodes)
            obsClusters.(direction){c} = {};
            obsClusterStats.(direction){c} = [];
            continue;
        end

        subAdj = adjacency(sigElectrodes, sigElectrodes);
        G = graph(subAdj);
        comp = conncomp(G);
        clusters = arrayfun(@(x) sigElectrodes(comp == x), 1:max(comp), 'UniformOutput', false);
        obsClusters.(direction){c} = clusters;

        clusterStats = cellfun(@(cl) sum(rho_rs(cl, c)), clusters);
        obsClusterStats.(direction){c} = clusterStats;
    end
end

% === Step 2: Permutation loop (shared)
h = waitbar(0, 'Running permutations...');
for perm = 1:nPerms
    waitbar(perm/nPerms, h, sprintf('Permutation %d/%d', perm, nPerms));
    [rho_perm, p_perm] = compute_corr_func(perm);

    for c = 1:nConds
        % Positive
        sigPos = find((p_perm(:, c) < alpha) & (rho_perm(:, c) > 0));
        if ~isempty(sigPos)
            subAdj = adjacency(sigPos, sigPos);
            G = graph(subAdj);
            comp = conncomp(G);
            clusters = arrayfun(@(x) sigPos(comp == x), 1:max(comp), 'UniformOutput', false);
            stats = cellfun(@(cl) sum(rho_perm(cl, c)), clusters);
            maxClusterStatsPerm.pos(perm, c) = max(stats);
        end

        % Negative
        sigNeg = find((p_perm(:, c) < alpha) & (rho_perm(:, c) < 0));
        if ~isempty(sigNeg)
            subAdj = adjacency(sigNeg, sigNeg);
            G = graph(subAdj);
            comp = conncomp(G);
            clusters = arrayfun(@(x) sigNeg(comp == x), 1:max(comp), 'UniformOutput', false);
            stats = cellfun(@(cl) sum(rho_perm(cl, c)), clusters);
            maxClusterStatsPerm.neg(perm, c) = min(stats);
        end
    end
end
close(h);

% === Step 3: Corrected p-values
for direction = ["pos", "neg"]
    for c = 1:nConds
        obsStats = obsClusterStats.(direction){c};
        nClusters = length(obsStats);
        correctedPvals.(direction){c} = zeros(size(obsStats));

        for cl = 1:nClusters
            obsVal = obsStats(cl);
            if direction == "pos"
                pval = mean(maxClusterStatsPerm.pos(:, c) >= obsVal);
            else
                pval = mean(maxClusterStatsPerm.neg(:, c) <= obsVal);
            end

            % Store in struct
            correctedPvals.(direction){c}(cl) = pval;

            % Assign to combined corrected_p_matrix
            cluster_inds = obsClusters.(direction){c}{cl};
            corrected_p_matrix(cluster_inds, c) = pval;
        end
    end
end

end