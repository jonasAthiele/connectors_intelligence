function [rho_mat, p_mat] = compute_partial_spearman(mse, rpm, confounds)
% Computes partial Spearman correlations between MSE and RPM controlling for confounds.
%
% Inputs:
%   mse       - nSubjects x nChannels x nTimescales matrix
%   rpm       - nSubjects x 1 vector
%   confounds - nSubjects x nConfounds matrix
%
% Outputs:
%   rho_mat   - nChannels x nTimescales matrix of correlation coefficients
%   p_mat     - nChannels x nTimescales matrix of p-values

[nSubjects, nChannels, nTimescales] = size(mse);

rho_mat = zeros(nChannels, nTimescales);
p_mat = ones(nChannels, nTimescales);

for ch = 1:nChannels
    for ts = 1:nTimescales
        x = mse(:, ch, ts);
        [r, p] = partialcorr(x, rpm, confounds, 'type', 'Spearman');
        rho_mat(ch, ts) = r;
        p_mat(ch, ts) = p;
    end
end

end