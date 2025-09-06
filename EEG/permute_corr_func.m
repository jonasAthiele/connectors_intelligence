function varargout = permute_corr_func(cmd, mse_in, rpm_in, confounds_in)
    persistent mse rpm confounds

    if nargin > 0 && ischar(cmd) && strcmp(cmd, 'set_data')
        mse = mse_in;
        rpm = rpm_in;
        confounds = confounds_in;
        return;
    end

    permIndex = cmd;
    
    rng(permIndex, 'twister'); % for reproducibility
    permOrder = randperm(length(rpm));
    rpm_perm = rpm(permOrder);

    [rho_perm, p_perm] = compute_partial_spearman(mse, rpm_perm, confounds);

    varargout{1} = rho_perm;
    varargout{2} = p_perm;
end