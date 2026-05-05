
function Y_filtered = get_driftremoved_fMRI(spm_mat_path, nii_path)

    % Load design matrix and filter struct
    load(spm_mat_path, 'SPM');

    % --- Step 1: Read 4D fMRI image ---
    V = spm_vol(nii_path);             % Get header for all volumes
    Y_raw = spm_read_vols(V);          % 4D: [x y z t]
    Y_2D = reshape(Y_raw, [], size(Y_raw, 4));  % [nVoxels x nTime]

    % --- Step 2: High-pass filter (drift removal) ---
    Y_filtered = spm_filter(SPM.xX.K, Y_2D')';  % [nVox x nTime]


end