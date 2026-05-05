function [CV_score_list, best_h, SystemTable_by_h] = CV_LocalCCA_Wald_3D( ...
    Ytilde, Stilde, Sigma22, final_mask, idx_mask, ...
    pseudo_inactive_all, hs_list, n, B, ...
    Affmat, atlas, funcmap_file)
% =========================================================
% CV_LocalCCA_Wald_3D
%
% Cross-validation procedure for selecting the spatial
% bandwidth h_s in the 3D local CCA-based Wald test,
% performed separately for each subject.
%
% ---------------------------------------------------------
% OVERVIEW
% ---------------------------------------------------------
% For each candidate spatial bandwidth h_s in hs_list, this
% function:
%
%   (1) Repeats a bootstrap-based Wald test B times,
%       where each bootstrap randomly subsamples the
%       pseudo-inactive voxel set to estimate the null
%       covariance of the local CCA estimator.
%
%   (2) For each bootstrap realization:
%       - Computes the local CCA-based Wald statistic
%         using a 3x3x3 spatial neighborhood and Gaussian
%         kernel with bandwidth h_s;
%       - Thresholds the Wald statistic to obtain a binary
%         activation map;
%       - Maps detected active voxels to anatomical regions
%         and functional systems using an atlas;
%       - Accumulates voxel counts per predefined functional
%         system.
%
%   (3) Defines a *functional mass CV score* for each
%       bootstrap as the total number of detected voxels
%       belonging to valid functional systems.
%
%   (4) Averages the functional mass scores across bootstrap
%       runs to obtain a CV score for each h_s.
%
% The optimal spatial bandwidth h_s is selected as the one
% that maximizes the average functional mass score.
%
% ---------------------------------------------------------
% INPUTS
% ---------------------------------------------------------
% Ytilde : [n x Vmask]
%   Drift-removed and prewhitened fMRI time series restricted
%   to brain mask voxels.
%
% Stilde : [n x m]
%   Drift-removed stimulus design matrix (after temporal
%   smoothing removal).
%
% Sigma22 : [m x m]
%   Regularized second-moment matrix of the design regressors.
%
% final_mask : [Nx x Ny x Nz] (logical)
%   3D brain mask in image space.
%
% idx_mask : [Vmask x 1]
%   Linear indices of voxels inside the brain mask.
%
% pseudo_inactive_all : [K x 1]
%   Mask-space indices (subset of 1:Vmask) corresponding to
%   voxels pre-screened as pseudo-inactive.
%
% hs_list : [H x 1]
%   List of candidate spatial bandwidths to be cross-validated.
%
% n : scalar
%   Number of time points.
%
% B : scalar
%   Number of bootstrap repetitions used for CV.
%
% Affmat : [4 x 4]
%   Affine matrix mapping voxel indices to MNI coordinates.
%
% atlas : struct
%   Loaded SPM atlas structure (e.g., Neuromorphometrics).
%
% funcmap_file : string
%   Path to Excel file mapping anatomical region names to
%   functional system labels (possibly multi-label).
%
% ---------------------------------------------------------
% OUTPUTS
% ---------------------------------------------------------
% CV_score_list : [H x 1]
%   Average functional mass CV score for each candidate h_s.
%
% best_h : scalar
%   Selected spatial bandwidth maximizing CV_score_list.
%
% SystemTable_by_h : {H x 1} cell array
%   Each cell contains a table summarizing voxel counts per
%   functional system for a given h_s.
%
%   For h_s = hs_list(id_h), SystemTable_by_h{id_h} has:
%       - Rows   : predefined functional systems
%       - Columns:
%           * Boot1 ... BootB  : voxel counts per bootstrap
%           * MeanVoxelCount  : average voxel count across B
%
% ---------------------------------------------------------
% FUNCTIONAL SYSTEMS USED
% ---------------------------------------------------------
% The CV score only counts voxels assigned to the following
% predefined functional systems:
%
%   {'Attention', 'EmotionControl', 'LanguageCognition',
%    'Auditory', 'Visual', 'Motor', 'Memory'}
%
% Voxels labeled as 'Unknown' or belonging to other systems
% are excluded from the CV score.
%
% ---------------------------------------------------------
% NOTES
% ---------------------------------------------------------
% - CV is performed independently for each subject.
% - The bootstrap resampling targets the null covariance
%   estimation step of the Wald test.
% - This procedure selects h_s based on stability and
%   functional interpretability rather than likelihood.
%
% =========================================================

ValidSystems = { ...
    'Attention', ...
    'EmotionControl', ...
    'LanguageCognition', ...
    'Auditory', ...
    'Visual', ...
    'Motor', ...
    'Memory' };

numSys = numel(ValidSystems);
num_h  = numel(hs_list);

CV_score_list    = zeros(num_h,1);
SystemTable_by_h = cell(num_h,1);      

% =========================================================
% Loop over h
% =========================================================
for id_h = 1:num_h

    fprintf('=== current h: %.3f ===\n', hs_list(id_h));
    hs = hs_list(id_h);

    CV_score_vec = zeros(B,1);

    % -----------------------------------------------------
    % Store voxel mass per system per bootstrap
    % size: [numSys x B]
    % -----------------------------------------------------
    SystemVoxelCount_boot = zeros(numSys, B);   %%%% NEW

    % =====================================================
    % Bootstrap loop
    % =====================================================
    for b = 1:B

        rng(b);   % reproducible CV

        % ---------------------------------------------
        % Subsample pseudo-inactive voxels
        % ---------------------------------------------
        num_pick = floor(numel(pseudo_inactive_all)/2);
        perm = randperm(numel(pseudo_inactive_all), num_pick);
        pseudo_inactive = pseudo_inactive_all(perm);

        % ---------------------------------------------
        % Local CCA + Wald
        % ---------------------------------------------
        [~, Act_map_3D] = Compute_LocalCCA_Wald_3D( ...
            Ytilde, Stilde, Sigma22, ...
            final_mask, idx_mask, pseudo_inactive, ...
            hs, n);

        % ---------------------------------------------
        % Atlas → System table
        % ---------------------------------------------
        [~, SystemTable] = AtlasSystemFromActivation( ...
            Act_map_3D, Affmat, atlas, funcmap_file);

        % ---------------------------------------------
        % Accumulate voxel mass per functional system
        % ---------------------------------------------
        for s = 1:numSys
            sys_name = ValidSystems{s};
            idx_sys = strcmpi(SystemTable.FunctionalSystem, sys_name);
            SystemVoxelCount_boot(s, b) = sum(SystemTable.VoxelCount(idx_sys));
        end

        % ---------------------------------------------
        % CV score (your functional mass score)
        % ---------------------------------------------
        CV_score_vec(b) = sum(SystemVoxelCount_boot(:, b));

        fprintf('  Bootstrap %d: CV score = %.2f\n', b, CV_score_vec(b));

    end

    % =====================================================
    % Aggregate for this h
    % =====================================================
    CV_score_list(id_h) = mean(CV_score_vec);

    % -----------------------------------------------------
    % Build SystemTable for this h
    % -----------------------------------------------------
    T = table;
    T.FunctionalSystem = ValidSystems(:);

    for b = 1:B
        T.(sprintf('Boot%d', b)) = SystemVoxelCount_boot(:, b);
    end

    T.MeanVoxelCount = mean(SystemVoxelCount_boot, 2);

    SystemTable_by_h{id_h} = T;   %%%% NEW

end

% =========================================================
% Select best h
% =========================================================
[~, best_idx] = max(CV_score_list);
best_h = hs_list(best_idx);

fprintf('=== Best h = %.3f (CV score = %.2f) ===\n', ...
        best_h, CV_score_list(best_idx));

end