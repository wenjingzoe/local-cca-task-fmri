clear; close all; clc;

%% =====================================================
% Add paths & initialize SPM
%% =====================================================
% ===== Spartan paths =====
project_root = '/home/wenjingy/SpatialLocalCCAForTask-related-fMRI';
addpath(genpath(project_root));
addpath(fullfile(project_root, ...
    'OpenNeuro_Emotional_ds000171', 'spm_25'));

root = fullfile(project_root, ...
    'OpenNeuro_Emotional_ds000171', 'NegativeMusic');

spm('defaults','fmri');
spm_jobman('initcfg');

atlas = spm_atlas('load','Neuromorphometrics');



sub_list = { ...
    'sub-control01', 'sub-control04', 'sub-control06', ...
    'sub-control07', 'sub-control09', 'sub-control11', ...
    'sub-control15', 'sub-control16', 'sub-control17', 'sub-control18' };

trainrun     = 'trainrun';
funcmap_file = 'Atlas_FunctionMap_Multilabel.xlsx';

% ---- Local CCA hyperparameters ----
hs_list = 0.5 : 0.25 : 3.0;    % spatial bandwidth candidates
B       = 10;                % number of bootstrap replicates

% ---- Design / temporal parameters ----
m  = 12;     % HRF length (in scans)
b  = 0.1;    % temporal smoothing bandwidth
TR = 3;      % seconds
n_design = 105;

%% =====================================================
% Build stimulus & design matrix (shared across subjects)
%% =====================================================
fprintf('Building design matrix...\n');

t_sec = (0:n_design-1) * TR;
events = readtable( ...
    fullfile(root, 'NegativeStart_task_events.tsv'), ...
    'FileType','text');

s = zeros(n_design,1);
for e = 1:height(events)
    if ismember(events.trial_type{e}, {'negative_music','positive_music'})
        onset = events.onset(e);
        dur   = events.duration(e);
        s( (t_sec >= onset) & (t_sec < onset + dur) ) = 1;
    end
end

% Convolution-style design matrix
S = zeros(n_design, m);
for j = 1:m
    S(j:end, j) = s(1:n_design-j+1);
end

% Drift removal operator
t  = (1:n_design)' / n_design;
Sd = build_local_linear_smoother(t, b);
A  = eye(n_design) - Sd;

Stilde = A * S;

Gamma   = (Stilde' * Stilde) / n_design;
Sigma22 = Gamma + 0.1 * eye(m);
P_S     = Stilde * ((Stilde' * Stilde) \ Stilde');


best_h_results = nan(numel(sub_list),1);

%% =====================================================
% Loop over subjects (TRAINING)
%% =====================================================
for sj = 1:numel(sub_list)

    subj = sub_list{sj};
    fprintf('\n================ TRAINING: %s ================\n', subj);

    subj_dir = fullfile(root, subj);
    func_dir = fullfile(subj_dir, 'func');
    anat_dir = fullfile(subj_dir, 'anat');
    spm_path = fullfile(subj_dir, 'first_level_analysis', 'SPM.mat');

    out_xlsx = fullfile(root, ...
        sprintf('CV_SystemTables_%s_train.xlsx', subj));

    %% -------------------------------------------------
    % Load fMRI data
    %% -------------------------------------------------
    files = spm_select('FPList', func_dir, ...
        ['^swar' subj '_task-music_' trainrun '_bold\.nii$']);

    if isempty(files)
        warning('No training fMRI found for %s. Skipping.', subj);
        continue;
    end

    nii_path = strtrim(files(1,:));
    V = spm_vol(nii_path);
    Affmat = V(1).mat;

    load(spm_path, 'SPM');

    %% -------------------------------------------------
    % Drift removal & prewhitening
    %% -------------------------------------------------
    Y = get_driftremoved_fMRI(spm_path, nii_path);   % [V x T]
    Vw = sqrtm(inv(full(SPM.xVi.V)));
    Y  = Y * Vw;

    n = size(Y,2);   % should be 105

    %% -------------------------------------------------
    % Brain mask
    %% -------------------------------------------------
    mask_file = fullfile(SPM.swd, 'mask.nii');

    rc1 = reslice_to_mask(mask_file, anat_dir, subj, 'c1');
    rc2 = reslice_to_mask(mask_file, anat_dir, subj, 'c2');
    rc3 = reslice_to_mask(mask_file, anat_dir, subj, 'c3');

    opts.exclude_labels = { ...
        'Brain Stem', '4th Ventricle', '3rd Ventricle', ...
        'Left Lateral Ventricle', 'Right Lateral Ventricle' };

    [final_mask, idx_mask] = generate_clean_mask( ...
        mask_file, rc1, rc2, rc3, opts);

    vol_size = size(final_mask);
    num_vox  = prod(vol_size);

    %% -------------------------------------------------
    % Rebuild mask-space data
    %% -------------------------------------------------
    Y4D = nan([vol_size, n]);
    for tt = 1:n
        tmp = zeros(num_vox,1);
        tmp(idx_mask) = Y(idx_mask, tt);
        Y4D(:,:,:,tt) = reshape(tmp, vol_size);
    end

    Y1D    = reshape(Y4D, [], n)';
    Ymask = Y1D(:, idx_mask);
    Ytilde = A * Ymask;

    %% -------------------------------------------------
    % Pseudo-inactive voxels (energy screening)
    %% -------------------------------------------------
    task_energy = mean((P_S * Ytilde).^2, 1);
    thr = quantile(task_energy, 0.95);

    pseudo_active_id   = find(task_energy > thr);
    pseudo_inactive_id = find(task_energy <= thr);

    %% -------------------------------------------------
    % Cross-validation over hs
    %% -------------------------------------------------
    fprintf('Running CV over spatial bandwidths...\n');

    [CV_score_list, best_h, SystemTable_by_h] = CV_LocalCCA_Wald_3D(Ytilde, Stilde, Sigma22, ...
            final_mask, idx_mask, ...
            pseudo_inactive_id, ...
            hs_list, n, B, ...
            Affmat, atlas, funcmap_file);

    best_h_results(sj) = best_h;
    %% -------------------------------------------------
    % Save CV results (one sheet per h)
    %% -------------------------------------------------
    for id_h = 1:numel(hs_list)

        T = SystemTable_by_h{id_h};

        % add subject & h info
        T.Subject = repmat({subj}, height(T), 1);
        T.h = repmat(hs_list(id_h), height(T), 1);

        sheet_name = sprintf('h_%s', ...
            strrep(num2str(hs_list(id_h),'%.3f'), '.', 'p'));

        writetable(T, out_xlsx, ...
            'Sheet', sheet_name, ...
            'WriteMode', 'overwrite');
    end

    fprintf('✓ Training CV finished for %s\n', subj);
    fprintf('  Best h = %.3f\n', best_h);
    fprintf('  Saved to: %s\n', out_xlsx);

end

fprintf('\n=== ALL SUBJECT TRAINING FINISHED ===\n');

 %% =====================================================
 % Save best h results with date
 %% =====================================================
 date_str = datestr(now, 'yyyy-mm-dd_HHMM');
 save_file = fullfile(root,sprintf('Best_h_results_%s.mat', date_str));
 save(save_file, 'best_h_results', 'sub_list', 'hs_list');
 fprintf('\nSaved best_h_results to:\n%s\n', save_file);
