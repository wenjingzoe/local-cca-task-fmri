clear; close all; clc;

addpath('/Users/wenjingy/Desktop/OpenNeuro_Emotional_ds000171/spm_25');
spm('defaults','fmri');
spm_jobman('initcfg');

atlas = spm_atlas('load','Neuromorphometrics');

%% =====================================================
% Paths & global settings
%% =====================================================
root = '/Users/wenjingy/Desktop/OpenNeuro_Emotional_ds000171/NegativeMusic_mdd';

%sub_list = {'sub-mdd02','sub-mdd03','sub-mdd05', ...
 %           'sub-mdd07','sub-mdd09','sub-mdd14', ...
  %          'sub-mdd15','sub-mdd18'};


sub_list = {'sub-mdd02'};

trainrun = 'trainrun';
funcmap_file = 'Atlas_FunctionMap_Multilabel.xlsx';

hs = 1;        % spatial bandwidth
m  = 12;       % HRF length (scans)
b  = 0.1;      % temporal smoothing bandwidt
hs_list=0.5:0.5:1;

%region_xlsx = fullfile(root,'AllSubjects_Region.xlsx');
%system_xlsx = fullfile(root,'AllSubjects_System.xlsx');

%% =====================================================
% Build stimulus & design matrix (shared)
%% =====================================================
n_design = 105;          % fixed by dataset
TR       = 3;            % seconds
t_sec    = (0:n_design-1) * TR;

events = readtable('NegativeStart_task_events.tsv','FileType','text');

s = zeros(n_design,1);
for e = 1:height(events)
    if ismember(events.trial_type{e},{'negative_music','positive_music'})
        onset = events.onset(e);
        dur   = events.duration(e);
        s( (t_sec>=onset) & (t_sec<onset+dur) ) = 1;
    end
end

S = zeros(n_design,m);
for j = 1:m
    S(j:end,j) = s(1:n_design-j+1);
end

t  = (1:n_design)' / n_design;
Sd = build_local_linear_smoother(t,b);
A  = eye(n_design) - Sd;
Stilde = A*S;

Gamma   = (Stilde'*Stilde) / n_design;
Sigma22 = Gamma + 0.1 * eye(m);
P_S     = Stilde * ((Stilde'*Stilde)\Stilde');

%% =====================================================
% Loop over subjects
%% =====================================================
for sj = 1:numel(sub_list)

    subj = sub_list{sj};
    fprintf('\n================ %s ================\n',subj);

    subj_dir = fullfile(root,subj);
    func_dir = fullfile(subj_dir,'func');
    anat_dir = fullfile(subj_dir,'anat');
    spm_path = fullfile(subj_dir,'first_level_analysis','SPM.mat');
    out_xlsx = fullfile(root, sprintf('CV_SystemTables_%s_train.xlsx', subj));

    files = spm_select('FPList',func_dir, ...
        ['^swar' subj '_task-music_' trainrun '_bold\.nii$']);

    if isempty(files)
        warning('No fMRI found for %s',subj);
        continue;
    end

    nii_path = strtrim(files(1,:));
    V = spm_vol(nii_path);
     Affmat = V(1).mat;   % 4x4 affine matrix

    load(spm_path,'SPM');

    %% ================================================
    % Drift removal + whitening
    %% ================================================
    Y = get_driftremoved_fMRI(spm_path,nii_path);   % [V x T]
    Vw = sqrtm(inv(full(SPM.xVi.V)));
    Y  = Y * Vw;

    n = size(Y,2);   % actual time points (should be 105)

    %% ================================================
    % Mask
    %% ================================================
    mask_file = fullfile(SPM.swd,'mask.nii');

    rc1 = reslice_to_mask(mask_file,anat_dir,subj,'c1');
    rc2 = reslice_to_mask(mask_file,anat_dir,subj,'c2');
    rc3 = reslice_to_mask(mask_file,anat_dir,subj,'c3');

    opts.exclude_labels = {'Brain Stem','4th Ventricle','3rd Ventricle', ...
        'Left Lateral Ventricle','Right Lateral Ventricle'};

    [final_mask,idx_mask] = generate_clean_mask( ...
        mask_file,rc1,rc2,rc3,opts);

    vol_size = size(final_mask);
    num_vox  = prod(vol_size);

    %% ================================================
    % Mask-space data
    %% ================================================
    Y4D = nan([vol_size,n]);
    for tt = 1:n
        tmp = zeros(num_vox,1);
        tmp(idx_mask) = Y(idx_mask,tt);
        Y4D(:,:,:,tt) = reshape(tmp,vol_size);
    end

    Y1D    = reshape(Y4D,[],n)';      % [n x V]
    Ymask = Y1D(:,idx_mask);          % [n x Vmask]
    Ytilde = A * Ymask;

    %% ================================================
    % Pseudo-active voxels (energy screening)
    %% ================================================
    task_energy = mean((P_S * Ytilde).^2,1);
    thr = quantile(task_energy,0.95);
    pseudo_active_id  = find(task_energy > thr);
    pseudo_inactive_id= find(task_energy < thr);

    B=5;
 

   [CV_score_list, best_h, SystemTable_by_h] = CV_LocalCCA_Wald_3D( ...
    Ytilde, Stilde, Sigma22, final_mask, idx_mask, ...
    pseudo_inactive_id, hs_list, n, B, ...
    Affmat, atlas, funcmap_file);

   num_h = numel(hs_list);

for id_h = 1:num_h

    % ---------------------------------------------
    % Get table for this h
    % ---------------------------------------------
    T = SystemTable_by_h{id_h};

    % ---------------------------------------------
    % Add subject column
    % ---------------------------------------------
    T.Subject = repmat({subj}, height(T), 1);

    % ---------------------------------------------
    % Safe sheet name (no dots)
    % e.g. h_1p000 instead of h_1.000
    % ---------------------------------------------
    h_val = hs_list(id_h);
    sheet_name = sprintf('h_%s', strrep(num2str(h_val,'%.3f'), '.', 'p'));

    % ---------------------------------------------
    % Write to Excel
    % ---------------------------------------------
    writetable(T, out_xlsx, ...
        'Sheet', sheet_name, ...
        'WriteMode', 'overwrite');

end

fprintf('✓ CV (train) system tables saved to:\n  %s\n', out_xlsx);


    %{
    %% ================================================
    % Local CCA + Wald
    %% ================================================
    [~, Act_map_3D] = ...
        Compute_LocalCCA_Wald_3D( ...
            Ytilde,Stilde,Sigma22, ...
            final_mask,idx_mask, ...
            pseudo_inactive_id, ...
            hs,n);

    %% ================================================
    % Atlas → Region + System summary
    %% ================================================
    [RegionTable,SystemTable] = AtlasSystemFromActivation( ...
        Act_map_3D, V(1).mat, atlas, funcmap_file);

    RegionTable.Subject = repmat({subj},height(RegionTable),1);
    SystemTable.Subject = repmat({subj},height(SystemTable),1);

    %% ================================================
    % Write Excel (ONE SHEET PER SUBJECT)
    %% ================================================
    writetable(RegionTable,region_xlsx, ...
        'Sheet',subj,'WriteMode','overwrite');

    writetable(SystemTable,system_xlsx, ...
        'Sheet',subj,'WriteMode','overwrite');

    %}


end


%% ================== Helper Function ==================

function rc_file = reslice_to_mask(mask_file, anat_dir, subj, tissue_tag)
    % Reslice tissue probability map to match fMRI mask space
    native_file = fullfile(anat_dir, sprintf('%s%s_T1w.nii', tissue_tag, subj));
    rc_file = fullfile(anat_dir, sprintf('r%s%s_T1w.nii', tissue_tag, subj));
    if ~isfile(native_file)
        error('Missing file: %s', native_file);
    end
    spm_reslice(char(mask_file, native_file), struct('which',1,'interp',1,'mean',false));
    if ~isfile(rc_file)
        error('Resliced %s not found: %s', tissue_tag, rc_file);
    end
end

