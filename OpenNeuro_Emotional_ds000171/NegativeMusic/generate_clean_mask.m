
function [final_mask, idx_mask, idx_non, info] = generate_clean_mask( ...
    mask_file, rc1_file, rc2_file, rc3_file, opts)
% GENERATE_CLEAN_MASK: 生成用于 CCA 的干净脑掩模（保留 GM/WM，排除 CSF 和脑室/脑干等）

% ---------- 默认参数 ----------
if nargin < 5, opts = struct(); end
opts = set_default(opts, 'gm_thr', 0.60);
opts = set_default(opts, 'wm_thr', 0.60);
opts = set_default(opts, 'use_wm', false);
opts = set_default(opts, 'csf_thr', 0.10);
opts = set_default(opts, 'exclude_labels', {'Brain Stem','4th Ventricle','3rd Ventricle', ...
     'Left Lateral Ventricle','Right Lateral Ventricle'});
opts = set_default(opts, 'atlas_name', 'Neuromorphometrics');
opts = set_default(opts, 'min_cluster', 0);

% ---------- 读取数据 ----------
Vmask = spm_vol(mask_file);
refMask = spm_read_vols(Vmask) > 0;      % raw brain mask
rc1 = spm_read_vols(spm_vol(rc1_file));  % GM
rc2 = spm_read_vols(spm_vol(rc2_file));  % WM
rc3 = spm_read_vols(spm_vol(rc3_file));  % CSF

% ---------- 基础掩模 ----------
gm_mask = rc1 > opts.gm_thr;
wm_mask = opts.use_wm * (rc2 > opts.wm_thr);
csf_mask = rc3 > opts.csf_thr;
brain_mask = (gm_mask | wm_mask) & ~csf_mask;
base_mask = refMask & brain_mask;

% ---------- Atlas 屏蔽 ----------
atlas_excl = false(size(base_mask));
try
    A = spm_atlas('load', opts.atlas_name);
    label_names = {A.labels.name};
    for i = 1:numel(opts.exclude_labels)
        label = opts.exclude_labels{i};
        label_idx = find(strcmp(label_names, label));
        if isempty(label_idx)
            warning('[!] Label not found: "%s"', label);
            continue;
        end

        
        label_val = A.labels(label_idx).index;
        Vm = spm_atlas('mask', A, label);
        if isstruct(Vm) && strcmp(Vm.fname, 'Neuromorphometrics_mask.nii')
            Vm.fname = fullfile(spm('Dir'), 'atlas', 'Neuromorphometrics', 'Neuromorphometrics.nii');
        end

        [p, n, e] = fileparts(Vm.fname);
        rfile = fullfile(p, ['r' n e]);

        if ~isfile(rfile)
            spm_reslice(char(Vmask.fname, Vm.fname), struct('which',1,'interp',0,'mean',false));
        end

        Vr = spm_vol(rfile);
        Y  = spm_read_vols(Vr);
        mLbl = (round(Y) == label_val);

        atlas_excl = atlas_excl | mLbl;
        fprintf('[✓] Excluded label "%s": %d voxels masked.\n', label, nnz(mLbl));
    end
catch ME
    warning('[!] Atlas exclusion skipped: %s', ME.message);
end

% ---------- 最终掩模 ----------
mask_no_atlas = base_mask & ~atlas_excl;
final_mask = mask_no_atlas; %3d logical mask

% ---------- 小簇移除 ----------
if opts.min_cluster > 0
    try
        CC = bwconncomp(final_mask, 6);
        sizes = cellfun(@numel, CC.PixelIdxList);
        small = sizes < opts.min_cluster;
        final_mask(cat(1, CC.PixelIdxList{small})) = false;
    catch
        warning('[!] Skipped cluster removal (Image Toolbox required).');
    end
end

% --------------------
idx_mask = find(final_mask);         % linear indice in MIN150
idx_non  = setdiff(1:numel(final_mask), idx_mask);

% ---------- inf0 structure----------
info = struct();
info.opts          = opts;
info.refMask       = refMask;
info.gm_mask       = gm_mask;
info.wm_mask       = wm_mask;
info.csf_mask      = csf_mask;
info.brain_mask    = brain_mask;
info.base_mask     = base_mask;
info.atlas_excl    = atlas_excl;
info.mask_no_atlas = mask_no_atlas;

end

function s = set_default(s, field, val)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = val;
end
end