function [RegionTable, SystemTable] = AtlasSystemFromActivation( ...
            Active_mask_3D, Affmat, atlas, funcmap_file)
% =========================================================
% AtlasSystemFromActivation
%
% Map active voxels to atlas regions, then aggregate voxel
% counts to functional systems (multi-label allowed).
%
% INPUTS:
%   Active_mask_3D : [Nx x Ny x Nz] logical activation mask
%   nii_path       : NIfTI path (for MNI transform)
%   atlas          : spm_atlas object (e.g. Neuromorphometrics)
%   funcmap_file   : Excel file mapping Region -> System(s)
%
% OUTPUTS:
%   RegionTable : table with columns
%       - RegionName
%       - VoxelCount
%       - FunctionalSystem
%
%   SystemTable : table with columns
%       - FunctionalSystem
%       - VoxelCount
% =========================================================

%% ---------------------------------------------------------
% Load functional-system mapping
% ---------------------------------------------------------
funcmap = readtable(funcmap_file);
label_names   = funcmap.LabelName;
label_systems = funcmap.FunctionalSystem;

SystemList = unique( ...
    strtrim(strsplit(strjoin(label_systems, ';'), ';')) );
SystemList = SystemList(~cellfun(@isempty, SystemList));

%% ---------------------------------------------------------
% Get active voxel indices
% ---------------------------------------------------------
[idx_x, idx_y, idx_z] = ind2sub(size(Active_mask_3D), find(Active_mask_3D));
NumActive = numel(idx_x);


%% ---------------------------------------------------------
% Query atlas for each active voxel
% ---------------------------------------------------------
all_region_names = cell(NumActive,1);

for k = 1:NumActive

    mni_coords = Affmat * [idx_x(k); idx_y(k); idx_z(k); 1];

    xY = struct( ...
        'xyz', mni_coords(1:3), ...
        'def', 'sphere', ...
        'spec', 1);

    region = spm_atlas('query', atlas, xY);

    if isempty(region)
        all_region_names{k} = 'Unknown';
    else
        all_region_names{k} = strjoin(region, ', ');
    end
end

%% ---------------------------------------------------------
% Count voxel frequency per region
% ---------------------------------------------------------
[regions_unique, ~, ic] = unique(all_region_names);
region_counts = accumarray(ic, 1);

% sort by frequency
[region_counts_sorted, sort_idx] = sort(region_counts, 'descend');
regions_sorted = regions_unique(sort_idx);

%% ---------------------------------------------------------
% Build RegionTable + accumulate system counts
% ---------------------------------------------------------
region_rows = {};
SystemVoxelCount = zeros(numel(SystemList), 1);

for i = 1:numel(regions_sorted)

    region = regions_sorted{i};
    count  = region_counts_sorted(i);

    % match region -> functional system(s)
    match_idx = strcmpi(label_names, region);

    if any(match_idx)
        systems_raw = label_systems{find(match_idx,1)};
        matched_systems = strtrim(strsplit(systems_raw, ';'));
    else
        matched_systems = {'Unknown'};
    end

    % accumulate system voxel counts (THIS is your logic)
    for s = 1:numel(SystemList)
        if any(strcmpi(matched_systems, SystemList{s}))
            SystemVoxelCount(s) = SystemVoxelCount(s) + count;
        end
    end

    region_rows(end+1,:) = { ...
        region, ...
        count, ...
        strjoin(matched_systems, ', ') ...
    };
end

RegionTable = cell2table(region_rows, ...
    'VariableNames', {'RegionName','VoxelCount','FunctionalSystem'});

%% ---------------------------------------------------------
% Build SystemTable
% ---------------------------------------------------------
SystemTable = table( ...
    SystemList(:), ...
    SystemVoxelCount(:), ...
    'VariableNames', {'FunctionalSystem','VoxelCount'});

end