function CV_score = Compute_FunctionalMassScore(SystemTable, ValidSystems)
% Compute functional mass score:
% sum of voxel counts over valid functional systems

CV_score = 0;

for i = 1:height(SystemTable)

    system_name = SystemTable.FunctionalSystem{i};
    voxel_count = SystemTable.VoxelCount(i);

    % skip invalid systems
    if isempty(system_name)
        continue;
    end

    if any(strcmpi(system_name, ValidSystems))
        CV_score = CV_score + voxel_count;
    end
end

end