
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
