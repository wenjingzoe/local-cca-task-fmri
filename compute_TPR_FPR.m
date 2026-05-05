function metrics = compute_TPR_FPR(rej_map, true_active_mask, interior_mask)

    % only evaluate interior voxels
    idx = interior_mask(:) == 1;

    R = rej_map(:);
    A = true_active_mask(:);

    R = R(idx);
    A = A(idx);

    TP = sum(R == 1 & A == 1);
    FP = sum(R == 1 & A == 0);
    FN = sum(R == 0 & A == 1);
    TN = sum(R == 0 & A == 0);

    metrics.TP  = TP;
    metrics.FP  = FP;
    metrics.FN  = FN;
    metrics.TN  = TN;

    metrics.TPR = TP / (TP + FN);
    metrics.FPR = FP / (FP + TN);

    % additional useful metrics
    metrics.Precision = TP / (TP + FP);
    metrics.Dice      = 2*TP / (2*TP + FP + FN);
    metrics.IoU       = TP / (TP + FP + FN);
end