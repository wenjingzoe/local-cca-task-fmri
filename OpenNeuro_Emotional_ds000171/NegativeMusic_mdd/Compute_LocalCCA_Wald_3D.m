
function [active_id, Act_map_3D, Tn_map_3D, ...
          Tn_map_1D, Tn_mask, beta_raw_hat_all] = ...
    Compute_LocalCCA_Wald_3D(Ytilde, Stilde, Sigma22, ...
                             final_mask, idx_mask, pseudo_inactive, ...
                             hs, n)
% =
% =========================================================
% Compute_LocalCCA_Wald_3D
%
% Compute a 3D local CCA–based Wald test statistic for
% task-fMRI data using spatially weighted second-order
% moment estimation.
%
% This function implements the LOCAL CCA inference step
% at a single spatial scale (bandwidth hs), treating each
% voxel as the center of a 3×3×3 spatial neighborhood.
%
% ---------------------------------------------------------
% MODEL OVERVIEW
% ---------------------------------------------------------
% At each voxel v0, local second-order moments are estimated
% using neighboring voxels v in a 3×3×3 cube:
%
%   sigma_11(v0) = (1/n) ∑_v w(v,v0) || y_v ||^2
%   sigma_21(v0) = (1/n) ∑_v w(v,v0) S̃ᵀ y_v
%
% where:
%   - y_v are drift-removed fMRI time series
%   - S̃ is the drift-removed design matrix
%   - w(v,v0) are Gaussian spatial kernel weights
%
% The raw local CCA estimator is:
%
%   \hat{\beta}_raw(v0) = sigma_11(v0)^(-1/2) Σ_22^{-1} sigma_21(v0)
%
% A Wald statistic is then constructed using an empirical
% null covariance estimated from pseudo-inactive voxels.
%
% ---------------------------------------------------------
% INPUTS
% ---------------------------------------------------------
% Ytilde          : [n × Vmask]
%   Drift-removed and (optionally) prewhitened fMRI data,
%   restricted to voxels inside the brain mask.
%
% Stilde          : [n × m]
%   Drift-removed design matrix (m = HRF length in scans).
%
% Sigma22         : [m × m]
%   Regularized second-moment matrix of the design:
%     Sigma22 = (1/n)\widetide{S}^T \widetide{S} + 0.01*I
%
% final_mask      : [Nx × Ny × Nz] (logical)
%   Brain mask in 3D volume space.
%
% idx_mask        : [Vmask × 1]
%   Linear indices of voxels inside final_mask.
%
% pseudo_inactive : vector of indices in {1,…,Vmask}
%   Voxels assumed to be inactive, used to estimate the
%   null covariance of β̂_raw.
%
% hs              : scalar
%   Spatial kernel bandwidth (in voxel units).
%
% n               : scalar
%   Number of time points.
%
% ---------------------------------------------------------
% OUTPUTS
% ---------------------------------------------------------
% active_id       : vector (linear indices)
%   Linear indices (full-volume) of detected active voxels.
%
% Act_map_3D      : [Nx × Ny × Nz] (logical)
%   Binary activation map after Wald thresholding.
%
% Tn_map_3D       : [Nx × Ny × Nz]
%   Wald statistic mapped back to 3D volume.
%
% Tn_map_1D       : [Nx*Ny*Nz × 1]
%   Wald statistic in full linear index space.
%
% Tn_mask         : [1 × Vmask]
%   Wald statistic restricted to mask voxels.
%
% beta_raw_hat_all: [m × Vmask]
%   Raw local CCA estimators at each masked voxel.
%
% =========================================================
% ---------------------------------------------------------
% Dimensions
% ---------------------------------------------------------
[Nx, Ny, Nz] = size(final_mask);
Vmask = numel(idx_mask);
m = size(Stilde, 2);
p = m;   % r = 1

% ---------------------------------------------------------
% Precompute 3D coordinates of mask voxels
% ---------------------------------------------------------
[x_all, y_all, z_all] = ind2sub([Nx, Ny, Nz], idx_mask);
coords_mask = [x_all, y_all, z_all];   % [Vmask x 3]

% ---------------------------------------------------------
% Map: full linear index -> mask index
% ---------------------------------------------------------
full2mask = zeros(Nx*Ny*Nz, 1);
full2mask(idx_mask) = 1:Vmask;

% ---------------------------------------------------------
% Storage
% ---------------------------------------------------------
beta_raw_hat_all = zeros(p, Vmask);

% ---------------------------------------------------------
% Loop over mask voxels
% ---------------------------------------------------------
for iv = 1:Vmask

    % Center voxel (full-space coordinates)
    x0 = coords_mask(iv,1);
    y0 = coords_mask(iv,2);
    z0 = coords_mask(iv,3);

    % ---- 3x3x3 neighborhood in full space ----
    x_range = max(1,x0-1):min(Nx,x0+1);
    y_range = max(1,y0-1):min(Ny,y0+1);
    z_range = max(1,z0-1):min(Nz,z0+1);

    [XX,YY,ZZ] = ndgrid(x_range, y_range, z_range);
    neigh_full = sub2ind([Nx,Ny,Nz], XX(:), YY(:), ZZ(:));

    % ---- keep only mask voxels ----
    neigh_full = neigh_full(final_mask(neigh_full));
    neigh_mask = full2mask(neigh_full);   % mask-space indices

    Kloc = numel(neigh_mask);

    % ---- spatial distances ----
    neigh_coords = coords_mask(neigh_mask, :);
    dists = sqrt( ...
        (neigh_coords(:,1)-x0).^2 + ...
        (neigh_coords(:,2)-y0).^2 + ...
        (neigh_coords(:,3)-z0).^2 );

    % ---- Gaussian kernel weights ----
    w = exp(-(dists.^2)/(2*hs^2));
    w = w / sum(w);

    % ---- local data ----
    Yloc = Ytilde(:, neigh_mask);   % [n x Kloc]

    % ---- sigma11 ----
    sigma11 = (1/n) * sum( w' .* sum(Yloc.^2, 1) );

    % ---- sigma21 ----
    sigma21 = zeros(p,1);
    for k = 1:Kloc
        sigma21 = sigma21 + w(k) * (Stilde' * Yloc(:,k));
    end
    sigma21 = sigma21 / n;

    % ---- raw CCA estimator ----
    beta_raw_hat_all(:,iv) = sigma11^(-1/2) * (Sigma22 \ sigma21);
end

% ---------------------------------------------------------
% Estimate null covariance (mask space)
% ---------------------------------------------------------
 num_pick = floor(length(pseudo_inactive)/2);
 perm = randperm(length(pseudo_inactive), num_pick); %**
 inactive_pick = pseudo_inactive(perm); %**

 beta0_hat = mean(beta_raw_hat_all(:, inactive_pick), 2);

 Delta_inact = sqrt(n) * ...
    (beta_raw_hat_all(:, inactive_pick) - beta0_hat);

Lambda_inact = cov(Delta_inact');

% numerical stabilization
Lambda_inact = (Lambda_inact + Lambda_inact')/2 + 1e-6*eye(p);

% ---------------------------------------------------------
% Wald statistic
% ---------------------------------------------------------
Delta_all = sqrt(n) * (beta_raw_hat_all - beta0_hat);
L = chol(Lambda_inact, 'lower');
Z = L \ Delta_all;

%% =========================================================
% Wald statistic in mask space
% =========================================================

% Z: [p x Vmask] whitened local CCA scores
Tn_mask = sum(Z.^2, 1);          % [1 x Vmask], Wald statistic per masked voxel

% Empirical threshold (e.g. top 5%)
thr = quantile(Tn_mask, 0.99);

% Detected active voxels (mask-space indices)
maskactive_id = find(Tn_mask > thr);

% Convert back to full-brain linear indices
active_id = idx_mask(maskactive_id);

%% =========================================================
% Map Wald statistic back to 3D brain volume
% =========================================================

% Initialize full-brain 1D map
Tn_map_1D = zeros(Nx * Ny * Nz, 1);

% Fill in mask voxels only
Tn_map_1D(idx_mask) = Tn_mask;

% Reshape to 3D volume
Tn_map_3D = reshape(Tn_map_1D, [Nx, Ny, Nz]);

%% =========================================================
% Binary activation map (thresholded)
% =========================================================

Act_map_1D = false(Nx * Ny * Nz, 1);
Act_map_1D(active_id) = true;

Act_map_3D = reshape(Act_map_1D, [Nx, Ny, Nz]);

end

