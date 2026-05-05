%% =========================================================
% MC_2D_NeighborhoodSensitivity_LocalCCA.m
%
% Sensitivity analysis of proposed local CCA-based Wald test
% w.r.t. spatial neighborhood size (n_y, n_x)
%
% Empirical-quantile-based Wald rejection
%
% Author: Wenjing ; last updated: 13-JAN-2026-correct version 
% =========================================================
close all; clear; clc; rng(0);
%% -------------------------
% Monte Carlo settings
% -------------------------
MC = 500;

%% -------------------------
% Spatial grid & activation mask
% -------------------------
Nx = 40; Ny = 40;
V  = Nx * Ny;

active_mask = false(Ny, Nx);

% ---- Four activation blocks (2 x 5) ----
h = 2; w = 5;

% top
active_mask(5:6, 9:13) = true;
% bottom
active_mask(Ny-h-3:Ny-h-2, 9:13) = true;
% left
active_mask(9:10, 3:7) = true;
% right
active_mask(13:14, Nx-w-5:Nx-1) = true;

%% -------------------------
% Linear indices
% -------------------------
active_idx   = find(active_mask(:));
inactive_idx = find(~active_mask(:));

figure;
imagesc(active_mask);
axis image;
colormap(gray);
colorbar;
title('True activation mask ');

%% -------------------------
% Data-generation parameters
% -------------------------
n      = 200;
m      = 18;
r      = 1;
p      = r*m;
rho    = 0.5;
sigma2 = 1;
hs     = 1;     % spatial bandwidth
%% -------------------------
% HRF
% -------------------------
h0 = zeros(m,1);
h0(3:10) = 1;
h0 = h0 / norm(h0);

%% -------------------------
% Stimulus & design
% -------------------------
s = block_stimulus(n,10);
S = zeros(n,m);
for j = 1:m
    S(j:end,j) = s(1:n-j+1);
end

%% -------------------------
% Drift & smoothing
% -------------------------
t  = (1:n)'/n;
d  = 10*sin(pi*(t-0.21));
Sd = build_local_linear_smoother(t, 0.1);
A  = eye(n) - Sd;
Stilde = A * S;

Gamma   = (Stilde' * Stilde) / n;
Sigma22 = Gamma + 0.1 * eye(p);

%% =========================================================
% Neighborhood sizes to test
%% =========================================================
neigh_list = [ ...
    1 1;
    3 3;
    5 5;
    7 7 ];

n_neigh = size(neigh_list,1);

TPR_all_quantile   = zeros(n_neigh,1);
FPR_all_quantile   = zeros(n_neigh,1);
Dice_all_quantile  = zeros(n_neigh,1);

TPR_all_chi  = zeros(n_neigh,1);
FPR_all_chi  = zeros(n_neigh,1);
Dice_all_chi = zeros(n_neigh,1);

%% =========================================================
% COMMON interior mask (based on largest neighborhood)
%% =========================================================

edge_margin_y = floor(7/2);
edge_margin_x = floor(7/2);

interior_mask = false(Ny, Nx);
interior_mask( ...
    (1+edge_margin_y):(Ny-edge_margin_y), ...
    (1+edge_margin_x):(Nx-edge_margin_x) ) = true;

interior_idx = find(interior_mask(:));


%% =========================================================
% Loop over neighborhood sizes
%% =========================================================
for ii = 1:n_neigh

    ny = neigh_list(ii,1);
    nx = neigh_list(ii,2);

    fprintf('\nNeighborhood size = (%d,%d)\n', ny, nx);

    half_y = floor(ny/2);
    half_x = floor(nx/2);

    TPR_MC_quantile  =zeros(MC,1);
    FPR_MC_quantile  =zeros(MC,1);
    Dice_MC_quantile =zeros(MC,1);

    TPR_MC_chi  = zeros(MC,1);
    FPR_MC_chi  = zeros(MC,1);
    Dice_MC_chi = zeros(MC,1);


    for mc = 1:MC

     fprintf('Running mc = %d\n', mc);


    %% -----------------------
    % Generate AR(1) noise for all voxels
    % ------------------------
    eps_all = zeros(n, V);
    
     for v = 1:V
       eps_all(:,v) = generate_ar1(n, rho, sigma2);
     end

     eps_tilde = A * eps_all;
    
     %% -----------------------
     % Generate voxelwise data
     % ------------------------
     Y = zeros(n,V);   % each column = one voxel time series
     Y(:,active_idx)=S*h0+ d + eps_all(:,active_idx);
     Y(:,inactive_idx)=d+eps_all(:,inactive_idx);

     Ytilde = A * Y;   % drift removed data
     P_S = Stilde * ((Stilde' * Stilde) \ Stilde');   % n x n
     task_energy = zeros(1, V);

     for v = 1:V
       yv = Ytilde(:,v);              % drift removed signal
       task_energy(v) = mean((P_S * yv).^2);  
     end

     E_map = reshape(task_energy, Ny, Nx);


     % threshold (top 5% as pseudo-active)
     thr = quantile(task_energy, 0.95);
     pseudo_active   = find(task_energy > thr);    % linear indices
     pseudo_inactive = find(task_energy <= thr);
   

     % pseudo-active voxel coordinates
     [y_act, x_act] = ind2sub([Ny, Nx], pseudo_active);
     pseudo_active_2d = [y_act, x_act];   % N_act x 2  [row, col]

     % pseudo-inactive voxel coordinates
     [y_inact, x_inact] = ind2sub([Ny, Nx], pseudo_inactive);
     pseudo_inactive_2d = [y_inact, x_inact];

     
    %% =====================================================
    % Loop over all voxels
    % =====================================================
     beta_raw_hat_all = zeros(p,V);   % p = dim(beta)


    for v0 = 1:V

    % -----------------------------
    %  (row=y, col=x)
    % -----------------------------
    [y0, x0] = ind2sub([Ny, Nx], v0);

    % -----------------------------
    % 截断邻域（自动处理边界）
    % -----------------------------
    x_range = max(1, x0-half_x) : min(Nx, x0+half_x); 
    y_range = max(1, y0-half_y) : min(Ny, y0+half_y);

    % -----------------------------
    % 所有邻域 2D 坐标
    % -----------------------------
    [YY, XX] = ndgrid(y_range, x_range);      % row-major

    neigh_coords = [YY(:), XX(:)];            % [y, x]=[row, col]

    % -----------------------------
    % 转成 linear index（关键一步）
    % -----------------------------
    neigh_idx = sub2ind([Ny, Nx], YY(:), XX(:));

    Kloc = numel(neigh_idx);

   % -----------------------------
   % Gaussian spatial weights
   % -----------------------------
    dists = sqrt( ...
    (neigh_coords(:,1) - y0).^2 + ...
    (neigh_coords(:,2) - x0).^2 );

    w = exp(-(dists.^2) / (2*hs^2));
 
   % normalize
    w = w / sum(w);
    % -----------------------------
    % Collect Spatially weighted signal and noise
    % -----------------------------
  
    epstilde_all= eps_tilde(:,neigh_idx);
    ybar = Ytilde(:,neigh_idx) * w; % [n x 1]

   %% compute Sigma
    Yloc=Ytilde(:,neigh_idx); % [ n x Kloc]
    sigma11 = (1/n) * sum( w' .* sum(Yloc.^2, 1) );
    sigma21 = zeros(p,1);

    for k = 1:Kloc
       sigma21 = sigma21 + w(k) * (Stilde' * Yloc(:,k));
    end
   
    sigma21 = sigma21 / n;
    beta_raw_hat=sigma11^(-1/2) * (Sigma22 \ sigma21);
    beta_raw_hat_all(:,v0) = beta_raw_hat;
  


   end



   %% =====================================================
   % Estimate null covariance from pseudo-inactive voxels
   % =====================================================

   num_pick = floor(length(pseudo_inactive)/2);
   perm = randperm(length(pseudo_inactive), num_pick);
   inactive_pick = pseudo_inactive(perm);

   % null mean (beta0)
   beta0_hat = mean(beta_raw_hat_all(:, inactive_pick), 2);   % p x 1

   % scaled deviations (time-scale!)
   Delta_inact = sqrt(n) * ...
    (beta_raw_hat_all(:, inactive_pick) - beta0_hat);      % p x Npick

   %covariance estimate
   Lambda_inact = cov(Delta_inact');    % p x p

   % =====================================================
   % Wald statistic for all voxels
   % =====================================================

   Delta_all = sqrt(n) * ...
    (beta_raw_hat_all - beta0_hat);   % p x V

   L = chol(Lambda_inact, 'lower');
   Z = L \ Delta_all;                    % p x V

   Tn = sum(Z.^2, 1);                    % 1 x V
   
   % =====================================================
   % quantile--empirical method
   % =====================================================

   thr = quantile(Tn(inactive_pick), 0.99);  % empirical null
   rej = false(V,1);
   rej(interior_idx) = Tn(interior_idx) > thr;
   rej_map = reshape(rej, Ny, Nx);

   % =====================================================
   %  chi_square(rm)-therotical method
   % =====================================================
   thr_theory = chi2inv(0.99, p);   % p = 18
   rej_chi_vec= zeros(V,1);
   rej_chi_vec(Tn>thr_theory)=1;
   rej_chi_map = reshape(rej_chi_vec, Ny, Nx);

   %====================================================
   %% metric
   %=====================================================
    metrics_chi   = compute_TPR_FPR(rej_chi_map, ...
                                active_mask, ...
                                interior_mask);
    TPR_MC_chi(mc)  = metrics_chi.TPR;
    FPR_MC_chi(mc)  = metrics_chi.FPR;
    Dice_MC_chi(mc) = metrics_chi.Dice;


    metrics_emp   = compute_TPR_FPR(rej_map, ...
                                active_mask, ...
                                interior_mask);

    TPR_MC_quantile(mc)  = metrics_emp.TPR;
    FPR_MC_quantile(mc)  = metrics_emp.FPR;
    Dice_MC_quantile(mc) = metrics_emp.Dice;
     
  end

    %% -------------------------
    % Average over MC
    % -------------------------

    TPR_all_quantile(ii)  = mean(TPR_MC_quantile);
    FPR_all_quantile(ii)  = mean(FPR_MC_quantile);
    Dice_all_quantile(ii) = mean(Dice_MC_quantile);

    TPR_all_chi(ii)  = mean(TPR_MC_chi);
    FPR_all_chi(ii)  = mean(FPR_MC_chi);
    Dice_all_chi(ii) = mean(Dice_MC_chi);


end


%% =========================================================
% Results tables
%% =========================================================

% ---- Empirical (quantile-based) Wald test ----
Results_empirical = table( ...
    neigh_list(:,1), neigh_list(:,2), ...
    TPR_all_quantile, FPR_all_quantile, Dice_all_quantile, ...
    'VariableNames', {'n_y','n_x','TPR','FPR','Dice'});

disp('Empirical (quantile-based) Wald test results:');
disp(Results_empirical);

% ---- Theoretical (chi-square-based) Wald test ----
Results_theoretical = table( ...
    neigh_list(:,1), neigh_list(:,2), ...
    TPR_all_chi, FPR_all_chi, Dice_all_chi, ...
    'VariableNames', {'n_y','n_x','TPR','FPR','Dice'});

disp('Theoretical (chi-square-based) Wald test results:');
disp(Results_theoretical);