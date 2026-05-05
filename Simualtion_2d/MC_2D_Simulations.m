% Build fixed activation mask on a 40x40 grid
% 4 blocks of size 2x5 placed at top, bottom, left, right
% -------------------------------------------------
%% =========================================================

close all; 
clear; clc; rng(0);

MC     = 500;       % Monte Carlo runs
Nx     = 40;
Ny     = 40;
active_mask = false(Ny, Nx);   % logical mask (row=y, col=x)

%% -------------------------
% 4 Block size
% -------------------------
h = 2;   % height
w = 5;   % width

%% -------------------------
% TOP block (centered)
% -------------------------
row_top = 5;                          % near top
col_top = 9 : 13;                     % centered
active_mask(row_top:row_top+h-1, col_top) = true;

%% -------------------------
% BOTTOM block (centered)
% -------------------------
row_bot = Ny - h - 4;                 % near bottom
col_bot = 9 : 13;
active_mask(row_bot:row_bot+h-1, col_bot) = true;

%% -------------------------
% LEFT block (centered)
% -------------------------
col_left = 3;                         % near left
row_left = 9 : 10;
active_mask(row_left, col_left:col_left+w-1) = true;

%% -------------------------
% RIGHT block (centered)
% -------------------------
col_right = Nx - w - 6;               % near right
row_right =13:14;
active_mask(row_right, col_right+1:col_right+w) = true;

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

%% ----------------
% data generation  20 x 20 x 200(time pt)

%------------------
%% -----------------------
% Global parameters
% ------------------------
V      = Nx * Ny;    % number of voxels
n      = 200;       % num of time pts
m      = 18;        % HRF length
r      = 1;         % stimulus types
rho    = 0.5;       % AR(1)
sigma2 = 1;         % noise variance
hs     = 1;         % spatial bandwidth
p      =r*m;

%% -----------------------
%% -----------------------
% HRF (true)
% ------------------------
h0 = zeros(m,1);
h0(3:10) = 1;
h0=h0/norm(h0);

%% -----------------------
% Stimulus & design matrix
% ------------------------
s = block_stimulus(n,10);
S = zeros(n,m);

for j = 1:m
    S(j:end,j) = s(1:n-j+1);
end


%% -----------------------
% Drift
% ------------------------
t = (1:n)'/n; d = 10*sin(pi*(t-0.21));
%% -----------------------
% Drift smoother
% ------------------------
Sd = build_local_linear_smoother(t, 0.1);
A  = eye(n) - Sd;
Stilde = A * S;
dtilde = A * d;


%% -----------------------
% Spatial neighborhood size
% ------------------------
nvx = 3;   % neighborhood width  (x direction)
nvy = 3;   % neighborhood height (y direction)

half_x = floor(nvx/2);   % = 2
half_y = floor(nvy/2);   % = 2

%% =====================================================
% Monte Carlo loop
% =====================================================

Method = {'Zhang K_{bc}'; 'Proposed T_n (Chi-square)'; 'Proposed T_n (Empirical)'};
n_method=numel(Method);

TPR_MC=zeros(MC,n_method);
FPR_MC=zeros(MC,n_method);
Dice_MC=zeros(MC,n_method);

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

%%step1:
% threshold (top 5% as pseudo-active)
thr = quantile(task_energy, 0.95);
pseudo_mask     =zeros(V,1);
pseudo_active   = find(task_energy > thr);    % linear indices
pseudo_inactive = find(task_energy <= thr);
pseudo_mask(pseudo_active)=1;

% pseudo-active voxel coordinates
[y_act, x_act] = ind2sub([Ny, Nx], pseudo_active);
pseudo_active_2d = [y_act, x_act];   % N_act x 2  [row, col]

% pseudo-inactive voxel coordinates
[y_inact, x_inact] = ind2sub([Ny, Nx], pseudo_inactive);
pseudo_inactive_2d = [y_inact, x_inact];



%% =====================================================
% Loop over all voxels
% =====================================================
beta_raw_hat_all = zeros(p,V);    % p = dim(beta), V = # of voxels
beta_raw_norm_all= zeros(V,1);

h_0_all=zeros(p, V);              % p = dim(beta)
h_0_norm_all= zeros(V,1);


h_hat_bc_zhang = zeros(p, V);
K_bc=zeros(V,1);

Gamma = (Stilde'*Stilde)/n;
Sigma22 = Gamma+0.1*eye(p);


for v0 = 1:V

% -----------------------------
%  (row=y, col=x)
% -----------------------------
[y0, x0] = ind2sub([Ny, Nx], v0);

% -----------------------------
% neighbour range of v0
% -----------------------------
x_range = max(1, x0-half_x) : min(Nx, x0+half_x); 
y_range = max(1, y0-half_y) : min(Ny, y0+half_y);

% -----------------------------
% neighbour 2d-coordinate
% -----------------------------
[YY, XX] = ndgrid(y_range, x_range);      % row-major

neigh_coords = [YY(:), XX(:)];            % [y, x]=[row, col]

% -----------------------------
% convert to linear index
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

    %% Feasible GLS for h
    [h_hat,~] = estimate_h_FGLS(Stilde,ybar);
    h_0_all(:,v0)=h_hat;   % p = dim(beta)
    h_0_norm_all(v0)=vecnorm(h_hat);

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
    beta_raw_norm_all(v0) = vecnorm( beta_raw_hat_all(v0) );


    %% zhang's method
    yv = Ytilde(:, v0);          % drift-removed signal
    % GLS
    [h_hat_zhang, R_hat] = estimate_h_FGLS(Stilde, yv);
    % drift estimate
    d_hat       = Sd * (yv - S*h_hat_zhang);
    d_hat_tilde = A * d_hat;

   % residuals
    r_hat=yv-Stilde*h_hat_zhang;
    r_hat_bc=r_hat-d_hat_tilde;
    % matrices
    V_tilde=Stilde' * (R_hat \ Stilde);   % rm x rm

    h_hat_bc=h_hat_zhang - (V_tilde \ (Stilde' * (R_hat \ d_hat_tilde)));
    h_hat_bc_zhang(:,v0)=h_hat_bc;

    % numerator
    num1 = (h_hat_bc)'*(V_tilde)*h_hat_bc;

   % denominator (with DOF correction!)
    den1 = (r_hat_bc' * (R_hat \ r_hat_bc)) / (n - r*m);

    K_bc(v0) = num1 / den1;


end


%% =====================================================
% Estimate null covariance from pseudo-inactive voxels (*)
% =====================================================

num_pick = floor(length(pseudo_inactive)/2);
perm = randperm(length(pseudo_inactive), num_pick);
inactive_pick = pseudo_inactive(perm);

% null mean (beta0)
beta0_hat = mean(beta_raw_hat_all(:, inactive_pick), 2);   % p x 1

% scaled deviations (time-scale!)
Delta_inact = sqrt(n) * ...
    (beta_raw_hat_all(:, inactive_pick) - beta0_hat);      % p x Npick

% covariance estimate
Lambda_inact = cov(Delta_inact');    % p x p

% =====================================================
% Wald statistic for all voxels
% =====================================================

Delta_all = sqrt(n) * ...
    (beta_raw_hat_all - beta0_hat);   % p x V

L = chol(Lambda_inact, 'lower');
Z = L \ Delta_all;                    % p x V

Tn = sum(Z.^2, 1);                    % 1 x V
Tn_map = reshape(Tn, Ny, Nx);

% define interior mask (same radius as neighborhood)
edge_margin_x = half_x;
edge_margin_y = half_y;

interior_mask = false(Ny, Nx);
interior_mask( ...
    (1+edge_margin_y):(Ny-edge_margin_y), ...
    (1+edge_margin_x):(Nx-edge_margin_x) ) = true;

interior_idx = find(interior_mask(:));

Tn_interior = Tn(interior_idx);

thr = quantile(Tn(inactive_pick), 0.99);  % empirical null
rej = false(V,1);
rej(interior_idx) = Tn(interior_idx) > thr;
rej_map = reshape(rej, Ny, Nx);


%use chi_square(rm) at 1% signifcant level
thr_theory = chi2inv(0.99, p);   % p = 18
rej_chi_vec= zeros(V,1);
rej_chi_vec(Tn>thr_theory)=1;
rej_chi_map = reshape(rej_chi_vec, Ny, Nx);


% =====================================================
%% zhang 
% =====================================================

thr = chi2inv(0.99, r*m);    % Theorem 4.2
rej_zhang = reshape(K_bc > thr, Ny, Nx);



%====================================================
%% metric
%=====================================================


metrics_zhang = compute_TPR_FPR(rej_zhang, ...
                                active_mask, ...
                                interior_mask);

metrics_chi   = compute_TPR_FPR(rej_chi_map, ...
                                active_mask, ...
                                interior_mask);

metrics_emp   = compute_TPR_FPR(rej_map, ...
                                active_mask, ...
                                interior_mask);

TPR_MC(mc,:)  = [metrics_zhang.TPR ; 
                metrics_chi.TPR; 
                metrics_emp.TPR]';

FPR_MC(mc,:)  = [metrics_zhang.FPR;
        metrics_chi.FPR;
        metrics_emp.FPR]';

Dice_MC(mc,:) = [metrics_zhang.Dice;
        metrics_chi.Dice;
        metrics_emp.Dice]';

end

TPR_mean  = mean(TPR_MC, 1);
FPR_mean  = mean(FPR_MC, 1);
Dice_mean = mean(Dice_MC, 1);


ResultsTable = table( ...
    TPR_mean(:), ...
    FPR_mean(:), ...
    Dice_mean(:), ...
    'VariableNames', {'TPR','FPR','Dice'}, ...
    'RowNames', Method);

disp(ResultsTable);

