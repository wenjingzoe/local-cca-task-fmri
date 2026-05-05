% =========================================================
% Monte Carlo verification of CLT:
% sqrt(n)(theta_tilde - theta_0) => N(0, Lambda_0)
% Lemma: asymptotic normality of \widetilde{\theta}_0
% Last update: 5-MAY-2026
% =========================================================
close all;
clear; clc; rng(1);

%% -----------------------
% Global parameters
% ------------------------
m      = 18;        % HRF length
r      = 1;         % number of stimuli
p      = r*m + 1;   % dim(theta_tilde)
rho    = 0.5;       % AR(1)
sigma2 = 1;         % noise variance
MC     = 500;       % Monte Carlo runs

n_list = [100 200 400 800];

%% -----------------------
% Spatial grid & kernel
% ------------------------
Nx = 5; Ny = 5;
[xg, yg] = meshgrid(1:Nx,1:Ny);
coords = [yg(:), xg(:)];
v0 = [3,3];
hs = 1.0;

dist = sqrt(sum((coords - v0).^2,2));
w_full = exp(-(dist.^2)/(2*hs^2));
w_full = w_full / sum(w_full);

idx_neigh = find(w_full > 1e-3);
w = w_full(idx_neigh);
Kloc = length(w);

%% -----------------------
% True HRF
% ------------------------
h0 = zeros(m,1);
h0(1) = 1;
h0 = h0 / norm(h0);

%% =========================================================
% Loop over n
% =========================================================
for ii = 1:length(n_list)

    n = n_list(ii);
    fprintf('Running n = %d\n', n);

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
    t = (1:n)'/n;
    d = 10*sin(pi*(t-0.21));

    %% -----------------------
    % Drift smoother
    % ------------------------
    Sd = build_local_linear_smoother(t, 0.05);
    A  = eye(n) - Sd;
    Stilde = A*S;

    %% -----------------------
    % Storage
    % ------------------------
    theta_tilde = zeros(p, MC);

    %% =====================================================
    % Monte Carlo loop
    % =====================================================
    for mc = 1:MC

        %% Generate AR(1) noise
        eps_all = zeros(n, Kloc);
        for k = 1:Kloc
            eps_all(:,k) = generate_ar1(n, rho, sigma2);
        end

        %% Generate data
        y_all = zeros(n, Kloc);
        for k = 1:Kloc
            y_all(:,k) = S*h0 + d + eps_all(:,k);
        end

        %% Drift removal
        ytilde_all = A * y_all;

        %% -----------------------
        % Compute theta_tilde
        % ------------------------

        sigma11_square= (1/n) * sum( w' .* sum(ytilde_all.^2, 1) );

        % sigma21
        sigma21 = zeros(r*m,1);
        for k = 1:Kloc
            sigma21 = sigma21 + w(k) * (1/n) * (Stilde' * ytilde_all(:,k));
        end

        theta_tilde(:,mc) = [sigma11_square; sigma21];

    end

    %% -----------------------
    % Centering & scaling
    % ------------------------
    theta_bar = mean(theta_tilde, 2);
    Delta = sqrt(n) * (theta_tilde - theta_bar);

    %% -----------------------
    % Estimate Lambda_0
    % ------------------------
    Lambda_hat = cov(Delta');

    %% -----------------------
    % Whitening
    % ------------------------
    L = chol(Lambda_hat, 'lower');
    Delta_white = L \ Delta;

    %% =====================================================
    % QQ plots
    % =====================================================
    num_fig = 2;
    comp_per_fig = ceil(p / num_fig);   


  
    for f=1:num_fig
       comp_start = (f-1)*comp_per_fig + 1;
       comp_end   = min(f*comp_per_fig, p); 
       figure('Position',[100 100 1500 400]);
       num_comp = comp_end - comp_start + 1;
       nrow = 2;
       ncol = ceil(num_comp / nrow);

        for j = comp_start:comp_end
        subplot(nrow, ncol, j - comp_start + 1)
        qqplot(Delta_white(j,:));
        xlim([-4 4]);
        ylim([-4 4]);
        axis square;
        xlabel('Standard Normal Quantiles');
        ylabel('Sample Quantiles');
        title(sprintf('Comp %d', j), 'FontSize', 11);
        grid on
        end


       sgtitle(sprintf( ...
        'QQ plots of whitened $\\sqrt{n}(\\tilde{\\theta}-\\bar{\\theta})$, comps %d--%d, n=%d', ...
        comp_start, comp_end, n),  ...
        'Interpreter','latex', 'FontSize', 20);

       fname = sprintf( ...
       'Fig_MC_Verify_CLT_theta_titlde_n%d_comps_%d_%d.png', ...
        n_list(ii), comp_start, comp_end);
    
       exportgraphics(gcf, fname, 'Resolution', 300);

       %close(gcf);
    end

 

end