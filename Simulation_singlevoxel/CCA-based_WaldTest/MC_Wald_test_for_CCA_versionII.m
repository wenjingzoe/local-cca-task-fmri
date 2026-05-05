%% =========================================================
% MC_verify_Wald_beta_raw.m
% Wald test for raw plug-in CCA estimator version ii--thesis version
% T_n -> chi^2_{rm}
%% Last update: 12-JAN-2026 

%% =========================================================
close all;
clear; clc; rng(1);

%% -----------------------
% Global parameters
% ------------------------
m      = 18;
r      = 1;
p      = r*m;
rho    = 0.5;
sigma2 = 1;
MC     = 500;

n_list = [100 200 400 800];

Nx = 5; Ny = 5;
hs = 1.0;

%% -----------------------
% Spatial grid & weights
% ------------------------
[xg, yg] = meshgrid(1:Nx,1:Ny);
coords_full = [yg(:), xg(:)];
v0 = [3,3];

dist = sqrt(sum((coords_full - v0).^2,2));
w_full = exp(-(dist.^2)/(2*hs^2));
w_full = w_full / sum(w_full);

idx = find(w_full > 1e-3);
w = w_full(idx);
coords = coords_full(idx,:);
Kloc = length(w);

%% -----------------------
% True HRF
% ------------------------
h0 = zeros(m,1);
%h0(1)=1; % H1
%h0=h0/norm(h0);

%% =========================================================
% Loop over n
%% =========================================================
for ii = 1:length(n_list)

    n = n_list(ii);
    fprintf('\nRunning n = %d\n', n);

    %% -----------------------
    % Design matrix
    % ------------------------
    s = block_stimulus(n,10);
    S = zeros(n,m);
    for j = 1:m
        S(j:end,j) = s(1:n-j+1);
    end

    %% Drift & smoother
    t = (1:n)'/n;
    d = 10*sin(pi*(t-0.21));
    Sd = build_local_linear_smoother(t,0.05);
    A = eye(n) - Sd;
    Stilde = A*S;
    sigma22=1/n*transpose(Stilde)*Stilde;
  
    %% -----------------------
    % Storage
    % ------------------------
    beta_raw_hat = zeros(p,MC);
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
        sigma11_hat=sqrt(sigma11_square);
        % sigma21
        sigma21 = zeros(r*m,1);
        for k = 1:Kloc
            sigma21 = sigma21 + w(k) * (1/n) * (Stilde' * ytilde_all(:,k));
        end

       
        beta_raw_hat(:,mc) = (1/sigma11_hat) * (sigma22 \ sigma21);

    end


%% =====================================================
% Wald statistic (H0: beta_raw = 0)
%% =====================================================
Delta = sqrt(n) * beta_raw_hat;          % null hypothesis
Lambda_beta = cov(Delta');               % plug-in covariance

L = chol(Lambda_beta, 'lower');
Z = L \ Delta;
Tn = sum(Z.^2, 1);

%% =====================================================
% Chi-square plot 
%% =====================================================
figure('Position',[100 100 800 600]);
histogram(Tn,50,'Normalization','pdf','FaceAlpha',0.5);
hold on;


x = linspace(0,60,400);
plot(x,chi2pdf(x,p),'r','LineWidth',2);

crit = chi2inv(0.99, p);
size_emp = mean(Tn > crit);

xline(crit, 'r--', 'LineWidth', 2);
xlim([0 60]);

xlabel('$T_n$','Interpreter','latex');
ylabel('Density','Interpreter','latex');

title(sprintf(['Wald statistic $T_n$ under $H_0$, $n=%d$, df=%d\n' ...
               'Empirical size (1\\%% level) = %.3f'], ...
               n, p, size_emp), ...
      'Interpreter','latex');
grid on;

fname = sprintf('Fig_Wald_null_rawCCA_n%d.png', n_list(ii));
exportgraphics(gcf, fname, 'Resolution', 300);

end


%{
 %% =====================================================
 % Wald statistic histogram with empirical power in title (H1)
 %% =====================================================
   
Delta = sqrt(n) * (beta_raw_hat);         % alternative hypothesis
Lambda_beta = cov(Delta');                % plug-in covariance

L = chol(Lambda_beta, 'lower');
Z = L \ Delta;
Tn = sum(Z.^2, 1);
figure('Position',[100 100 800 600]);

% 1. empirical histogram
histogram(Tn, 50, ...
    'Normalization','pdf', ...
    'FaceAlpha',0.5, ...
    'FaceColor',[0.6 0.6 0.9], ...
    'EdgeColor','none');
hold on;

% 2. chi-square pdf
x = linspace(0, max(Tn)*1.1, 400);
plot(x, chi2pdf(x, p), 'r-', 'LineWidth', 2);

% 3. critical value
crit = chi2inv(0.99, p);
xline(crit, 'r--', 'LineWidth', 2,'Color','black');

% 4. empirical power
power_emp = mean(Tn > crit);

% 5. labels
xlabel('$T_n$', 'Interpreter','latex', 'FontSize', 14);
ylabel('Density', 'Interpreter','latex', 'FontSize', 14);

% 6. title (thesis-style)
title(sprintf([ ...
    'Wald statistic under $H_1$, $n=%d$, df=%d\n' ...
    'Empirical rejection probability (1\\%% level) = %.3f' ], ...
    n, p, power_emp), ...
    'Interpreter','latex', ...
    'FontSize', 16);

% 7. save figure (thesis-style filename)
fname = sprintf( ...
    'Fig_Wald_raw_CCA_H1_chi2_df%d_alpha01_n%d.png', ...
    p, n);
exportgraphics(gcf, fname, 'Resolution', 300);   

 

end
%}
