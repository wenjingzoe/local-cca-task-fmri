function [h_hat, Rhat] = estimate_h_FGLS(Stilde, y)
% =========================================================
% estimate_h_FGLS
% Two-step Feasible GLS estimator for h
% Assumes AR(1) noise
%
% INPUT:
%   Stilde : n x m  (design after drift removal)
%   y      : n x 1  (spatially weighted signal)
%
% OUTPUT:
%   h_hat  : m x 1
% =========================================================

[n,~] = size(Stilde);

%% -------- Step 1: OLS --------
h_ols = (Stilde' * Stilde) \ (Stilde' * y);
res   = y - Stilde * h_ols;

%% -------- Estimate rho (AR(1)) --------
rho_hat = (res(2:end)' * res(1:end-1)) / ...
          (res(1:end-1)' * res(1:end-1));
rho_hat = max(min(rho_hat,0.95),-0.95);

%% -------- Build Rhat --------
Rhat = toeplitz(rho_hat.^(0:n-1));

%% -------- Step 2: GLS --------
h_hat1 = (Stilde' * (Rhat \ Stilde)) \ ...
        (Stilde' * (Rhat \ y));
res   = y - Stilde * h_hat1;
%% -------- Estimate rho (AR(1)) --------
rho_hat = (res(2:end)' * res(1:end-1)) / ...
          (res(1:end-1)' * res(1:end-1));
rho_hat = max(min(rho_hat,0.95),-0.95);

%% -------- Build Rhat --------
Rhat = toeplitz(rho_hat.^(0:n-1));

%% -------- Step 2: GLS --------
h_hat = (Stilde' * (Rhat \ Stilde)) \ ...
        (Stilde' * (Rhat \ y));
err=vecnorm(h_hat-h_hat1);


while err>10e-4

    h_hat_pre=h_hat;
    res   = y - Stilde * h_hat_pre;
    %-------- Estimate rho (AR(1)) --------
    rho_hat = (res(2:end)' * res(1:end-1)) / ...
          (res(1:end-1)' * res(1:end-1));
    rho_hat = max(min(rho_hat,0.95),-0.95);
    Rhat = toeplitz(rho_hat.^(0:n-1));

    h_hat = (Stilde' * (Rhat \ Stilde)) \ ...
        (Stilde' * (Rhat \ y));
    err=vecnorm(h_hat-h_hat_pre);
end

%disp(err);

end