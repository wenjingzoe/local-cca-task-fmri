
function Sd = build_local_linear_smoother(t, b)

% Local linear smoothing matrix S_d (Zhang 2008 style)
% t : n x 1 time points (e.g. t_i = i/n)
% b : bandwidth

n = length(t);
Sd = zeros(n,n);

for i = 1:n
    
    % Differences
    u = (t - t(i)) / b;
    % Gaussian kernel weights
    K = exp(-0.5*u.^2) / sqrt(2*pi);
    w = K / b;   % scaling by 1/b    
    % Design matrix for local linear fit
    X = [ones(n,1), t - t(i)];
    
    % Weight matrix
    W = diag(w);
    
    % Local linear projection
    XtWX = X' * W * X;
    
    % Numerical safety
    if rcond(XtWX) < 1e-10
        Sd(i,:) = 0;
    else
        A = [1 0] / XtWX;
        Sd(i,:) = A * X' * W;
    end
    
end
end
