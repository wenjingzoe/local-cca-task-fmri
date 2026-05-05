
function e = generate_ar1(n,rho,sigma2)
e = zeros(n,1);
e(1) = sqrt(sigma2/(1-rho^2))*randn;
for t = 2:n
    e(t) = rho*e(t-1) + sqrt(sigma2)*randn;
end
end
