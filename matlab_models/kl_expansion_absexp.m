% KL expansion function with squared exponential covariance
function [phi, lambda] = kl_expansion_absexp(N, variance, x)
    phi = zeros(length(x), N);
    lambda = zeros(N, 1);
    L = max(x) - min(x); % Domain length
    for i = 1:N
        for j = 1:length(x)
            phi(j, i) = sin(i * pi * x(j) / L);
        end
        lambda(i) = variance * exp(-i^2 * pi^2 / (2 * L^2));  % Eigenvalues for squared exponential covariance
    end
end