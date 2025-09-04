% % KL expansion function with squared exponential covariance

function [phi_KL, lambda_KL] = kl_expansion_SE(N_KL,sigma2, L, x)
% Squared exponential covariance function
N = length(x);
% C = zeros(N, N);

C = sigma2 .* exp(-(x-x').^2 ./ (2 .* L^2));

%% other ways of generating covariance functions
% cov_func = @(x1,x2) 2 .* exp(-abs(x1 - x2) / L);
% C3 = zeros(N,N);
% for i=1:N
%     for j=1:N
%         C3(i,j) = sigma2 .* exp(-(x(i) - x(j)).^2 /(2 .* L^2)); %option 1
%         C3(i,j) = cov_func(x(i),x(j)); %option 2

%
%     end
% end

%%
% Eigen decomposition of the covariance matrix
[V, D] = eig(C);
[lambda, idx] = sort(diag(D), 'descend'); % Sort eigenvalues in descending order
phi = V(:, idx);                              % Corresponding eigenvectors

lambda_KL = lambda(1:N_KL);
phi_KL = phi(:, 1:N_KL);
end
