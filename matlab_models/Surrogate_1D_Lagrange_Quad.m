% Surrogate_1D_Lagrange.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: April, 2025
%
clear; clc; close all;
%% Try to build a surrogate for the following function
% using Lagrange polynomials
f_HF = @(q) (6.*q-2).^2 .* sin(12.*q - 4);
q_space = sort(unifrnd(0,1,1000,1));
rho_q = 1; %
MC = 20;
q_data_MC = unifrnd(0,1,MC,1);
[q_data_GL, w_GL] = lgwt(8, 0, 1);

f_HF_mean = mean(f_HF(q_space));
y_data_MC = f_HF(q_data_MC);
y_data_GL = f_HF(q_data_GL);

f_MC_mean = sum(y_data_MC)./M;
f_quad_mean = sum(y_data_GL.*w_GL.*rho_q);

%% Now build a lagrange polynomial
M = length(y_data_MC);
lag = @(x) 0;
for i=1:M
    f = make_lagrange(q_data_MC,i,M);
    lag = @(x) lag(x) + y_data_MC(i).*f(x);
end
f_LF = lag;
mean_MC = mean(f_LC);
%% Buil
%% Evaluate the LF and HF models
figure(1);clf; hold on;
plot(q_space,f_HF(q_space),'k','LineWidth',3);
plot(q_space,f_LF(q_space),'--r','LineWidth',3);
plot(q_data_MC,y_data_MC,'ko');
grid on; set(gca,'FontSize',20);
legend('HF','LF','Data')
%%
function f = make_lagrange(q_points,m_curr,M)
%% Make a recursive function call for Lagrange Polynomials
%
% q_points = vector of q values from your sampling
% m_curr   = number of values in the q_vector
% q        = current parameter value you are interested in

ids = 1:M;
ids = ids(ids~=m_curr);
q_eval = q_points(ids);
q_m    = q_points(m_curr);
f = @(q,q_eval) 1;
for i=1:M-1
    f = @(q,q_eval) f(q,q_eval) .* (q - q_eval(i))./(q_m - q_eval(i));
end
f = @(q) f(q,q_eval);
end


% lgwt subfunction for calculating points and weights
% Written by Greg von Winckel - 02/25/2004
%https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
function [x, w] = lgwt(N, a, b)
    % lgwt computes the points and weights for Gaussian-Legendre quadrature
    %
    % Input:
    %   N - number of points
    %   a - lower limit of the interval
    %   b - upper limit of the interval
    %
    % Output:
    %   x - quadrature points
    %   w - weights associated with the points

    x = zeros(N, 1); % Preallocate points
    w = zeros(N, 1); % Preallocate weights

    % Initial setup
    n = (1:N)';
    beta = 0.5 ./ sqrt(1 - (2*n).^(-2)); % Recurrence coefficients
    T = diag(beta, 1) + diag(beta, -1); % Tridiagonal matrix

    % Compute eigenvalues (roots) and weights
    [V, D] = eig(T);
    x = diag(D); % Eigenvalues are the roots
    w = 2 * V(1, :).^2'; % Weights

    % Scale roots and weights to the interval [a, b]
    x = a + (x + 1) * (b - a) / 2; % Map roots to [a, b]
    w = w * (b - a) / 2; % Scale weights to [a, b]
end