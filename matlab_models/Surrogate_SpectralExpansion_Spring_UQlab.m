%% Mitchel J. Colebank
%  Project 5: Problem 3
clear; clc; close all;
%%
N_MC = 5e3;      % Number of Monte Carlo samples
M = 8;          % Number of Gauss-Hermite quadrature points used to approximate discrete projection
n_poly = 3;
n_param = 2;
%%
% Initial conditions
z0 = 2; v0 = -2;
IC = [z0;v0];

% Two parameter spring model with no forcing term
tend = 20;
tspace = linspace(0,tend,20);
f_HF = @(q) min(call_spring_model(q,IC,tspace));

%%
% Input the parameter means and standard deviations.
qbar = [1 1];
UB = [5 5];
LB = [0 0];
%%
% Compute N_MC samples for each of the parameters and construct the grid for omega_F.
q_MC = zeros(N_MC,2);
q_MC(:,1) = unifrnd(LB(1),UB(1),N_MC,1);
q_MC(:,2) = unifrnd(LB(2),UB(2),N_MC,1);

%%
% Input the M=10 Gauss-Legendre quadrature points and weights from a table.
[points_leg, weights_leg] = gaussLegendreQuadrature(M);
weights_leg = 0.5.*weights_leg;
% points_leg = [-0.9739 -0.86506 -0.6794 -0.4334 -0.14887 0.14887 0.4334 0.6794 0.86506 0.9739]';
% weights_leg = 0.5*[0.06667 0.14945 0.21908 0.2693 0.2955 0.2955 0.2693 0.21908 0.14945 0.06667]';
% Get corresponding polynomials (note that normalization factors are
% already included here
P0 = @(i) 1;
P1 = @(i) points_leg(i);
P2 = @(i) (0.5)*(3.*points_leg(i).^2 - 1);
P3 = @(i) (0.5)*(5.*points_leg(i).^3 - 3.*points_leg(i));
h0 = 1;
h1 = 1/3;
h2 = 1/5;
h3 = 1/7;
% p_func = @(i) [P0(i) P1(i) P2(i) P3(i)];
% % Make multiindex
% multi_cell

%% We scale the uniform [-1,1] to our upper and lower bound
q_leg = 0.5.*(UB-LB).*points_leg + 0.5*(UB+LB);
% We will need nchoosk(poly_order + num_param, num_param) = 10 here
K = nchoosek(n_poly+n_param, n_poly);           % Number of tensored Hermite basis functions Psi(q)
PCE_coeff = zeros(1,K);

%% Construct MC estimates
y_MC = zeros(1,N_MC);
for i=1:N_MC
    y_MC(i) = f_HF(q_MC(i,:));
end

%% Spectral surrogate
counter = 0;
for j1 = 1:M
    for j2 = 1:M
        q_eval = [q_leg(j1,1) q_leg(j2,2)];
        yval = f_HF(q_eval);
        PCE_coeff(1)  = PCE_coeff(1) + yval*P0(j1)*P0(j2)*weights_leg(j1)*weights_leg(j2)./(h0*h0);
        PCE_coeff(2)  = PCE_coeff(2) + yval*P1(j1)*P0(j2)*weights_leg(j1)*weights_leg(j2)./(h1*h0);
        PCE_coeff(3)  = PCE_coeff(3) + yval*P2(j1)*P0(j2)*weights_leg(j1)*weights_leg(j2)./(h2*h0);
        PCE_coeff(4)  = PCE_coeff(4) + yval*P3(j1)*P0(j2)*weights_leg(j1)*weights_leg(j2)./(h3*h0);
        PCE_coeff(5)  = PCE_coeff(5) + yval*P0(j1)*P1(j2)*weights_leg(j1)*weights_leg(j2)./(h0*h1);
        PCE_coeff(6)  = PCE_coeff(6) + yval*P1(j1)*P1(j2)*weights_leg(j1)*weights_leg(j2)./(h1*h1);
        PCE_coeff(7)  = PCE_coeff(7) + yval*P2(j1)*P1(j2)*weights_leg(j1)*weights_leg(j2)./(h2*h1);
        PCE_coeff(8)  = PCE_coeff(8) + yval*P0(j1)*P2(j2)*weights_leg(j1)*weights_leg(j2)./(h0*h2);
        PCE_coeff(9)  = PCE_coeff(9) + yval*P1(j1)*P2(j2)*weights_leg(j1)*weights_leg(j2)./(h1*h2);
        PCE_coeff(10) = PCE_coeff(10) + yval*P0(j1)*P3(j2)*weights_leg(j1)*weights_leg(j2)./(h0*h3);
        counter = counter+1;
    end
end

%%

f = build_PCE_surrogate(PCE_coeff,UB,LB);
%% Show the means and variances
gamma = [h0*h0 h1*h0 h2*h0 h3*h0 h0*h1 h1*h1 h2*h1 h0*h2 h1*h2 h0*h3];
MC_mean = mean(y_MC);
MC_var  = var(y_MC);
PCE_mean = PCE_coeff(1);
PCE_var = sum((PCE_coeff(2:end).^2) .* gamma(2:end));
format bank
disp('MC Estimate & Number Eval')
disp([MC_mean MC_var N_MC]);
disp('PCE Estimate')
disp([PCE_mean PCE_var counter]);


%%

function f = build_PCE_surrogate(coeff,UB,LB)
q_leg = @(x) 0.5.*(UB-LB).*x + 0.5*(UB+LB);
P0 = @(x) 1;
P1 = @(x) q_leg(x);
P2 = @(x) (0.5)*(3.*q_leg(x).^2 - 1);
P3 = @(x) (0.5)*(5.*q_leg(x).^3 - 3.*q_leg(x));

Pfunc_vec = @(x,i) (i==0)*P0(x)+(i==1)*P1(x)+(i==2)*P2(x)+(i==3)*P3(x);

f = @(x) 0;
max_poly = 3;
for i=0:max_poly
    for j=0:max_poly
        if i+j<=max_poly
            f = @(x) f(x) + coeff(i+j+1).*Pfunc_vec(x,j).*Pfunc_vec(x,i);
        end
    end

end
end


function y = call_spring_model(q,IC,tspace)
[~,y] = ode45(@spring_model_2par,tspace,IC,[],q);
y = y(:,1)';
end

function dy = spring_model_2par(t,y,params)

% Unpack parameters
C      = params(1);
K      = params(2);

% Redefine state variables for convenience
z = y(1);
v = y(2);


% RHS equations
accel_eq = (- C.*v - K.*z);
dy = [v;
    accel_eq];

end



%%
function [points, weights] = gaussLegendreQuadrature(n)
% Generate Gauss-Legendre quadrature points and weights
% Input:
%   n - Number of quadrature points
% Output:
%   points - Quadrature points
%   weights - Quadrature weights

% Initialization
points = zeros(n, 1);
weights = zeros(n, 1);

% Initial guesses for the roots (using Chebyshev nodes)
m = (1:n)';
initial_guesses = cos(pi * (m - 0.25) / (n + 0.5));

% Compute roots using Newton-Raphson method
for i = 1:n
    x = initial_guesses(i); % Initial guess
    for iter = 1:500 % Iteration limit for convergence
        % Evaluate the Legendre polynomial and its derivative
        [P, dP] = legendrePolynomial(x, n);
        % Update guess using Newton's method
        dx = -P / dP;
        x = x + dx;
        % Check for convergence
        if abs(dx) < 1e-15
            break;
        end
    end
    points(i) = x;
end

% Compute weights
for i = 1:n
    [~, dP] = legendrePolynomial(points(i), n);
    weights(i) = 2 / ((1 - points(i)^2) * (dP^2));
end

% Sort points and weights (optional, for consistency)
[points, order] = sort(points);
weights = weights(order);
end

function [P, dP] = legendrePolynomial(x, n)
% Compute Legendre polynomial P_n(x) and its derivative P'_n(x)
% Input:
%   x - Evaluation point
%   n - Polynomial degree
% Output:
%   P - Value of P_n(x)
%   dP - Value of P'_n(x)

P0 = 1; % P_0(x)
P1 = x; % P_1(x)

if n == 0
    P = P0;
    dP = 0;
    return;
elseif n == 1
    P = P1;
    dP = P0;
    return;
end

for k = 2:n
    P2 = ((2*k - 1) * x * P1 - (k - 1) * P0) / k; % Recurrence relation
    P0 = P1;
    P1 = P2;
end

P = P1;
dP = n * (P0 - x * P1) / (1 - x^2); % Derivative using recurrence relation
end