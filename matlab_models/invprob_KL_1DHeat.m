clear; clc; %close all;
rng(900);
N_kl_terms = 19;
kl_L = .15;
meas_noise = 0.1;
%%
T = 5; % End time for numerical solution
L = 1; % Length of the 1D element
dx = 0.05; % Spatial discretization
dt = 1e-3; % Temporal discretization
x_domain = 0:dx:L; % x points
t_domain = 0:dt:T; % t points
Nx = length(x_domain); % Number of x points
Nt = length(t_domain); % Number of t points

% initial and boundary conditions (feel free to make more complicated)
IC = 0.1.*ones(Nx,1);
BC_0 = 4; % BC at x=0
BC_L = 5; % BC at x=L

%% Thermal conductivity parameter
% Generate true diffusivity using KL expansion
alpha_mean = 0.01;
alpha_variance = 0.0001;
% [phi_true, lambda_true] = kl_expansion_SE(N_kl_terms, alpha_variance, x_domain);
[phi_true, lambda_true] = kl_expansion_SE(N_kl_terms, alpha_variance, kl_L, x_domain);
%%
coeff_true = randn(N_kl_terms, 1);  % Random coefficients for true diffusivity
alpha_true = abs(cos(x_domain.*pi)).*alpha_mean/2;% + phi_true * (sqrt(lambda_true) .* coeff_true);
%%

stab_number = max(alpha_true).*dt./(dx.^2);
if stab_number >= 0.5
    warning('Numerical Instability may occur.')
end
u_true = heateq_1D(alpha_true,t_domain,x_domain,IC,BC_0,BC_L);
u_data = u_true(end,:) + normrnd(0,meas_noise,1,Nx);

%% With observations and the true alpha, we now want to try and learn the random coefficients
func_call = @(alpha_est) heateq_1D(alpha_est,t_domain,x_domain,IC,BC_0,BC_L);
% u_obs = 
% Inverse problem: estimate diffusivity from observations
% Optimization setup
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
coeff_est = fminunc(@(coeff) obj_func_1Dheat_KL(coeff, u_data,alpha_mean, phi_true, lambda_true, func_call), coeff_true, options);

% Estimated diffusivity
alpha_est = alpha_mean + phi_true * (sqrt(lambda_true) .* coeff_est);

%% Plot the final temperature distribution
u_opt = heateq_1D(alpha_est,t_domain,x_domain,IC,BC_0,BC_L);
figure; 
subplot(1,2,1); hold on;
plot(x_domain, u_true(end,:), 'k','LineWidth',2);
plot(x_domain, u_data(end,:), 'ko','LineWidth',2);
plot(x_domain,u_opt(end,:),'-or','LineWidth',2)
xlabel('Position (m)');
ylabel('Temperature (C)');
set(gca,'FontSize',20); grid on;

subplot(1,2,2); hold on;
plot(x_domain, alpha_true, 'k','LineWidth',2);
plot(x_domain, alpha_est, '--r','LineWidth',2);
xlabel('Position (m)');
ylabel('Diffusivity Parameter');
set(gca,'FontSize',20); grid on;