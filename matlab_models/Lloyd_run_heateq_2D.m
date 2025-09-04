% run_heateq_2D.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the 1D heat equation. Note that, as of now, this script ONLY
% considers Dirichlet boundary conditions. You can use a constant
% conductivity, alpha, or describe a vector valued alpha.
close all; clear; clc;
%%
T  = 3; % End time for numerical solution
Lx = 1; % Length of the 1D element (x-dir)
Ly = 1; % Length of the 1D element (y-dir)
dx = 0.05; % Spatial discretization
dy = 0.05; % Spatial discretization
dt = 1e-3; % Temporal discretization
x_domain = 0:dx:Lx; % x points
y_domain = 0:dy:Ly; % y points
t_domain = 0:dt:T; % t points
Nx = length(x_domain); % Number of x points
Ny = length(y_domain); % Number of y points
Nt = length(t_domain); % Number of t points

% initial and boundary conditions (feel free to make more complicated)
IC = 0.25.*ones(Nx,Ny);% a matrix of all 1 entries
% u(round(Nx/2), round(Ny/2)) = 1; % Initial heat pulse in the center
% u(round(Nx/2), round(Ny/4)) = 1; 
% u(round(Nx/4), round(Ny/2)) = 1;
BCx_0 = 4+2*y_domain; % BC at x=0
BCy_0 = 5; % BC at y=0
BCx_L = 4+2*y_domain; % BC at x=L
BCy_L = 4+2*x_domain; % BC at y=L

%% Thermal conductivity parameter
% Note, we will consider this as a fixed constant and, later, as a random
% field. Also note that the numerical stability of the problem requires
% that alpha*dt/(dx^2)+alpha*dt/(dx^2)<0.5

% Constant, uniform thermal conductivity
alpha = 0.01;         % Thermal diffusivity constant (m^2/s)

% Non-uniform conductivity
%alpha = 0.01.*cos(x_domain./10)+0.06.*sin(y_domain'./10);

stab_number = max(alpha).*dt./(dx.^2)+max(alpha).*dt./(dy.^2);
if stab_number >= 0.5
    warning('Numerical Instability may occur.')
end

u_save = heateq_2D(alpha,t_domain,x_domain,y_domain,IC,...
                    BCx_0,BCy_0,BCx_L,BCy_L);


%%
% Plot the final temperature distribution
[X, Y] = meshgrid(x_domain, y_domain);

figure;
subplot(1,2,1);
surf(X, Y, u_save(:,:,1));
xlabel('x (m)');
ylabel('y (m)');
zlabel('Temperature (°C)');
title('Initial temperature distribution');
colorbar;
clim([0 5]); zlim([0 5])

subplot(1,2,2);
surf(X, Y, u_save(:,:,end));
xlabel('x (m)');
ylabel('y (m)');
zlabel('Temperature (°C)');
title('Final temperature distribution');
colorbar;
clim([0 5]); zlim([0 5])

%%
figure;
subplot(1,3,1);
heatmap(u_save(:,:,1)); %colormap('gray'); 
xlabel('x (m)');
ylabel('y (m)');
set(gca,'FontSize',20); grid on; clim([0 5]);
title('Temperature at t=0');

subplot(1,3,2);
heatmap(u_save(:,:,floor(end/2))); %colormap('gray');  
xlabel('x (m)');
set(gca,'FontSize',20); grid on;  clim([0 5]);
title('Temperature at t=T/2');

subplot(1,3,3);
heatmap(u_save(:,:,end)); %colormap('gray');  
xlabel('x (m)');
set(gca,'FontSize',20); grid on;  clim([0 5]);
title('Temperature at t=T');
colorbar;

%% Complex Step Method

h = 1e-6; % Step Size

% THERMAL DIFFUSIVITY (alpha)

sens_alpha = zeros(Nx, Ny); % Preallocate for sensitivity
alpha_perturbed = alpha + 1i*h; % Perturb by imaginary step
u_perturbed_alpha = heateq_2D(alpha_perturbed, t_domain, x_domain, y_domain, IC, BCx_0, BCy_0, BCx_L, BCy_L); % Plug in pertubation and evaluate model
final_temp_perturbed_alpha = u_perturbed_alpha(:,:,end); % Evaluate final temperature 
sens_alpha = imag(final_temp_perturbed_alpha) / h; % Take imaginary part as sensitivity

% Plot sensitivity
figure; 
surf(x_domain, y_domain, sens_alpha');
xlabel('x');
ylabel('y');
zlabel('Sensitivity of Temperature to Alpha');
title('Complex Step Sensitivity: Final Time');
colorbar;


% BOUNDARIES
% LEFT BOUNDARY (X=0)

sens_BCx_0 = zeros(Nx, Ny); % Preallocate for sensitivity
BCx_0_perturbed = BCx_0 + 1i*h; % Perturb by imaginary step
u_perturbed_BCx_0 = heateq_2D(alpha, t_domain, x_domain, y_domain, IC, BCx_0_perturbed, BCy_0, BCx_L, BCy_L); % Plug in pertubation and evaluate model
final_temp_perturbed_BCx_0 = u_perturbed_BCx_0(:,:,end); % Evaluate final temperature 
sens_BCx_0 = imag(final_temp_perturbed_BCx_0) / h; % Take imaginary part as sensitivity

% LOWER BOUNDARY (Y=0)

sens_BCy_0 = zeros(Nx, Ny); % Preallocate for sensitivity
BCy_0_perturbed = BCy_0 + 1i*h; % Perturb by imaginary step
u_perturbed_BCy_0 = heateq_2D(alpha, t_domain, x_domain, y_domain, IC, BCx_0, BCy_0_perturbed, BCx_L, BCy_L); % Plug in pertubation and evaluate model
final_temp_perturbed_BCy_0 = u_perturbed_BCy_0(:,:,end); % Evaluate final temperature 
sens_BCy_0 = imag(final_temp_perturbed_BCy_0) / h; % Take imaginary part as sensitivity

% RIGHT BOUNDARY (X=L)

sens_BCx_L = zeros(Nx, Ny); % Preallocate for sensitivity
BCx_L_perturbed = BCx_L + 1i*h; % Perturb by imaginary step
u_perturbed_BCx_L = heateq_2D(alpha, t_domain, x_domain, y_domain, IC, BCx_0, BCy_0, BCx_L_perturbed, BCy_L); % Plug in pertubation and evaluate model
final_temp_perturbed_BCx_L = u_perturbed_BCx_L(:,:,end); % Evaluate final temperature 
sens_BCx_L = imag(final_temp_perturbed_BCx_L) / h; % Take imaginary part as sensitivity

% UPPER BOUNDARY (Y=L)

sens_BCy_L = zeros(Nx, Ny); % Preallocate for sensitivity
BCy_L_perturbed = BCy_L + 1i*h; % Perturb by imaginary step
u_perturbed_BCy_L = heateq_2D(alpha, t_domain, x_domain, y_domain, IC, BCx_0, BCy_0, BCx_L, BCy_L_perturbed); % Plug in pertubation and evaluate model
final_temp_perturbed_BCy_L = u_perturbed_BCy_L(:,:,end); % Evaluate final temperature 
sens_BCy_L = imag(final_temp_perturbed_BCy_L) / h; % Take imaginary part as sensitivity


% Plot boundary sensitivities
figure;
subplot (2,2,1);
surf(x_domain, y_domain, sens_BCx_0');
xlabel('x');
ylabel('y');
zlabel('Sensitivity');
title('Left Boundary');
colorbar;

subplot (2,2,2);
surf(x_domain, y_domain, sens_BCy_0');
xlabel('x');
ylabel('y');
zlabel('Sensitivity');
title('Lower Boundary');
colorbar;

subplot(2,2,3);
surf(x_domain, y_domain, sens_BCx_L');
xlabel('x');
ylabel('y');
zlabel('Sensitivity');
title('Right Boundary');
colorbar;

subplot(2,2,4);
surf(x_domain, y_domain, sens_BCy_L');
xlabel('x');
ylabel('y');
zlabel('Sensitivity');
title('Upper Boundary');
colorbar;

sgtitle({'\bf Sensitivity of the 2-D Heat Equation to Boundary Condition Changes'}, 'FontSize', 10);

%% MCMC to infer constant alpha


% Noise Variance under different error models
% Constant Noise Variance:
%noise_var = 0.01;

% Spatially-dependent noise
noise_var =  0.01 + 0.05 * X; 


%define true alpha as param_star
param_star = 0.01;
theta0 = 0.012;             % Initial guess
M = 10000;                   % Number samples
s_theta = 1e-4;
covar = s_theta^2;          % proposal variance
UB = 0.02;
LB = 0.005; 

u_true = heateq_2D(param_star,t_domain,x_domain,y_domain,IC,...
                    BCx_0,BCy_0,BCx_L,BCy_L);

u_final = u_true(:,:,end);

data = u_final + sqrt(noise_var) * randn(size(u_final));

f_mod = @(alpha) get_final_temp(heateq_2D(alpha, t_domain, x_domain, y_domain, IC, ...
                                   BCx_0, BCy_0, BCx_L, BCy_L));

% do i need to do an if statement to make sure that our parameter is
% between these values?
prior_F = @(theta) 1 / (UB - LB);

% MJC: before you pass the data in, lets make it a single vector
data = data(:);
[chain, s2] = Metropolis_errupdate(f_mod, data, prior_F, theta0, M, covar, UB, LB);



%%
figure(1);clf

hold on;
plot(chain(1));
yline(param_star,'--k','LineWidth',2);
title('Alpha');

figure(2);clf

hold on;
ksdensity(chain(1));
xline(param_star,'--k','LineWidth',2);
title('Alpha');


%%
% Plot trace and histogram
figure;
subplot(2,1,1); 
plot(chain); 
title('MCMC Path of \\alpha'); 
ylabel('\\alpha');

subplot(2,1,2); 
histogram(chain, 40); 
title('Posterior of \\alpha'); 
xlabel('\\alpha');

figure;
subplot(2,1,1); 
plot(s2); 
title('MCMC Path of \\sigma^2'); 
ylabel('\\sigma^2');

subplot(2,1,2); 
histogram(s2, 40); 
title('Posterior of \\sigma^2'); 
xlabel('\\sigma^2');

%% Temperature plots of posterior vs true

num_samples = 100;
u_post = zeros(Nx, Ny, num_samples);
for i = 1:num_samples
    alpha_i = chain(1,end - i + 1);
    u_i = heateq_2D(alpha_i, t_domain, x_domain, y_domain, IC, ...
                    BCx_0, BCy_0, BCx_L, BCy_L);
    u_post(:,:,i) = u_i(:,:,end);
end

posterior_mean = mean(u_post, 3);
posterior_std  = std(u_post, 0, 3);

% Plot posterior predictive mean
figure;
surf(x_domain, y_domain, posterior_mean'); 
shading interp
xlabel('x'); 
ylabel('y'); 
zlabel('Temperature');
title('Posterior Predictive Mean'); 
colorbar;

% Plot posterior predictive uncertainty
figure;
surf(x_domain, y_domain, posterior_std'); 
shading interp
xlabel('x'); 
ylabel('y'); 
zlabel('Std Dev');
title('Posterior Predictive Std Dev'); 
colorbar;

% Compare true and inferred final temperature
figure;
subplot(1,3,1);
surf(x_domain, y_domain, u_true(:,:,end)');
title('True Final Temperature'); 
shading interp;

subplot(1,3,2);
surf(x_domain, y_domain, posterior_mean');
title('Posterior Predictive Mean'); 
shading interp;

subplot(1,3,3);
surf(x_domain, y_domain, posterior_std');
title('Predictive Std Dev'); 
shading interp;


function u_end = get_final_temp(u_full)
    u_end = u_full(:,:,end);
    % MJC: we need to provide a vector instead of matrix
    u_end = u_end(:);
end