% run_heateq_1D.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the 1D heat equation. Note that, as of now, this script ONLY
% considers Dirichlet boundary conditions. You can use a constant
% conductivity, alpha, or describe a vector valued alpha.

T = 5; % End time for numerical solution
L = 1; % Length of the 1D element
dx = 0.05; % Spatial discretization
dt = 1e-3; % Temporal discretization
x_domain = 0:dx:L; % x points
t_domain = 0:dt:T; % t points
Nx = length(x_domain); % Number of x points
Nt = length(t_domain); % Number of t points

% initial and boundary conditions (feel free to make more complicated)
IC = 1.*ones(Nx,1);
BC_0 = 4; % BC at x=0
BC_L = 5; % BC at x=L

%% Thermal conductivity parameter
% Note, we will consider this as a fixed constant and, later, as a random
% field. Also note that the numerical stability of the problem requires
% that alpha*dt/(dx^2)

% Constant, uniform thermal conductivity
alpha = 0.01;         % Thermal diffusivity constant (m^2/s)

% Non-uniform conductivity
% alpha = linspace(0.05,0.002,Nx);

stab_number = max(alpha).*dt./(dx.^2);
if stab_number >= 0.5
    warning('Numerical Instability may occur.')
end
u_save = heateq_1D(alpha,t_domain,x_domain,IC,BC_0,BC_L);


% Plot the final temperature distribution
figure; hold on;
plot(x_domain, u_save(1,:), '-o');
plot(x_domain, u_save(round(Nt/2),:), '-o');
plot(x_domain, u_save(end,:), '-o');
legend('t=0','t=T/2','t=T')
xlabel('Position (m)');
ylabel('Temperature (C)');
title('Temperature distribution along the rod');
set(gca,'FontSize',20); grid on;



