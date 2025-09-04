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
T  = 5; % End time for numerical solution
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
IC = 1.*ones(Nx,Ny);
% u(round(Nx/2), round(Ny/2)) = 1; % Initial heat pulse in the center
% u(round(Nx/2), round(Ny/4)) = 1; 
% u(round(Nx/4), round(Ny/2)) = 1;
BCx_0 = 4; % BC at x=0
BCy_0 = 4; % BC at y=0
BCx_L = 5; % BC at x=L
BCy_L = 5; % BC at y=L

%% Thermal conductivity parameter
% Note, we will consider this as a fixed constant and, later, as a random
% field. Also note that the numerical stability of the problem requires
% that alpha*dt/(dx^2)+alpha*dt/(dx^2)<0.5

% Constant, uniform thermal conductivity
alpha = 0.01;         % Thermal diffusivity constant (m^2/s)

% Non-uniform conductivity
% alpha = 0.01.*cos(x_domain./10)+0.06.*sin(y_domain'./10);

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

