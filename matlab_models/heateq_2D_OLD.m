% heateq_2D.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Inputs:
% alpha: thermal diffusivity (scalar or vector)
% t: time (used in ODE solver)
% X: x-spatial dimension
% Y: y-spatial dimension
% IC: initial condition (should be a matrix for 2D)
% BCx_0, BCy_0: bounday condition (dirichlet) at x=0 or y=0
% BCx_L,BCy_L: bounday condition (dirichlet) at x=L or y=L
%
% Outputs:
% u_save: the output
function u_save = heateq_2D(alpha,t,X,Y,IC,...
                    BCx_0,BCy_0,BCx_L,BCy_L)
% Parameters
Nx = length(X);              % Number of spatial grid points (x-dir)
Ny = length(Y);              % Number of spatial grid points (y-dir)
Nt = length(t);            % Number of time steps
dx = X(2) - X(1);          % Spatial step size (x-dir)
dy = Y(2) - Y(1);          % Spatial step size (y-dir)  
dt = t(2) - t(1);          % Time step size


if isscalar(alpha)
    alpha_xy = alpha.*ones(Nx,Ny);
else
    alpha_xy = alpha;
end
rx = dt / dx^2;       % Stability parameter in the x direction
ry = dt / dy^2;       % Stability parameter in the y direction

% Initial condition
u = IC;

% Boundary conditions (Dirichlet)
u(1, :)   = BCx_0;       % Left boundary
u(end, :) = BCx_L;     % Right boundary
u(:, 1)   = BCy_0;       % Bottom boundary
u(:, end) = BCy_L;     % Top boundary
u_save = zeros(Nx,Nx,Nt);
% Time-stepping loop
for n = 1:Nt
    u_new = u;
    for i = 2:Nx-1
        for j = 2:Ny-1
            u_new(i, j) = u(i, j) + alpha_xy(i).*(rx * (u(i+1, j) - 2 * u(i, j) + u(i-1, j)) ...
                                   + ry * (u(i, j+1) - 2 * u(i, j) + u(i, j-1)));
        end
    end
    u = u_new;
    u_save(:,:,n) = u;
end


