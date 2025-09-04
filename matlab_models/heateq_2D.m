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
            da_dx = 0.5.*(alpha_xy(i+1,j) - alpha_xy(i-1,j))./dx;
            da_dy = 0.5.*(alpha_xy(i,j+1) - alpha_xy(i,j-1))./dy;
            du_dx = 0.5.*(u(i+1,j)-u(i-1,j))./dx;
            du_dy = 0.5.*(u(i,j+1)-u(i,j-1))./dy;
            d2u_dx = (u(i+1,j) - 2*u(i,j) + u(i-1,j))./(dx.^2);
            d2u_dy = (u(i,j+1) - 2*u(i,j) + u(i,j-1))./(dy.^2);
            u_new(i, j) = u(i, j) + ...
                dt.*(da_dx.*du_dx + alpha_xy(i,j).*d2u_dx + ...
                da_dy.*du_dy + alpha_xy(i,j).*d2u_dy);
        end
    end
    u = u_new;
    u_save(:,:,n) = u;
end


