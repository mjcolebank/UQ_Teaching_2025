% heateq_1D.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Inputs:
% alpha: thermal diffusivity (scalar or vector)
% t: time (used in ODE solver)
% X: spatial dimension
% IC: initial condition (should be vector)
% BC_0: bounday condition (dirichlet) at x=0
% BC_L: bounday condition (dirichlet) at x=L
%
% Outputs:
% u_save: the output

function u_save = heateq_1D_V2(alpha,t,X,IC,BC_0,BC_L)
% Parameters
Nx = length(X);              % Number of spatial grid points
Nt = length(t);            % Number of time steps
dx = X(2) - X(1);          % Spatial step size
dt = t(2) - t(1);          % Time step size



if isscalar(alpha)
    alpha_x = alpha.*ones(Nx,1);
else
    alpha_x = alpha;
end

% Initial condition
u = IC;

% Boundary conditions
u(1) = BC_0;    % Left boundary (Dirichlet condition)
u(end) = BC_L;  % Right boundary (Dirichlet condition)

u_save = zeros(Nt,Nx);
%% Discretize alpha

%%
% Time-stepping loop
for t_i = 1:Nt
    u_new = u;

    for i = 2:Nx-1
        d_alpha = 0.5.*(alpha_x(i+1) - alpha_x(i-1))./dx;
        du_dx = 0.5.*(u(i+1)-u(i-1))./dx;
        d2u_dx = (u(i+1) - 2*u(i) + u(i-1))./(dx.^2);
        u_new(i) = u(i) + dt.*(d_alpha.*du_dx + alpha_x(i).*d2u_dx);
    end
    u = u_new;
    u_save(t_i,:) = u;
end


end
