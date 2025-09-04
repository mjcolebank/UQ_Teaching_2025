% SIR_model.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Inputs:
% t: time (used in ODE solver)
% y: current states of the system
% params: parameters of the system
% N: the total population considered in the model
%
% Outputs:
% dy: the vector containing RHS equations

function dy = SIR_model(t,y,params,N)

% Unpack parameters
gamma = params(1);
k     = params(2);
r     = params(3);
mu    = params(4);

% Redefine state variables for convenience
S = y(1);
I = y(2);
R = y(3);

% RHS equations
dy = [mu*N - mu*S - gamma*k*I*S;
    gamma*k*I*S - (r + mu)*I;
    r*I - mu*R];

end
