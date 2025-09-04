% spring_model.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Inputs:
% t: time (used in ODE solver)
% y: current states of the system
% params: parameters of the system
%
% Outputs:
% dy: the vector containing RHS equations

function dy = spring_model(t,y,params)

% Unpack parameters
m      = params(1);
c      = params(2);
k      = params(3);
f0     = params(4);
omegaF = params(5);

% Redefine state variables for convenience
x = y(1);
v = y(2);


% RHS equations
accel_eq = (- c.*v - k.*x)./m + f0.*cos(omegaF.*t)/m;
dy = [v;
      accel_eq];

end