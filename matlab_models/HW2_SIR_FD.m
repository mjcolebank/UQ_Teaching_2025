% run_SIR.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the SIR model

% Initial conditions
S0 = 999; I0 = 1; R0  = 0;
X0 = [S0;I0;R0];

% Define the total population based on the initial conditions
N = sum(X0);

% The SIR model has four parameters, as described below
gamma = 0.1; % infection coefficient, [0,1]
k     = 0.2; % contact rate for infection, [0,1]
r     = 0.4; % recovery rate, [0,1]
mu    = 0.2; % combined birth and death rate, [0,1]

% gamma = 0.4; % infection coefficient, [0,1]
% k     = 0.2; % contact rate for infection, [0,1]
% r     = 0.4; % recovery rate, [0,1]
% mu    = 0.2; % combined birth and death rate, [0,1]

% Put all the parameters together
param_star = [gamma; k; r; mu];
num_param = 4;

% Define the end time for numerical solution and the required step size
% (not necessary)
n_xpts = 200;
tend = 10;
tspace = linspace(0,tend,n_xpts);

% Solve the system
[t,y] = ode45(@SIR_model,tspace,X0,[],param_star,N);
y_base = y(:,2);
%% Loop over the parameter values and perturb one element at a time
f_sens = zeros(num_param,n_xpts); % Declare storage
h = 1e-6; % STEP SIZE
for i=1:num_param
    param_step = param_star;
    param_step(i) = param_star(i)+h;
    [t,y] = ode45(@SIR_model,tspace,X0,[],param_step,N);
    f_sens(i,:) = (y(:,2) - y_base)./h;
end

%% loop over initial conditions
f_IC = zeros(3,n_xpts); % Declare storage
h = 1e-6; % STEP SIZE
for i=1:3
    x0_step = X0;
    x0_step(i) = x0_step(i)+h;
    [t,y] = ode45(@SIR_model,tspace,x0_step,[],param_star,N);
    f_IC(i,:) = (y(:,2) - y_base)./h;
end


%% Plotting routine
figure;
plot(tspace,f_sens);
legend({'$S_\gamma$','$S_k$','$S_r$','$S_{\mu}$'},'Interpreter','latex')

figure;
plot(tspace,f_IC);
legend({'$S_{S_0}$','$S_{I_0}$','$S_{R_0}$'},'Interpreter','latex')