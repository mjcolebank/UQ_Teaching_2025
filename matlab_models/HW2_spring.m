% run_nonhomogenous_spring.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the spring equation (Mass-Spring-Damper) with nonhomogenous
% RHS, i.e. F(t) = f0*cos(omega*t)

% Initial conditions
x0 = 2; v0 = -2;
IC = [x0;v0];

% The spring model has three parameters
m = 3;        % mass, kg
c = 2;        % damper friction, kg/s
k = 5;        % spring resistance, kg/s^2
f0 = 3.0;     % forcing term, kg cm / s^2
omegaF = 2.0; % frequency of forcing term, rad/s

% Put all the parameters together
param_star = [m;c;k;f0;omegaF];
num_param = 5;

% Define the end time for numerical solution and the required step size
% (latter not necessary)
tend = 20;
x_pts = 200;
tspace = linspace(0,tend,x_pts);

% Solve the system
[t,y] = ode45(@spring_model,tspace,IC,[],param_star);
y_base = y(:,1);

%% Loop over the parameter values and perturb one element at a time
f_sens = zeros(num_param,n_xpts); % Declare storage
h = 1e-6; % STEP SIZE
for i=1:num_param-1
    param_step = param_star;
    param_step(i) = param_star(i)+h;
    [t,y] = ode45(@spring_model,tspace,IC,[],param_step);
    f_sens(i,:) = (y(:,1) - y_base)./h;
end

%% loop over initial conditions
f_IC = zeros(2,n_xpts); % Declare storage
h = 1e-6; % STEP SIZE
for i=1:2
    x0_step = IC;
    x0_step(i) = param_star(i)+h;
    [t,y] = ode45(@spring_model,tspace,x0_step,[],param_star);
    f_IC(i,:) = (y(:,2) - y_base)./h;
end


%% Plotting routine
figure;
plot(tspace,f_sens);
legend({'$S_m$','$S_c$','$S_k$','$S_{f0}$','$S_{\omega}$'},'Interpreter','latex')

figure;
plot(tspace,f_IC);
legend({'$S_{S_0}$','$S_{I_0}$','$S_{R_0}$'},'Interpreter','latex')

