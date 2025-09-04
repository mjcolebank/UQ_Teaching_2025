% run_homogenous_spring.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the spring equation (Mass-Spring-Damper) with homogenous
% RHS, i.e. F(t) = 0

% Initial conditions
x0 = 2; v0 = -2;
IC = [x0;v0];

% The spring model has three parameters
m = 3;        % mass, kg
c = 2;        % damper friction, kg/s
k = 5;        % spring resistance, kg/s^2
f0 = 0.0;     % forcing term, kg cm / s^2
omegaF = 0.0; % frequency of forcing term, rad/s

% Put all the parameters together
param = [m;c;k;f0;omegaF];


% Define the end time for numerical solution and the required step size
% (latter not necessary)
tend = 20;
tspace = linspace(0,tend,1000);

% Solve the system
[t,y] = ode45(@spring_model,tspace,IC,[],param);


%% Plotting routine
figure(1); clf; hold on;
plot(t,y(:,1),':r','LineWidth',2);
plot(t,y(:,2),'--b','LineWidth',2);
legend('Displacement','Velocity');
grid on; set(gca,'FontSize',20);
xlabel('Time (s)')


