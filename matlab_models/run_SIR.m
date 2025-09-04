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
gamma = 0.4; % infection coefficient, [0,1]
k     = 0.2; % contact rate for infection, [0,1]
r     = 0.4; % recovery rate, [0,1]
mu    = 0.2; % combined birth and death rate, [0,1]

% Put all the parameters together
param = [gamma; k; r; mu];


% Define the end time for numerical solution and the required step size
% (not necessary)
tend = 10;
tspace = linspace(0,tend,1000);

% Solve the system
[t,y] = ode45(@SIR_model,tspace,X0,[],param,N);


%% Plotting routine
figure(1); clf; hold on;
plot(t,y(:,1),':r','LineWidth',2);
plot(t,y(:,2),'--b','LineWidth',2);
plot(t,y(:,3),'-.k','LineWidth',2);
legend('Susceptible','Infected','Recovered');
grid on; set(gca,'FontSize',20);
ylabel('Number of Individuals')
xlabel('Time (days)')
