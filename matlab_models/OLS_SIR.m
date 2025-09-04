% OLS_SIR.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the SIR model and the calibrate to data under the assumption of
% normal, gaussian iid errors

% Initial conditions
S0 = 999; I0 = 1; R0  = 0;
X0 = [S0;I0;R0];

% Define the total population based on the initial conditions
N = sum(X0);

% Initial guess on the parameter values
gamma = 0.1; % infection coefficient, [0,1]
k     = 0.2; % contact rate for infection, [0,1]
r     = 0.2; % recovery rate, [0,1]
mu    = 0.2; % combined birth and death rate, [0,1]

% Put all the parameters together
param = [gamma; k; r; mu];


%% Now, we can load in the data that we would like to calibrate to.
% The file "SIR_data" contains t_data and I_data (as well as the true
% parameter values for checking if our optimization was correct).
% I_data contains 34 equispaced points of noisy infected data
% t_data is the corresponding time points
load SIR_data.mat I_data t_data
% We specify the t values as the data values
tspace = t_data;

% Solve the system to see how our apriori or nominal fits are
[t,y] = ode45(@SIR_model,tspace,X0,[],param,N);
figure(1); clf; hold on;
plot(t,y(:,2),'--r','LineWidth',2);
plot(t_data,I_data,'ko','MarkerSize',8);
%% We want to find the parameters that best explain this data using an,
% ordinary least squares approach
% First, we need to define a cost functional, which will be a separate
% function call (defined at the bottom of this script

OLS_cost = @(q) get_SIR_cost(q,tspace,X0,N,I_data);

% We will start by using the "fminsearch" function, which is UNCONSTRAINED
% optimization. Later we will try using "fmincon" which can handle bounded
% parameter domains


[param_opt,J_final] = fminsearch(OLS_cost,param);

%% After running the optmization, we can construct confidence intervals for the estimates
% and provide output prediction intervals
% First, we need the covariance

[YPRED, DELTA] = nlpredci(MODELFUN,X,BETA,RESID,'covar',SIGMA)
CI = nlparci(BETA,RESID,'covar',SIGMA)

function J = get_SIR_cost(q,tspace,X0,N,I_data)
    % Take the parameter values passed in from the optmization routine and
    % generate SIR predictions
    [~,y_sim] = ode45(@SIR_model,tspace,X0,[],q,N);

    % Now calculate the sum of squared errors
    residual = I_data - y_sim(:,2);
    J = sum(residual.^2);
end