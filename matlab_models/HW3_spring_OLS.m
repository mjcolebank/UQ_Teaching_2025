% HW3_spring_OLS.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
% Uses the nonhomogenenous spring model and calibrates to data under the assumption of
% normal, gaussian iid errors
% Here we only consider displacement data
clear; clc; close all;
%%
% STUDENTS: define initial conditons
IC = [];

% STUDENTS: define the model parameters
C = [];        % damper friction, 1/s
K = [];        % spring resistance, 1/s^2
F0 = [];     % forcing term, cm / s^2
omegaF = []; % frequency of forcing term, rad/s

% STUDENTS: define the noise variance
noise_std = [];

% STUDENTS: define the true parameter (param_star) and an initial guess
% that is NOT the true value
param_star = []; % Get rid of mass term
param_guess = []; % We will fix f0 and omegaF

% STUDENTS: Chang the final time, number of sample points, etc.
n_data = []
tend = [];
tspace = linspace(0,tend,n_data);
n_xpts = length(tspace);

%%
num_param = 4;
true_data = call_spring_model(param_star,IC,tspace);
noisy_data = true_data+normrnd(0,noise_std,1,length(tspace));
initial_guess = call_spring_model(param_guess,IC,tspace);

% Solve the system to see how our apriori or nominal fits are
figure(1); clf; hold on;
plot(tspace,true_data,'-k','LineWidth',2);
plot(tspace,noisy_data,'ko','MarkerSize',8);
plot(tspace,initial_guess,'--r','LineWidth',2);


%% We want to find the parameters that best explain this data using an,
% ordinary least squares approach
% First, we need to define a cost functional, which will be a separate
% function call (defined at the bottom of this script

OLS_cost = @(q) get_spring_cost(q,tspace,noisy_data,IC);

% We will start by using the "f" function, which is UNCONSTRAINED
% optimization. Later we will try using "fmincon" which can handle bounded
% parameter domains


[param_opt,J_final] = fminsearch(OLS_cost,param_guess);
disp('Final estimate');
disp(param_opt');
disp('true value');
disp(param_star);


figure(1); hold on;
plot(tspace,call_spring_model(param_opt,IC,tspace),'--b','LineWidth',2);

%% After running the optmization, we can construct confidence intervals for the estimates
% and provide output prediction intervals
% First, we need the covariance, sigma^2 * inv(F^T*F)
% We approximate S via finite differences
S = zeros(n_xpts,num_param);
h = 1e-8;
for i=1:num_param
    param_step = param_opt;
    param_step(i) = param_opt(i)+h;
    S(:,i) = (call_spring_model(param_step,IC,tspace) - call_spring_model(param_opt,IC,tspace))./h;
end
%% Get estimates of noise variance, s2, and covariance, s2*(S^T*S)^-1
res = noisy_data' - call_spring_model(param_opt,IC,tspace)';
s2 = res'*res./(n_xpts-num_param);
covar = s2.*inv(S'*S);
% Confidence intervals in the parameters
CI_plus  = param_opt'+tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));
CI_minus = param_opt'-tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));

disp('Parameter CI')
disp([CI_minus CI_plus])

%% Now construct response confidence and prediction intervals using a similar formula
t_test = linspace(0,1.5.*max(tspace),100);
n_xtest = length(t_test);
% Need model sensitivity at test points
g = zeros(n_xtest,num_param);
h = 1e-8;
for i=1:num_param
    param_step = param_opt;
    param_step(i) = param_opt(i)+h;
    g(:,i) = (call_spring_model(param_step,IC,t_test) - call_spring_model(param_opt,IC,t_test))./h;
end

Y_pred = call_spring_model(param_opt,IC,t_test);
Y_CI = zeros(n_xtest,2);
Y_PI = zeros(n_xtest,2);
y_stderr_CI = zeros(n_xtest,1);
y_stderr_PI = zeros(n_xtest,1);

for i=1:n_xtest
    y_stderr_CI(i,1) = tinv(0.975,n_xpts-num_param).*sqrt(g(i,:)*covar*g(i,:)');
    Y_CI(i,1) = Y_pred(i) + tinv(0.975,n_xpts-num_param).*sqrt(g(i,:)*covar*g(i,:)');
     Y_CI(i,2) = Y_pred(i) - tinv(0.975,n_xpts-num_param).*sqrt(g(i,:)*covar*g(i,:)');

     y_stderr_PI(i,1) = tinv(0.975,n_xpts-num_param).*sqrt(s2+g(i,:)*covar*g(i,:)');
     Y_PI(i,1) = Y_pred(i) + tinv(0.975,n_xpts-num_param).*sqrt(s2+g(i,:)*covar*g(i,:)');
     Y_PI(i,2) = Y_pred(i) - tinv(0.975,n_xpts-num_param).*sqrt(s2+g(i,:)*covar*g(i,:)');
end

%% Plot the model response, CI, and PI
figure(2);clf;hold on;
plot(t_test,Y_pred,'-k','LineWidth',2);
plot(t_test,Y_CI,':b','LineWidth',2)
plot(t_test,Y_PI,':r','LineWidth',2);
plot(tspace,noisy_data,'ko','MarkerSize',8);
plot(t_test,call_spring_model(param_star,IC,t_test),'-c','LineWidth',2)
grid on;
set(gca,'FontSize',20);
axis tight;

% To show how the uncertainty is nonlinear because of a nonlinear model,
% plot the standard confidence and prediction error
figure(3);clf;hold on;
plot(t_test,y_stderr_CI,':b','LineWidth',2)
plot(t_test,y_stderr_PI,':r','LineWidth',2);
grid on;
set(gca,'FontSize',20);
axis tight;
%% Sub functions
function J = get_spring_cost(q,tspace,data,IC)
    % Take the parameter values passed in from the optmization routine and
    % generate simulations
   sim = call_spring_model(q,IC,tspace);


    % Now calculate the sum of squared errors
    residual = data - sim;
    J = sum(residual.^2);
end

function y = call_spring_model(q,IC,tspace)
[~,y] = ode45(@spring_model_4par,tspace,IC,[],q);
y = y(:,1)';
end

function dy = spring_model_4par(t,y,params)

% Unpack parameters
C      = params(1);
K      = params(2);
F0     = params(3);
omegaF = params(4);

% Redefine state variables for convenience
z = y(1);
v = y(2);


% RHS equations
accel_eq = (- C.*v - K.*z) + F0.*cos(omegaF.*t);
dy = [v;
      accel_eq];

end