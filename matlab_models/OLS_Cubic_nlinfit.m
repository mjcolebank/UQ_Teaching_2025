% OLS_Cubic.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines a cubic model and then calibrates to data under the assumption of
% normal, gaussian iid errors

xspace = linspace(-5,5, 20);
% xspace = linspace(-0.5,0.5, 20);

cubic_model = @(q,x) q(1) + q(2).*x ...
            + q(3).*x.^2 + q(4).*x.^3;

param_star = [5 0.1 0.05 0.01];
param_guess = [4 0.2 0.07 0.03];
% param_star = [0.1 2 1e-5 1e-5];

num_param = 4;
n_xpts = length(xspace);

true_data = cubic_model(param_star,xspace);
noisy_data = true_data+normrnd(0,0.25,1,20);
initial_guess = cubic_model(param_guess,xspace);

% Solve the system to see how our apriori or nominal fits are
figure(1); clf; hold on;
plot(xspace,true_data,'-k','LineWidth',2);
plot(xspace,noisy_data,'ko','MarkerSize',8);
plot(xspace,initial_guess,'--r','LineWidth',2);

%% We want to find the parameters that best explain this data using an,
% ordinary least squares approach


% We will start by using the "nlinfit" function, which is UNCONSTRAINED
% optimization. Later we will try using "fmincon" which can handle bounded
% parameter domains

% Note: the nlinfit function is nice because it calculates the optimal
% parameter value and returns the residual, jacobian, and covariance
% estimates
[param_opt,residual_final,Jacob_opt,COVB,MSE] = nlinfit(xspace,noisy_data,cubic_model,param_guess);
disp('Final estimate');
disp(param_opt');
disp('true value');
disp(param_star);

figure(1); hold on;
plot(xspace,cubic_model(param_opt,xspace),'--b','LineWidth',2);

%% After running the optmization, we can construct confidence intervals for the estimates
% and provide output prediction intervals
% The "nlpredci" function takes the model and independent variable and
% provides predicted values at the independent variable and the associated
% uncertainty in terms of the confidence interval
xtest = linspace(-7,7,30);
[YPRED, delta] = nlpredci(cubic_model,xtest,param_opt,residual_final,'Covariance',COVB);

figure(2);clf; hold on;
plot(xtest,YPRED,'--k','LineWidth',2);
plot(xtest,YPRED+delta,'--r','LineWidth',2);
plot(xtest,YPRED-delta,'--r','LineWidth',2);
plot(xspace,noisy_data,'ko','MarkerSize',8);

% We can also get parameter confidence intervals
CI = nlparci(param_opt,residual_final,'covariance',COVB);
disp('parameter lower/upper confidence intervals')
disp(CI)

