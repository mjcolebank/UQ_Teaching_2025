% OLS_Cubic.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines a cubic model and then calibrates to data under the assumption of
% normal, gaussian iid errors

% Region where parameters are identifiable.
xspace = linspace(-5,5,20 );
% Region where parameters are NOT identifiable.
% xspace = linspace(-0.5,0.5,20 );

cubic_model = @(q,x) q(1) + q(2).*x ...
            + q(3).*x.^2 + q(4).*x.^3;

param_star = [5 0.1 0.05 0.01];
param_guess = [4 0.2 0.07 0.03];

num_param = 4;
n_xpts = length(xspace);

true_data = cubic_model(param_star,xspace);
noisy_data = true_data+normrnd(0,0.3,1,20);
initial_guess = cubic_model(param_guess,xspace);

% Solve the system to see how our apriori or nominal fits are
figure(1); clf; hold on;
plot(xspace,true_data,'-k','LineWidth',2);
plot(xspace,noisy_data,'ko','MarkerSize',8);
plot(xspace,initial_guess,'--r','LineWidth',2);

%% We want to find the parameters that best explain this data using an,
% ordinary least squares approach
% First, we need to define a cost functional, which will be a separate
% function call (defined at the bottom of this script

OLS_cost = @(q) get_cubic_cost(q,xspace,noisy_data);

% We will start by using the "f" function, which is UNCONSTRAINED
% optimization. Later we will try using "fmincon" which can handle bounded
% parameter domains


[param_opt,J_final] = fminsearch(OLS_cost,param_guess);
disp('Final estimate');
disp(param_opt);
disp('true value');
disp(param_star);


figure(1); hold on;
plot(xspace,cubic_model(param_opt,xspace),'--b','LineWidth',2);

%% After running the optmization, we can construct confidence intervals for the estimates
% and provide output prediction intervals
% First, we need the covariance, sigma^2 * inv(F^T*F)
X_des = [ones(1,n_xpts); xspace; xspace.^2; xspace.^3]';
res = noisy_data' - cubic_model(param_opt,xspace)';
s2 = res'*res./(n_xpts-num_param);
covar = s2.*inv(X_des'*X_des);
% Confidence intervals in the parameters
CI_plus  = param_opt'+tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));
CI_minus = param_opt'-tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));

disp('Parameter CI')
disp([CI_minus CI_plus])

%% Now construct response confidence and prediction intervals using a similar formula
x_test = linspace(1.5.*min(xspace),1.5.*max(xspace),30); % Since there are negatives here, we use 50% on either side of the domain
n_xtest = length(x_test);
X_test_des = [ones(1,n_xtest); x_test; x_test.^2; x_test.^3]';

Y_pred = cubic_model(param_opt,x_test);
% Y_CI_plus = Y_pred + tinv(0.975,n_xpts-num_param).*sqrt(X_test_des*diag(covar)*X_test_des');
Y_CI = zeros(n_xtest,2);
Y_PI = zeros(n_xtest,2);

for i=1:n_xtest
    Y_CI(i,1) = Y_pred(i) + tinv(0.975,n_xpts-num_param).*sqrt(X_test_des(i,:)*covar*X_test_des(i,:)');
     Y_CI(i,2) = Y_pred(i) - tinv(0.975,n_xpts-num_param).*sqrt(X_test_des(i,:)*covar*X_test_des(i,:)');
     Y_PI(i,1) = Y_pred(i) + tinv(0.975,n_xpts-num_param).*sqrt(s2+X_test_des(i,:)*covar*X_test_des(i,:)');
     Y_PI(i,2) = Y_pred(i) - tinv(0.975,n_xpts-num_param).*sqrt(s2+X_test_des(i,:)*covar*X_test_des(i,:)');
end

figure(2);clf;hold on;
plot(x_test,Y_pred,'-k','LineWidth',2);
plot(x_test,Y_CI,':b','LineWidth',2)
plot(x_test,Y_PI,':r','LineWidth',2);
plot(xspace,noisy_data,'ko','MarkerSize',8);
grid on;
set(gca,'FontSize',20);
axis tight;
%%
function J = get_cubic_cost(q,xspace,data)
    % Take the parameter values passed in from the optmization routine and
    % generate simulations
   sim = q(1) + q(2).*xspace ...
            + q(3).*xspace.^2 + q(4).*xspace.^3;


    % Now calculate the sum of squared errors
    residual = data - sim;
    J = sum(residual.^2);
end