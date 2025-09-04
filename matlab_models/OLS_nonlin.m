% OLS_nonlin.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
% Defines a sine model and then calibrates to data under the assumption of
% normal, gaussian iid errors
xspace = linspace(-5,5,20 );

sin_model = @(q,x) q(1)+q(2).*sin(pi.*q(3).*x);

param_star = [5 1 0.3];
param_guess = [5 0.8 0.4];

num_param = 3;
n_xpts = length(xspace);

true_data = sin_model(param_star,xspace);
noisy_data = true_data+normrnd(0,0.5,1,20);
initial_guess = sin_model(param_guess,xspace);

% Solve the system to see how our apriori or nominal fits are
figure(1); clf; hold on;
plot(xspace,true_data,'-k','LineWidth',2);
plot(xspace,noisy_data,'ko','MarkerSize',8);
plot(xspace,initial_guess,'--r','LineWidth',2);

%% We want to find the parameters that best explain this data using an,
% ordinary least squares approach
% First, we need to define a cost functional, which will be a separate
% function call (defined at the bottom of this script

OLS_cost = @(q) get_sine_cost(q,xspace,noisy_data);

% We will start by using the "f" function, which is UNCONSTRAINED
% optimization. Later we will try using "fmincon" which can handle bounded
% parameter domains


[param_opt,J_final] = fminsearch(OLS_cost,param_guess);
disp('Final estimate');
disp(param_opt');
disp('true value');
disp(param_star);


figure(1); hold on;
plot(xspace,sin_model(param_opt,xspace),'--b','LineWidth',2);

%% After running the optmization, we can construct confidence intervals for the estimates
% and provide output prediction intervals
% First, we need the covariance, sigma^2 * inv(F^T*F)
% We approximate S via finite differences
S = zeros(n_xpts,num_param);
h = 1e-8;
for i=1:num_param
    param_step = param_opt;
    param_step(i) = param_opt(i)+h;
    S(:,i) = (sin_model(param_step,xspace) - sin_model(param_opt,xspace))./h;
end
%%
res = noisy_data' - sin_model(param_opt,xspace)';
s2 = res'*res./(n_xpts-num_param);
covar = s2.*inv(S'*S);
% Confidence intervals in the parameters
CI_plus  = param_opt'+tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));
CI_minus = param_opt'-tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));

disp('Parameter CI')
disp([CI_minus CI_plus])

%% Now construct response confidence and prediction intervals using a similar formula
x_test = linspace(1.5.*min(xspace),1.5.*max(xspace),100);
n_xtest = length(x_test);
% Need model sensitivity at test points
g = zeros(n_xtest,num_param);
h = 1e-8;
for i=1:num_param
    param_step = param_opt;
    param_step(i) = param_opt(i)+h;
    g(:,i) = (sin_model(param_step,x_test) - sin_model(param_opt,x_test))./h;
end

Y_pred = sin_model(param_opt,x_test);
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
plot(x_test,Y_pred,'-k','LineWidth',2);
plot(x_test,Y_CI,':b','LineWidth',2)
plot(x_test,Y_PI,':r','LineWidth',2);
plot(xspace,noisy_data,'ko','MarkerSize',8);
grid on;
set(gca,'FontSize',20);
axis tight;

% To show how the uncertainty is nonlinear because of a nonlinear model,
% plot the standard confidence and prediction error
figure(3);clf;hold on;
plot(x_test,y_stderr_CI,':b','LineWidth',2)
plot(x_test,y_stderr_PI,':r','LineWidth',2);
grid on;
set(gca,'FontSize',20);
axis tight;
%%
function J = get_sine_cost(q,xspace,data)
    % Take the parameter values passed in from the optmization routine and
    % generate simulations
   sim = q(1)+q(2).*sin(pi.*q(3).*xspace);


    % Now calculate the sum of squared errors
    residual = data - sim;
    J = sum(residual.^2);
end