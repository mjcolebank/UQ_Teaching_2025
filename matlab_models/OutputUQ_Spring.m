% OutputUQ_Spring
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
% Uses the homogenenous spring model and calibrates to data under the assumption of
% normal, gaussian iid errors
% Here we only consider displacement data
clear; clc; close all;
%%
% Initial conditions
x0 = 2; v0 = -2;
IC = [x0;v0];

% Solve the unforced spring model with 2 parameters
m = 12; %3, 5, 10        % mass, kg
c = 2;        % damper friction, kg/s
k = 5;        % spring resistance, kg/s^2

% Noise paramters; here we will examine things WITHOUT optimizing
noise_std = 0.2;
noise_var = noise_std.^2;


% Put all the parameters together
% Only look at C=c/m and K=k/m
param_star = [c/m;k/m];

tend = 20;
tspace = linspace(0,tend,30);
t_test = linspace(0,1.5.*max(tspace),100);

num_param = 2;
n_xpts = length(tspace);


%% Construct confidence intervals for the estimates
% and provide output prediction intervals
% First, we need the covariance, sigma^2 * inv(F^T*F)
% We approximate S via finite differences
S = zeros(n_xpts,num_param);
h = 1e-8;
for i=1:num_param
    param_step = param_star;
    param_step(i) = param_star(i)+h;
    S(:,i) = (call_spring_model(param_step,IC,tspace) - call_spring_model(param_star,IC,tspace))./h;
end
%%
s2 = noise_var;
covar = s2.*inv(S'*S);
% Confidence intervals in the parameters
CI_plus  = param_star+tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));
CI_minus = param_star-tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));

disp('Parameter CI')
disp([CI_minus param_star CI_plus])

%% Now construct response confidence and prediction intervals using a similar formula
n_xtest = length(t_test);
% Need model sensitivity at test points
g = zeros(n_xtest,num_param);
h = 1e-8;
for i=1:num_param
    param_step = param_star;
    param_step(i) = param_star(i)+h;
    g(:,i) = (call_spring_model(param_step,IC,t_test) - call_spring_model(param_star,IC,t_test))./h;
end

Y_pred = call_spring_model(param_star,IC,t_test);
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
figure(1);clf;
subplot(1,2,1); hold on;
plot(t_test,Y_pred,'-k','LineWidth',2);
plot(t_test,Y_CI,':b','LineWidth',2)
plot(t_test,Y_PI,':r','LineWidth',2);
grid on;
set(gca,'FontSize',20);
axis tight;

% To show how the uncertainty is nonlinear because of a nonlinear model,
% plot the standard confidence and prediction error
subplot(1,2,2); hold on;
plot(t_test,y_stderr_CI,':b','LineWidth',2)
plot(t_test,y_stderr_PI,':r','LineWidth',2);
grid on;
set(gca,'FontSize',20);
axis tight;

%% Now, consider a monte carlo estimate
ub = 1.2.*param_star;
lb = 0.8.*param_star;
M = 10000;
mu = zeros(1,n_xtest);
mu2 = zeros(1,n_xtest);
C_MC = normrnd(param_star(1),0.1.*param_star(1),M,1);%unifrnd(lb(1),ub(1),M,1);
K_MC = normrnd(param_star(2),0.1.*param_star(2),M,1);%unifrnd(lb(2),ub(2),M,1);

for i=1:M
    param_MC = [C_MC(i) K_MC(i)];
    y = call_spring_model(param_MC,IC,t_test);
    mu = mu+y;
    mu2 = mu2+y.^2;
end
mu = mu./M;
mu2 = mu2./M;
var = mu2 - mu.^2;
%%
figure(20);clf;
subplot(1,2,1); hold on;
plot(t_test,mu,'-k','LineWidth',2);
plot(t_test,mu+2.*(sqrt(var)),':r','LineWidth',2)
plot(t_test,mu-2.*(sqrt(var)),':r','LineWidth',2);
grid on;
set(gca,'FontSize',20);
axis tight;

% To show how the uncertainty is nonlinear because of a nonlinear model,
% plot the standard confidence and prediction error
subplot(1,2,2); hold on;
plot(t_test,sqrt(var),':b','LineWidth',2)
grid on;
set(gca,'FontSize',20);
axis tight;

%%
figure(40);clf; hold on;
plot(t_test,Y_pred,':r','LineWidth',2)
plot(t_test,Y_CI,':b','LineWidth',2)
plot(t_test,mu,'-.k','LineWidth',2);
plot(t_test,mu+2.*(sqrt(var)),':c','LineWidth',2)
plot(t_test,mu-2.*(sqrt(var)),':c','LineWidth',2);
grid on;
set(gca,'FontSize',20);
axis tight;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = call_spring_model(q,IC,tspace)
[~,y] = ode45(@spring_model_2par,tspace,IC,[],q);
y = y(:,1)';
end

function dy = spring_model_2par(t,y,params)

% Unpack parameters
C      = params(1);
K      = params(2);


% Redefine state variables for convenience
x = y(1);
v = y(2);


% RHS equations
accel_eq = (- C.*v - K.*x);
dy = [v;
      accel_eq];

end