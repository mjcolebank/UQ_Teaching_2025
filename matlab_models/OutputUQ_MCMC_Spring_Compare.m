% Metropolis_Regression.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
clear; clc; close all;
%% Compute the posterior distribution under the assumption of a uniform
% Initial conditions
z0 = 2; v0 = -2;
IC = [z0;v0];

% The spring model has three parameters
C = 0.66;        % damper friction, kg/s
K = 2;        % spring resistance, kg/s^2
F0 = 1;     % forcing term, kg cm / s^2
omegaF = 0.8; % frequency of forcing term, rad/s

% Put all the parameters together
param_star = [C;K;F0;omegaF]'; % Get rid of mass term

tend = 20;
tspace = linspace(0,tend,30);
num_param = 4;
n_xpts = length(tspace);


f_mod = @(q) call_spring_model(q,IC,tspace);
true_signal = f_mod(param_star);
%% Things that can be adjusted
% Measurement noise
noise_var = 0.01;

% MCMC algorithm parameters
M = 10000;
theta0 = param_star.*1.00;
data        = true_signal + normrnd(0,sqrt(noise_var),1,n_xpts);

s_theta = 1e-1.*ones(num_param,1); %1e-3, 0.1
covar_const = eye(num_param).*s_theta;

%% In this example, we provide an "estimated" covariance using sensitivity analysis
S = zeros(n_xpts,num_param);
h = 1e-6;
f_theta0 = f_mod(theta0);
for i=1:num_param
    param_step = theta0;
    param_step(i) = param_step(i)+h;
    S(:,i) = (f_mod(param_step) - f_theta0)./h;
end
F = S'*S;
s2_est = sum( (data-f_mod(theta0)).^2)./(n_xpts - num_param);

covar_est = s2_est.*inv(F);
%% First, construct the prior.
% We start by assuming uniform (consider Gaussian) on a +-20% range
UB_uni = [4; 4; 4; 4]; % Twenty percent above
LB_uni = [0.1; 0.1; 0.1; 0.1]; % Twenty percent above

prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes
prior_F = @(param) prior_unif(param,LB_uni,UB_uni);

%% Chain when covariance is from sensitivity
k0 = 1000;
[chain,s2] = AdaptiveMetropolis_s2update(f_mod,data,prior_F,theta0,k0,M,covar_est,UB_uni,LB_uni);

%% Get MAP points
MAP_const = zeros(num_param,1);
for i=1:num_param
    % Get density values
    [f,x] = ksdensity(chain(i,:));
    % Find the maximum of the PDF
    [who,where] = max(f);
    % Store it
    MAP_const(i) = x(where);
end

%% Construct Bayesian credible + Pred Intervals
M_CIPI = 2000;
y_CI = zeros(n_xpts,M_CIPI);
y_PI = zeros(n_xpts,M_CIPI);

for i=1:M_CIPI % We will draw the last 2000 samples from the chain
    par_i = chain(:,end-(i-1));
    y_CI(:,i) = f_mod(par_i); % Generate model response
    y_PI(:,i) = y_CI(:,i)+normrnd(0,sqrt(s2(end-(i-1)))); % perturb with measurement noise estimate
end
% Get quantiles
yCI_quant = quantile(y_CI',[0.025 0.5 0.975]);
yPI_quant = quantile(y_PI',[0.025 0.5 0.975]);

%%
% This is for plotting UQ bands in gray
which_t_slice = 8;
t_fill = [tspace fliplr(tspace)];
figure(1);clf; hold on;
fill(t_fill,[yPI_quant(3,:) fliplr(yPI_quant(1,:))],[0.55 0.55 0.55],'EdgeColor','none');
fill(t_fill,[yCI_quant(3,:) fliplr(yCI_quant(1,:))],[0.8 0.8 0.8],'EdgeColor','none');
plot(tspace,yCI_quant(2,:),'r','LineWidth',2);
plot(tspace,data,'ko','LineWidth',2);
xline(tspace(which_t_slice),'--k','LineWidth',3)
set(gca,'FontSize',20); grid on;

figure(2);clf; hold on;
[f,x] = ksdensity(y_PI(which_t_slice,:));
plot(x,f,'Linewidth',3);
f_quant = quantile(y_PI(which_t_slice,:),[0.025 0.975]);
xline(f_quant,'--k','LineWidth',3);
title(sprintf('Posterior Predictive PDF at t=%f',tspace(which_t_slice)))

%%


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