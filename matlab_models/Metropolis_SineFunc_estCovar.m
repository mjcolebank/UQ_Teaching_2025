% Metropolis_SineFunc_estCovar.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
clear; clc; close all;
rng(249)
%% Compute the posterior distribution under the assumption of a uniform
% prior using the nonlinear Sine model we investigated before
xspace = linspace(-5,5,20 );
sin_model = @(q,x) q(1)+q(2).*sin(pi.*q(3).*x);
param_star = [5 1 0.3];
num_param = 3;
n_xpts = length(xspace);
f_mod = @(param) sin_model(param,xspace);
true_signal = f_mod(param_star);
%% Things that can be adjusted
% Measurement noise
% noise_var = 0.01;
noise_var = 0.6;
% noise_var = 1.5;

% MCMC algorithm parameters
M = 10000;
theta0 = param_star.*0.8

data        = true_signal + normrnd(0,sqrt(noise_var),1,n_xpts);

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


s_theta = 0.1.*ones(3,1); %1e-3, 0.1
covar_const = eye(num_param).*s_theta;
%% First, construct the prior.
% We start by assuming uniform on a +-20% range
UB_uni = [10; 4; 1]; % Twenty percent above
LB_uni = [1; 0.1; 0.01]; % Twenty percent above

prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes

prior_F = @(param) prior_unif(param,LB_uni,UB_uni);

%% Chain when covariance is not from sensitivity
[chain_covar_const] = Metropolis_Algorithm(f_mod,data,prior_F,theta0,s2_est,M,covar_const,UB_uni,LB_uni);

%% Chain when covariance is from sensitivity
[chain_covar_est] = Metropolis_Algorithm(f_mod,data,prior_F,theta0,s2_est,M,covar_est,UB_uni,LB_uni);
%%
figure(1);clf
subplot(2,3,1); hold on;
plot(chain_covar_const(1,:));
yline(param_star(1),'--k','LineWidth',2);
subplot(2,3,2); hold on;
plot(chain_covar_const(2,:));
yline(param_star(2),'--k','LineWidth',2);
subplot(2,3,3); hold on;
plot(chain_covar_const(3,:));
yline(param_star(3),'--k','LineWidth',2);

subplot(2,3,4); hold on;
plot(chain_covar_est(1,:));
yline(param_star(1),'--k','LineWidth',2);
subplot(2,3,5); hold on;
plot(chain_covar_est(2,:));
yline(param_star(2),'--k','LineWidth',2);
subplot(2,3,6); hold on;
plot(chain_covar_est(3,:));
yline(param_star(3),'--k','LineWidth',2);
%%
figure(2);clf
subplot(1,3,1); hold on;
ksdensity(chain_covar_const(1,:));
ksdensity(chain_covar_est(1,:));
xline(param_star(1),'--k','LineWidth',2);
subplot(1,3,2); hold on;
ksdensity(chain_covar_const(2,:));
ksdensity(chain_covar_est(2,:));
xline(param_star(2),'--k','LineWidth',2);
subplot(1,3,3); hold on;
ksdensity(chain_covar_const(3,:));
ksdensity(chain_covar_est(3,:));
xline(param_star(3),'--k','LineWidth',2);

%%
figure(3);clf
subplot(1,3,1); hold on;
ksdensity(chain_covar_const([1 2],:)'); view(2);
subplot(1,3,2); 
ksdensity(chain_covar_const([1 3],:)'); view(2)
subplot(1,3,3);
ksdensity(chain_covar_const([2 3],:)'); view(2)

figure(4);clf
subplot(1,3,1); 
ksdensity(chain_covar_est([1 2],:)'); view(2)
subplot(1,3,2); 
ksdensity(chain_covar_est([1 3],:)'); view(2)
subplot(1,3,3);
ksdensity(chain_covar_est([2 3],:)'); view(2)
%%
clc
figure(5); clf;
for i=1:3
    for j=i+1:3
        disp(2*(i-1)+j);
        subplot(2,2,2*(i-1)+j-1); hold on;
        plot(chain_covar_const(j,:),chain_covar_const(i,:),'o');
    end
    
end