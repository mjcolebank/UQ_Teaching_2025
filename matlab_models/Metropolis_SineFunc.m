% Metropolis_Regression.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
clear; clc; close all;
%% Compute the posterior distribution under the assumption of a uniform
% prior using the nonlinear Sine model we investigated before
xspace = linspace(-5,5,20);

sin_model = @(q,x) q(1)+q(2).*sin(pi.*q(3).*x);
param_star = [5 1 0.3];
num_param = 3;
n_xpts = length(xspace);

f_mod = @(param) sin_model(param,xspace);
true_signal = f_mod(param_star);

%% Things that can be adjusted
% Gaussian prior terms
prior_mu = param_star.*1.0;
prior_var = param_star.*0.5;
% Measurement noise
% noise_var = 0.01;
noise_var = 0.5;
% noise_var = 1.5;

% Proposal distribution
s_theta = 0.005.*abs(param_star);
covar = eye(num_param).*s_theta;
% covar = [0.05 0.02 0; 0.02 0.01 0; 0 0 0.003];

% MCMC algorithm parameters
M = 10000;
theta0 = param_star.*1.0;
data        = true_signal + normrnd(0,sqrt(noise_var),1,n_xpts);

%% When we don't know s2 apriori (typically the case), we get an estimate
s2_est = sum( (data-f_mod(theta0)).^2)./(n_xpts - num_param);

%% First, construct the prior.
% We start by assuming uniform (consider Gaussian) on a +-20% range
UB_uni = [10; 4; 1]; % Twenty percent above
LB_uni = [1; 0.1; 0.01]; % Twenty percent above
UB_gauss = prior_mu+3.*sqrt(prior_var); % Twenty percent above
LB_gauss = prior_mu-3.*sqrt(prior_var); % Twenty percent above

prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes

% Second case, independent gaussian distributions with variance=10%

prior_std = sqrt(prior_var); % Standard deviation
prior_covar = eye(num_param,num_param).*prior_var; % Use a covariance since this is a 3D problem
prior_gauss = @(param,mu,covar) exp(-(param(:)-mu(:))'*inv(covar)*(param(:)-mu(:)))./sqrt(2.*pi.*det(covar));



%% Uniform first
prior_F = @(param) prior_unif(param,LB_uni,UB_uni);
chain_uni = Metropolis_Algorithm(f_mod,data,prior_F,theta0,s2_est,M,covar,UB_uni,LB_uni);

%% Now Gaussian
prior_F = @(param) prior_gauss(param,prior_mu,prior_covar);
chain_gauss = Metropolis_Algorithm(f_mod,data,prior_F,theta0,s2_est,M,covar,UB_gauss,LB_gauss);
%%
figure(1);clf
subplot(2,3,1); hold on;
plot(chain_uni(1,:));
yline(param_star(1),'--k','LineWidth',2);
subplot(2,3,2); hold on;
plot(chain_uni(2,:));
yline(param_star(2),'--k','LineWidth',2);
subplot(2,3,3); hold on;
plot(chain_uni(3,:));
yline(param_star(3),'--k','LineWidth',2);

subplot(2,3,4); hold on;
plot(chain_gauss(1,:));
yline(param_star(1),'--k','LineWidth',2);
subplot(2,3,5); hold on;
plot(chain_gauss(2,:));
yline(param_star(2),'--k','LineWidth',2);
subplot(2,3,6); hold on;
plot(chain_gauss(3,:));
yline(param_star(3),'--k','LineWidth',2);
%%
figure(2);clf
subplot(1,3,1); hold on;
ksdensity(chain_uni(1,:));
ksdensity(chain_gauss(1,:));
xline(param_star(1),'--k','LineWidth',2);
subplot(1,3,2); hold on;
ksdensity(chain_uni(2,:));
ksdensity(chain_gauss(2,:));
xline(param_star(2),'--k','LineWidth',2);
subplot(1,3,3); hold on;
ksdensity(chain_uni(3,:));
ksdensity(chain_gauss(3,:));
xline(param_star(3),'--k','LineWidth',2);
%%
figure(3);clf
subplot(1,3,1); 
ksdensity(chain_uni([1 2],:)'); view(2)
subplot(1,3,2); 
ksdensity(chain_uni([1 3],:)'); view(2)
subplot(1,3,3);
ksdensity(chain_uni([2 3],:)'); view(2)

figure(4);clf
subplot(1,3,1); 
ksdensity(chain_gauss([1 2],:)'); view(2)
subplot(1,3,2); 
ksdensity(chain_gauss([1 3],:)'); view(2)
subplot(1,3,3);
ksdensity(chain_gauss([2 3],:)'); view(2)
