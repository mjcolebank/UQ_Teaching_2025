% Metropolis_Regression.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
clear; clc; close all;
%% Compute the posterior distribution under the assumption of a uniform
% Initial conditions
% Initial conditions
S0 = 999; I0 = 1; R0  = 0;
IC = [S0;I0;R0];

flag = 0;

% The SIR model has four parameters, as described below
gamma = 0.1; % infection coefficient, [0,1]
k     = 0.1; % contact rate for infection, [0,1]
r     = 0.8; % recovery rate, [0,1]
mu    = 0.1; % combined birth and death rate, [0,1]

tend = 10;
tspace = linspace(0,tend,30);

if flag==0
    param_star = [gamma, k, r, mu]; % Get rid of mass term
    UB_uni = [1,1,1,1]; % Twenty percent above
    LB_uni = [0,0,0,0]; % Twenty percent above
    f_mod = @(q) call_SIR_model(q,IC,tspace);
else
    param_star = [gamma*k,r,mu];
    f_mod = @(q) call_SIR_model_3par(q,IC,tspace);
    UB_uni = [1,1,1]; % Twenty percent above
    LB_uni = [0,0,0]; % Twenty percent above
end

num_param = length(param_star);
n_xpts = length(tspace);



true_signal = f_mod(param_star);
%% Things that can be adjusted
% Measurement noise
noise_var = 60;

% MCMC algorithm parameters
M = 10000;
theta0 = param_star.*1.00;

data        = true_signal + normrnd(0,sqrt(noise_var),1,n_xpts);
figure; plot(true_signal)
hold on;
plot(data,'ko')

%% Constant covariance
s_theta = 1e-6.*ones(num_param,1); %1e-3, 0.1
covar = eye(num_param).*s_theta;

%% In this example, we provide an "estimated" covariance using sensitivity analysis
% S = zeros(n_xpts,num_param);
% h = 1e-6;
% f_theta0 = f_mod(theta0);
% for i=1:num_param
%     param_step = theta0;
%     param_step(i) = param_step(i)+h;
%     S(:,i) = (f_mod(param_step) - f_theta0)./h;
% end
% F = S'*S;
% s2_est = sum( (data-f_mod(theta0)).^2)./(n_xpts - num_param);
% 
% covar = s2_est.*inv(F);
%% First, construct the prior.

prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes
prior_F = @(param) prior_unif(param,LB_uni,UB_uni);

%% Chain when covariance is not from sensitivity
[chain,s2] = Metropolis_errupdate(f_mod,data,prior_F,theta0,M,covar,UB_uni,LB_uni);


%%
figure(1);clf
for i=1:num_param
subplot(2,3,i); hold on;
plot(chain(i,:));
yline(param_star(i),'--k','LineWidth',2);
end
%%
figure(2);clf
for i=1:num_param
subplot(2,3,i); hold on;
ksdensity(chain(i,:));
xline(param_star(i),'--k','LineWidth',2);
end
%% Plot pairwise and get MAP points
MAP_const = zeros(num_param,1);

figure(3); clf;
for i=1:num_param
    for j=i+1:num_param
        subplot(num_param-1,num_param-1,(num_param-1)*(i-1)+j-1); hold on;
        plot(chain(j,:),chain(i,:),'o');
    end
    % Get density values
    [f,x] = ksdensity(chain(i,:));
    % Find the maximum of the PDF
    [who,where] = max(f);
    % Store it
    MAP_const(i) = x(where);

end
%%

disp('MAP - const');
disp(MAP_const');
disp('True');
disp(param_star');

figure(10);clf; hold on;
plot(tspace,f_mod(MAP_const),'LineWidth',2);
plot(tspace,data,'ko');
plot(tspace,true_signal,'--k','LineWidth',2)
legend('MAP','Data','Truth')
%%

function y = call_SIR_model_IC(q,tspace)
IC = q(1:2);
param = q(3:end);
N = sum(IC);

[~,y] = ode45(@SIR_model,tspace,IC,[],param,N);
y = y(:,2)';
end

function y = call_SIR_model_3par(q,IC,tspace)
param = q;
N = sum(IC);
[~,y] = ode45(@SIR_model_3par,tspace,IC,[],param,N);
y = y(:,2)';
end

function y = call_SIR_model(q,IC,tspace)
param = q;
N = sum(IC);
[~,y] = ode45(@SIR_model,tspace,IC,[],param,N);
y = y(:,2)';
end

function dy = SIR_model(t,y,params,N)

% Unpack parameters
gamma = params(1);
k     = params(2);
r     = params(3);
mu    = params(4);

% Redefine state variables for convenience
S = y(1);
I = y(2);
R = y(3);

% RHS equations
dy = [mu*N - mu*S - gamma*k*I*S;
    gamma*k*I*S - (r + mu)*I;
    r*I - mu*R];

end

function dy = SIR_model_3par(t,y,params,N)

% Unpack parameters
GK = params(1);
r     = params(2);
mu    = params(3);

% Redefine state variables for convenience
S = y(1);
I = y(2);
R = y(3);

% RHS equations
dy = [mu*N - mu*S - GK*I*S;
    GK*I*S - (r + mu)*I;
    r*I - mu*R];

end
