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
param_star = [z0;v0;C;K;F0;omegaF]; % Get rid of mass term

tend = 20;
tspace = linspace(0,tend,30);
num_param = length(param_star);
n_xpts = length(tspace);


f_mod = @(q) call_spring_model(q,tspace);
true_signal = f_mod(param_star);
%% Things that can be adjusted
% Measurement noise
% noise_var = 0.01;
noise_var = 0.05;
% noise_var = 1.5;

% MCMC algorithm parameters
M = 10000;
theta0 = param_star.*1.00;

data        = true_signal + normrnd(0,sqrt(noise_var),1,n_xpts);
%% Constant covariance
s_theta = 1e-2.*ones(num_param,1); %1e-3, 0.1
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
% We start by assuming uniform (consider Gaussian) on a +-20% range
UB_uni = [5; 5; 4; 4; 4; 4]; % Twenty percent above
LB_uni = [-5; -5; 0.1; 0.1; 0.1; 0.1]; % Twenty percent above

prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes
prior_F = @(param) prior_unif(param,LB_uni,UB_uni);

%% Chain when covariance is not from sensitivity
[chain,s2] = Metropolis_errupdate(f_mod,data,prior_F,theta0,M,covar,UB_uni,LB_uni);


%%
figure(1);clf
subplot(2,3,1); hold on;
plot(chain(1,:));
yline(param_star(1),'--k','LineWidth',2);
subplot(2,3,2); hold on;
plot(chain(2,:));
yline(param_star(2),'--k','LineWidth',2);
subplot(2,3,3); hold on;
plot(chain(3,:));
yline(param_star(3),'--k','LineWidth',2);
subplot(2,3,4); hold on;
plot(chain(4,:));
yline(param_star(4),'--k','LineWidth',2);
subplot(2,3,5); hold on;
plot(chain(5,:));
yline(param_star(5),'--k','LineWidth',2);
subplot(2,3,6); hold on;
plot(chain(6,:));
yline(param_star(6),'--k','LineWidth',2);
%%
figure(2);clf
subplot(2,3,1); hold on;
ksdensity(chain(1,:));
xline(param_star(1),'--k','LineWidth',2);
subplot(2,3,2); hold on;
ksdensity(chain(2,:));
xline(param_star(2),'--k','LineWidth',2);
subplot(2,3,3); hold on;
ksdensity(chain(3,:));
xline(param_star(3),'--k','LineWidth',2);
subplot(2,3,4); hold on;
ksdensity(chain(4,:));
xline(param_star(4),'--k','LineWidth',2);
subplot(2,3,5); hold on;
ksdensity(chain(5,:));
xline(param_star(5),'--k','LineWidth',2);
subplot(2,3,6); hold on;
ksdensity(chain(6,:));
xline(param_star(6),'--k','LineWidth',2);
%%

%% Plot pairwise and get MAP points
MAP_const = zeros(num_param,1);
MAP_est = zeros(num_param,1);

%% Plot pairwise and get MAP points
MAP_const = zeros(num_param,1);
MAP_est = zeros(num_param,1);

figure(5); clf;
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


function y = call_spring_model(q,tspace)
IC = q(1:2);
param = q(3:end);
[~,y] = ode45(@spring_model_4par,tspace,IC,[],param);
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