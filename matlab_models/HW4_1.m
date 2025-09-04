% HW4_1 script
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: April, 2025
% PLEASE EDIT THIS CODE FOR YOUR HOMEWORK; THIS WILL NOT RUN ON ITS OWN
clear; clc; close all;
%% Compute the posterior distribution under the assumption of a uniform
% DEFINE THE INITIAL CONDITIONS HERE

% Five parameter non-homogeneous spring problem
% DEFINE THE FIVE PARAMETERS OF THE MODEL HERE

% PUT ALL THE PARAMETERS INTO THIS VECTOR
param_star = [X;X;X;X]; % Get rid of mass term

% CHANGE THE NOISE VARIANCE HERE
noise_var = XXX;

% CHANGE THE COVARIANCE PROPOSAL MAGNITUDE HERE (s_theta)
s_theta = XXX.*ones(num_param,1);

% CHANGE THE UPPER AND LOWER BOUNDS HERE
UB_uni = XXX;
LB_uni = XXX;


%% These things do not need to change
tend = 20;
tspace = linspace(0,tend,30);
num_param = length(param_star);
n_xpts = length(tspace);

% This defines the QoI, which returns just one entry of the spring problem
% solution, which has to come from an ODE solver
f_mod = @(q) call_spring_model(q,IC,tspace);
true_signal = f_mod(param_star);
M = 10000;
theta0 = param_star.*1.00;
data        = true_signal + normrnd(0,sqrt(noise_var),1,n_xpts);

%% Constant covariance
covar = eye(num_param).*s_theta;
%% Prior distribution
prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes
prior_F = @(param) prior_unif(param,LB_uni,UB_uni);

% Run MCMC
[chain,s2] = Metropolis_errupdate(f_mod,data,prior_F,theta0,M,covar,UB_uni,LB_uni);


%% PLOTTING
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

%%

%% Plot pairwise and get MAP points
MAP_const = zeros(num_param,1);
MAP_est = zeros(num_param,1);

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

figure(4);clf; hold on;
plot(tspace,f_mod(MAP_const),'LineWidth',2);
plot(tspace,data,'ko');
plot(tspace,true_signal,'--k','LineWidth',2)
legend('MAP','Data','Truth')
%%


function y = call_spring_model(q,IC,tspace)
[~,y] = ode45(@spring_model_5par,tspace,IC,[],q);
y = y(:,1)';
end

function dy = spring_model_5par(t,y,params)

% Unpack parameters
m      = params(1);
c      = params(2);
k      = params(3);
f0     = params(4);
omegaF = params(5);

% Redefine state variables for convenience
x = y(1);
v = y(2);


% RHS equations
accel_eq = (- c.*v - k.*x)./m + f0.*cos(omegaF.*t)/m;
dy = [v;
      accel_eq];

end