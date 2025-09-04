% run_HIVmodel.m
% Original Author: Ralph Smith, NC State University
% Editor: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the HIV model

clear; clc; close all;
%%

Names = {'$T_1$','$T_2$','$T^*_1$','$T^*_2$','$V$','$E$'};
lam1 = 1e+4; % T1 cell target rate
d1 = .01;    % T1 cell death rate
epsilon = 0.1; % T1 treatment efficacy
k1 = 8.0e-7; %T1 infection rate
lam2 = 31.98; %T2 cell target rate
d2 = 0.01; % T2 cell death rate
f = 0.34; %Treatment efficacy reduction
k2 = 1e-4; % T2 cell infection rate
delta = 0.7; %Infected cell death rate
m1 = 1.0e-5; %T1 clearance rate
m2 = 1.0e-5; % T2 clearance rate
NT = 100; %Virions produced per infected cell
c = 13; % Virus natural death rate
rho1 = 1; %average virons infecting T1 cell
rho2 = 1; %average virons infecting T2 cell
lamE = 1; %Immune effector production rate
bE = 0.3; % Max birth rate of immune effectors
Kb = 100; % Birth constant for immune effectors
dE = 0.25; % Max death rate of immune effectors
Kd = 500; %Death constant for immune effectors
deltaE = 0.1; % Natrual death rate of immune effectors

parameter_names = {...
    'lam1', ...  % T1 cell target rate
    'd1', ...    % T1 cell death rate
    'epsilon', ... % T1 treatment efficacy
    'k1', ...    % T1 infection rate
    'lam2', ...  % T2 cell target rate
    'd2', ...    % T2 cell death rate
    'f', ...     % Treatment efficacy reduction
    'k2', ...    % T2 cell infection rate
    'delta', ... % Infected cell death rate
    'm1', ...    % T1 clearance rate
    'm2', ...    % T2 clearance rate
    'NT', ...    % Virions produced per infected cell
    'c', ...     % Virus natural death rate
    'rho1', ...  % Average virions infecting T1 cell
    'rho2', ...  % Average virions infecting T2 cell
    'lamE', ...  % Immune effector production rate
    'bE', ...    % Max birth rate of immune effectors
    'Kb', ...    % Birth constant for immune effectors
    'dE', ...    % Max death rate of immune effectors
    'Kd', ...    % Death constant for immune effectors
    'deltaE' ... % Natural death rate of immune effectors
};

theta = [epsilon;k1;lam2;d2;f;NT;c;rho1;rho2;lamE;...
    dE;deltaE;m1;m2;Kd;bE; delta; d1; k2; lam1; Kb];

Y0 = [0.9e+6; 4000; 1e-1; 1e-1; 1; 12];
t_data = 0:6:200;
t_highres = 0:0.1:200;
noise_std = 1;
noise_var = noise_std.^2;
n_data = length(t_data);
%%
theta_star = theta;
% Identify which parameters to infer
q_ids = 1:21;
num_par = length(q_ids);

theta_nominal = theta_star(q_ids);
y_star = call_HIV_model(theta_nominal,q_ids,theta_star,Y0,t_highres);
data  = call_HIV_model(theta_nominal,q_ids,theta_star,Y0,t_data) + normrnd(0,noise_std,n_data,1);

figure(1);clf; hold on;
plot(t_highres,y_star,'k','LineWidth',3);
plot(t_data,data,'ro','LineWidth',2)


%% Run the morris screening algorithm
UB = 1.5.*theta_star;
LB = 0.5.*theta_star;
f = @(q) call_HIV_model(q,q_ids,theta_nominal,Y0,t_data);

% This will return a vector of sensitivites corresponding to the time points in the E state
[mu,mu_star,sigma] = run_morris(f,UB,LB,n_data);
%%
figure(1);clf;
subplot(1,2,1); hold on;
plot(t_data,mu_star','LineWidth',3); ylabel('$\mu^*$');
subplot(1,2,2); hold on;
plot(t_data,sigma','LineWidth',3); ylabel('$\sigma$');
%%
mu_star_scalar = mean(mu_star,2);
sigma_scalar = mean(sigma,2);
rank_scalar = sqrt(mu_star_scalar.^2 + sigma_scalar.^2);
figure(2);clf;hold on;
for i=1:num_par
    plot(mu_star_scalar(i),sigma_scalar(i),'k.','MarkerSize',20)
    text(mu_star_scalar(i).*1.01,sigma_scalar(i),parameter_names{i},'FontSize',16,'Interpreter','latex');
end
figure(3); clf;
h = bar(rank_scalar);
set(gca,'XTick',1:21,'XTickLabel',parameter_names)
%% MCMC
prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes
prior_F = @(param) prior_unif(param,LB,UB);

%% Chain when covariance is from sensitivity
% covar = eye(num_par).*1e-16;%.*(abs(theta_star));
% M = 10000;
% k0 = 1000;
% [chain,s2] = AdaptiveMetropolis_s2update(f,data,prior_F,theta_star,k0,M,covar,UB,LB);
%%






%%
function E = call_HIV_model(q,q_ids,q_nominal,Y0,t_data)
theta = q_nominal;
theta(q_ids) = q; % Only update parameters being inferred
ode_options = odeset('RelTol',1e-6);
[~,Y] = ode15s(@HIV_model,t_data,Y0,ode_options,theta);
E = Y(:,6);
end


%%
function [mu,mu_star,sigma] = run_morris(f,UB,LB,n_xpts)
%%
R =20;           % Number of samples we want
p = length(UB);   % Number of parameters
l = 60;                % Number of levels
delta = l./(2*(l-1));  % Step Size
upper = UB(:)';
lower = LB(:)';
d = zeros(R,p,n_xpts); % Store the elementary effects
%% Try to use the randomization algorithm
% Note that all parameters are scaled to be in the range 0,1 and then
% rescaled in the model evaluation.
A = zeros(p+1,p);
for i=1:p
    A(i+1:p+1,i) = ones((p-(i-1)),1);
end
X = zeros((p+1).*R,p.*R);
%% 
F_storage = cell(p+1,R);
qstar = unifrnd(0,1,R,p);
Jp = ones(p+1,p);
J1 = ones(p+1,1);
P = eye(p,p);
UL_MAT = eye(p).*(upper-lower);
func_evals = 1;
for i=1:R
    qcurr = qstar(i,:);
    pm1 = rand(p,1);
    Dstar = eye(p).*(pm1 > 0.5) - eye(p).*(pm1 <= 0.5);
    [who,where] = sort(rand(p,1));
    Pstar = P(where,:);
    Astar = J1*qcurr + (delta./2).*(( (2.*A - Jp)*Dstar + Jp))*Pstar;
    C = J1*(lower) + Astar*UL_MAT;
    fpast = f(C(1,:));
    F_storage{1,i} = fpast; % Store the entire solution in case we want to look at time series
    for j=1:p
        fstep = f(C(j+1,:));
        % Calculate the elementary effect with the QoI
        par_sign = sign(C(j+1,where(j)) - C(j,where(j))); % Determine whether this is a forward or backward difference
        d(i,where(j),:) =  par_sign.*(fstep - fpast)./delta; % Elementary effect
        fpast = fstep;
        F_storage{where(j)+1,i} = fpast;
    end
end

mu = squeeze(mean(d,1));
mu_star= squeeze(mean(abs(d),1));
sigma = zeros(p,n_xpts);
for i=1:n_xpts
    di = squeeze(d(:,:,i));
    sigma(:,i) = sqrt(sum((di-mu_star(:,i)').^2,1)./(R-1));
end

end