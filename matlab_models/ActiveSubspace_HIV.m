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


%% Construct avtive subspace representation
UB = 1.2.*theta_star;
LB = 0.8.*theta_star;
f = @(q) call_HIV_model(q,q_ids,theta_nominal,Y0,t_data);

M = 1000;
X = lhsdesign(M,num_par);
h = 1e-6;
par_LHS = LB+(UB-LB).*X';
par_LHS = par_LHS';
pert_dir = rand(M,num_par);
I_h = eye(num_par).*h;
grad_f = zeros(num_par,n_data);
C_f = zeros(num_par,num_par);
for i=1:M
    fpast = f(par_LHS(1,:));
    for j=1:num_par
        fstep = f(par_LHS(j+1,:)+I_h(j,:));
        % Calculate the elementary effect with the QoI
       
        grad_f(j,:) =  (fstep - fpast)./delta; % Elementary effect
       
    end
    C_f = C_f+grad_f*grad_f';

end
%%
C_f = C_f./M;
%%
[V,D,W] = eig(C_f,'vector');
[D, ind] = sort(D,'descend');
W = W(:, ind);
act_scores = sum(diag(D)*W'.^2);
%% Now build active subspace
r = 3;
n_train = 100;
act_subspace = @(theta) W(:,1:r)'*theta;
x_train = zeros(r,n_train);
y_train_8 = zeros(1,n_train);
y_train_15 = zeros(1,n_train);
X_train = lhsdesign(n_train,num_par);
par_train = LB+(UB-LB).*X_train';
for j=1:n_train
    y = f(par_train(:,j));
    y_train_8(j) = y(8);
    y_train_15(j) = y(15);

    x_train(:,j) = act_subspace(par_train(:,j));
end
X_poly = [ones(1,n_train); x_train; x_train.^2; x_train.^3; x_train.^4];
 B = inv(X_poly*X_poly')*X_poly*y_train_8';%regress(y_train_8,X_poly')

 %
 % [who,where] = sort(x_train);
 % y_train_8 = y_train_8(where);
 % for i=1:r
 %     [who,where] = sort(x_train(i,:));
 %     figure; plot(who,y_train_8(where),'o');
 %     hold on;
 %     plot(who,B'*X_poly);
 % end

 figure; plot(y_train_8-B'*X_poly,'o')
%%
function E = call_HIV_model(q,q_ids,q_nominal,Y0,t_data)
theta = q_nominal;
theta(q_ids) = q; % Only update parameters being inferred
ode_options = odeset('RelTol',1e-6);
[~,Y] = ode15s(@HIV_model,t_data,Y0,ode_options,theta);
E = Y(:,6);
end