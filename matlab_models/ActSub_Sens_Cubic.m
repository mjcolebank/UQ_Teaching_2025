% Morris_Cubic
% ORIGINAL AUTHOR: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: February, 2025
%
% Uses the trajectory algorithm to identify Morris' indices
%%
clear; clc; close all;
%% Initial conditions for the SIR model
% Region where parameters are identifiable.
cubic_model = @(q,x) q(1) + q(2).*x ...
            + q(3).*x.^2 + q(4).*x.^3;

% We will evaluate the model at either x=0.1, x=2, or x=4
X_loc =0.1;
f = @(q) cubic_model(q,X_loc);

param_star = [5 0.1 0.05 0.01];
param_guess = [4 0.2 0.07 0.03];

num_param = 4;

% Identify upper and lower bounds for [gamma, k, r, mu]
parameter_names = {'$q_1$','$q_2$','$q_3$','$q_4$'};
UB = [3 1 0.5 0.1];
LB = [0.1 0.1 0.05 0.01];
num_par = length(UB);

%% First, run active subspace analysis
M = 2000; h = 1e-6;
X = lhsdesign(M,num_par);
par_LHS = LB+(UB-LB).*X;
I_h = eye(num_par).*h.*1i;
grad_f = zeros(num_par,1);
C_f = zeros(num_par,num_par);
for i=1:M
    fpast = f(par_LHS(1,:));
    for j=1:num_par
        fstep = f(par_LHS(i,:)+I_h(j,:));
        % Calculate the elementary effect with the QoI
       
        grad_f(j,:) =  imag(fstep)./h; % Elementary effect
       
    end
    C_f = C_f+grad_f*grad_f';

end
C_f = C_f./M;
[V,D,W] = eig(C_f,'vector');
[D, ind] = sort(D,'descend');
W = W(:, ind);
act_scores = sum(diag(D)*W'.^2);

figure(1); clf;
bar(D);

%% Compare to Morris
 [mu,mu_star,sigma] = call_morris(f,UB,LB);

rank = sqrt(mu_star.^2 + sigma.^2);
figure(2); clf; 
subplot(1,2,1); hold on;
bar([rank; mu_star]);
xticks(1:4)
xticklabels(parameter_names);

subplot(1,2,2);
bar(act_scores);
xticklabels(parameter_names);
xticklabels(parameter_names);

%%



function [mu,mu_star,sigma] = call_morris(f,UB,LB)
%%
R = 100;           % Number of samples we want
p = 4;   % Number of parameters
l = 60;                % Number of levels
delta = l./(2*(l-1));  % Step Size
upper = UB(:)';
lower = LB(:)';
d = zeros(R,p); % Store the elementary effects
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
qstar = unifrnd(0,1-delta,R,p);
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
        d(i,where(j)) =  par_sign.*(fstep - fpast)./delta; % Elementary effect
        fpast = fstep;
        F_storage{where(j)+1,i} = fpast;
    end
end

mu = mean(d,1);
mu_star= mean(abs(d),1);
sigma = sqrt(sum((d-mu_star).^2,1)./(R-1));
end