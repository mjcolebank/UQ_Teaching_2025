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
X_loc =0.5;
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
subplot(1,2,1); bar(D); title('Eigenvalues')
subplot(1,2,2); bar(act_scores); xticks(1:4); xticklabels(parameter_names)

%% Construct active subspace variables
z = @(q) W(:,1)'*q';
n_train = 50;
n_test = 20;
q_train = LB + (UB-LB).*unifrnd(0,1,n_train,num_param);
q_test = LB + (UB-LB).*unifrnd(0,1,n_test,num_param);


f_train = zeros(n_train,1);
f_test = zeros(n_test,1);
z_train = zeros(n_train,1);
z_test = zeros(n_test,1);
for j=1:length(q_train)
    z_train(j) = z(q_train(j,:));
    f_train(j) = f(q_train(j,:));
end
for j=1:length(q_test)
    z_test(j) = z(q_test(j,:));
    f_test(j) = f(q_test(j,:));

end
act_poly = polyfit(z_train,f_train,1);

[z_train,where] = sort(z_train);
f_train = f_train(where);
[z_test,where] = sort(z_test);
f_test = f_test(where);

figure(2); clf; hold on;
plot(z_train,polyval(act_poly,z_train),'--r');
plot(z_test,f_test,'ko');
