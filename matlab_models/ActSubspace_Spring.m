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

% Five parameter non-homogeneous spring problem
% The spring model has three parameters
m = 3;        % mass, kg
c = 2;        % damper friction, kg/s
k = 5;        % spring resistance, kg/s^2
F0 = 1;     % forcing term, kg cm / s^2
omegaF = 0.8; % frequency of forcing term, rad/s
noise_std = 0.15;

% Put all the parameters together
param_star = [m;c;k;F0;omegaF]; % Get rid of mass term
par_names = {'m','c','k','F0','\omega'};
tend = 20;
tspace = linspace(0,tend,30);
num_param = length(param_star);
n_xpts = length(tspace);


f_mod = @(q) call_spring_model(q,IC,tspace);
true_signal = f_mod(param_star);
n_data = length(true_signal);
%% Construct active subspace representation
UB = 1.2.*param_star;
LB = 0.8.*param_star;

M = 2000;
X = lhsdesign(M,num_param);
h = 1e-8;
par_LHS = LB+(UB-LB).*X';
par_LHS = par_LHS';
pert_dir = rand(M,num_param);
I_h = eye(num_param).*h;
grad_f4 = zeros(num_param,n_data);
C_f4 = zeros(num_param,num_param);

grad_f14 = zeros(num_param,n_data);
C_f14 = zeros(num_param,num_param);

grad_f25 = zeros(num_param,n_data);
C_f25 = zeros(num_param,num_param);
for i=1:M
    % fbase = f_mod(par_LHS(1,:));
    for j=1:num_param
        fstepP = f_mod(par_LHS(i,:)+I_h(j,:));  
        fstepM = f_mod(par_LHS(i,:)-I_h(j,:));  
        grad_f4(j,:) =  0.5.*(fstepP(4)-fstepM(4))./h;
        grad_f14(j,:) =  0.5.*(fstepP(14)-fstepM(14))./h;
        grad_f25(j,:) =  0.5.*(fstepP(25)-fstepM(25))./h;
    end
    C_f4 = C_f4+grad_f4*grad_f4';
    C_f14 = C_f14+grad_f14*grad_f14';
    C_f25 = C_f25+grad_f25*grad_f25';

end
C_f4 = C_f4./M;
C_f14 = C_f14./M;
C_f25 = C_f25./M;

[~,D4,W4] = eig(C_f4,'vector');
[D4, ind] = sort(D4,'descend');
W4 = W4(:, ind);
act_scores4 = sum(diag(D4)*W4'.^2);

[~,D14,W14] = eig(C_f14,'vector');
[D14, ind] = sort(D14,'descend');
W14 = W14(:, ind);
act_scores14 = sum(diag(D14)*W14'.^2);

[~,D25,W25] = eig(C_f25,'vector');
[D25, ind] = sort(D25,'descend');
W25 = W25(:, ind);
act_scores25 = sum(diag(D25)*W25'.^2);

%%
figure(1);clf;
subplot(1,2,1);
bar([D4./max(D4) D14./max(D14) D25./max(D25)]);
xticks(1:num_param)
title('Eigenvalues'); grid on;
set(gca,'FontSize',20);

subplot(1,2,2);
bar([act_scores4./max(act_scores4); act_scores14./max(act_scores14); act_scores25./max(act_scores25)]')
xticks(1:num_param)
xticklabels(par_names);
title('Activity Scores'); grid on;
set(gca,'FontSize',20);

%% Now build active subspace
r = 1;
n_train = 100;
act_subspace = @(theta,W) W(:,1:r)'*theta;


x_train_4 = zeros(r,n_train);
x_train_14 = zeros(r,n_train);
x_train_25 = zeros(r,n_train);

y_train_4 = zeros(1,n_train);
y_train_14 = zeros(1,n_train);
y_train_25 = zeros(1,n_train);

lhs_samp = lhsdesign(n_train,num_param);
par_train = LB+(UB-LB).*lhs_samp';
for j=1:n_train
    y = f_mod(par_train(:,j));
    y_train_4(j) = y(4);
    y_train_14(j) = y(14);
    y_train_25(j) = y(25);

    x_train_4(:,j) = act_subspace(par_train(:,j),W4);
    x_train_14(:,j) = act_subspace(par_train(:,j),W14);
    x_train_25(:,j) = act_subspace(par_train(:,j),W25);
end
[x_train_4,sort_ID] = sort(x_train_4);
y_train_4 = y_train_4(sort_ID);

[x_train_14,sort_ID] = sort(x_train_14);
y_train_14 = y_train_14(sort_ID);

[x_train_25,sort_ID] = sort(x_train_25);
y_train_25 = y_train_25(sort_ID);

X_poly4 = [ones(1,n_train); x_train_4; x_train_4.^2; x_train_4.^3; x_train_4.^4];
X_poly14 = [ones(1,n_train); x_train_14; x_train_14.^2; x_train_14.^3; x_train_14.^4];
X_poly25 = [ones(1,n_train); x_train_25; x_train_25.^2; x_train_25.^3; x_train_25.^4];

 B4 = inv(X_poly4*X_poly4')*X_poly4*y_train_4';
 B14 = inv(X_poly14*X_poly14')*X_poly14*y_train_14';
 B25 = inv(X_poly25*X_poly25')*X_poly25*y_train_25';

 figure(3);clf;
 subplot(1,3,1); hold on;
 plot(x_train_4,y_train_4,'ko');
 plot(x_train_4,X_poly4'*B4,'--r');
 title('t=2.0690')
 grid on; set(gca,'FontSize',20);
 
 subplot(1,3,2); hold on;
 plot(x_train_14,y_train_14,'ko');
 plot(x_train_14,X_poly14'*B14,'--r');
  title('t=8.9655')
  grid on; set(gca,'FontSize',20);

 subplot(1,3,3); hold on;
 plot(x_train_25,y_train_25,'ko');
 plot(x_train_25,X_poly25'*B25,'--r');
   title('t=16.5517');
   grid on; set(gca,'FontSize',20);

 
 return
%% First, construct the prior.
% We start by assuming uniform (consider Gaussian) on a +-20% range
UB_uni = [10; 10; 10; 10; 10]; % Twenty percent above
LB_uni = [0.1; 0.1; 0.1; 0.1; 0.1]; % Twenty percent above

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