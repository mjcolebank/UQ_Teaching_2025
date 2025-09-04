% Bayes_Quad_Example.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
clear; clc; close all;
%% Compute the posterior distribution under the assumption of a uniform
% prior using the nonlinear Sine model we investigated before
xspace = linspace(-5,5,20 );

sin_model = @(q,x) q(1)+q(2).*sin(pi.*q(3).*x);

param_star = [5 1 0.3];

num_param = 3;
n_xpts = length(xspace);

%% Things that can be adjusted
% Gaussian prior terms
prior_mu = param_star.*1.1;
prior_var = param_star.*0.1;
% Measurement noise
meas_var = 0.5;

%% First, construct the prior.
% We start by assuming uniform (consider Gaussian) on a +-20% range
upper = param_star.*1.2; % Twenty percent above
lower = param_star.*0.8; % Twenty percent above

prior_unif = @(param,a,b) 1./prod(b-a); % The prior PDF is constant throughout the domain (equiprobable outcomes

% Second case, independent gaussian distributions with variance=10%

prior_std = sqrt(prior_var); % Standard deviation
prior_covar = eye(num_param,num_param).*prior_var; % Use a covariance since this is a 3D problem
prior_gauss = @(param,mu,covar) exp(-(param-mu)*inv(covar)*(param-mu)')./sqrt(2.*pi.*det(covar));

%% First, plot prior
n_parpts = 100;
% Uniform
 figure(1); clf; hold on;
for i=1:num_param
    param_range = linspace(lower(i),upper(i),n_parpts);
      prior_eval = zeros(n_parpts,1);

    for j=1:n_parpts
        prior_eval(j) = prior_unif(param_range(j),lower(i),upper(i));
    end
        subplot(2,2,i); hold on;
        plot(param_range,prior_eval,'--','LineWidth',2);
        xline(param_star(i),'--k')
end

% Gaussian
 figure(2); clf; hold on;
for i=1:num_param
    param_range = linspace(prior_mu(i)-3.*prior_std(i),prior_mu(i)+3.*prior_std(i),n_parpts);
   prior_eval = zeros(n_parpts,1);
   param_eval = param_star;
    
    for j=1:n_parpts
        param_eval(i) = param_range(j);
        prior_eval(j) = prior_gauss(param_eval,prior_mu,prior_covar);
    end
        subplot(2,2,i); hold on;
        plot(param_range,prior_eval,'--','LineWidth',2);
        xline(param_star(i),'--k')
end

%% Now, build a likelihood based on the assumption of iid, normal mean zero errors

f_mod = @(param) sin_model(param,xspace);
true_signal = f_mod(param_star);
data        = true_signal + normrnd(0,sqrt(meas_var),1,n_xpts);

figure(99);clf;hold on;
plot(true_signal,'r');
plot(data,'ko')
%%
lik_func = @(param,data) exp(-sum((f_mod(param)-data).^2) ./ (2.*meas_var)).* (2.*pi.*meas_var).^(-n_xpts/2);

% plot likelihood by fixing all parameters at true value, except for
% parameter i (similar to profile likelihood)
 figure(3); clf; hold on;
for i=1:num_param
    param_range = linspace(prior_mu(i)-3.*prior_std(i),prior_mu(i)+3.*prior_std(i),n_parpts);
   lik_eval_true = zeros(n_parpts,1); % True signal, no noise
   lik_eval_data = zeros(n_parpts,1); % Noisy data
   
   param_eval = param_star;
    for j=1:n_parpts
        param_eval(i) = param_range(j);
        lik_eval_true(j) = lik_func(param_eval,true_signal);
        lik_eval_data(j) = lik_func(param_eval,data);
    end
        subplot(2,2,i); hold on;
        plot(param_range,lik_eval_true,'b','LineWidth',2);
        ylabel('likelihood - true signal')
        yyaxis('right')
        plot(param_range,lik_eval_data,'-r','LineWidth',2);
        ylabel('likelihood - noisy data')
        
        xline(param_star(i),'--k');
end
legend('True','Noisy')
%% Now, consider a quadrature approach to constructing the posterior
% Use trapezoidal rule to start
n_quad = 20;
x1range = linspace(prior_mu(1)-2.*prior_std(1),prior_mu(1)+2.*prior_std(1),n_quad);
x2range = linspace(prior_mu(2)-2.*prior_std(2),prior_mu(2)+2.*prior_std(2),n_quad);
x3range = linspace(prior_mu(3)-2.*prior_std(3),prior_mu(3)+2.*prior_std(3),n_quad);



[X1unif,X2unif,X3unif] = meshgrid(x1range,x2range,x3range);
evid_quad_unif = zeros(n_quad,n_quad,n_quad);
evid_quad_gauss = zeros(n_quad,n_quad,n_quad);
lik_eval = zeros(n_quad,n_quad,n_quad);


% A triple for-loop for the three dimensions (yikes)
for i1=1:n_quad
    for i2=1:n_quad
        for i3=1:n_quad
            param_eval = [X1unif(i1,i2,i3),X2unif(i1,i2,i3),X3unif(i1,i2,i3)];
            lik_eval(i1,i2,i3) = lik_func(param_eval,data);
            evid_quad_unif(i1,i2,i3) =  lik_eval(i1,i2,i3).*prior_unif(param_eval,lower,upper); % Equal probability everywhere
            evid_quad_gauss(i1,i2,i3) = lik_eval(i1,i2,i3) .*prior_gauss(param_eval,prior_mu,prior_covar);
        end
    end
end
% Now, compute the evidence in the denominator

post_unif = evid_quad_unif./trapz(x3range,trapz(x2range,trapz(x1range,evid_quad_unif,1),2));

post_gauss = evid_quad_gauss./trapz(x3range,trapz(x2range,trapz(x1range,evid_quad_gauss,1),2));


%% Now look at marginal posteriors
marg_post_unif1 = squeeze(trapz(x2range,trapz(x3range,post_unif,3),2));%./evidence_tot_unif;
marg_post_unif2 = squeeze(trapz(x1range,trapz(x3range,post_unif,3),1));
marg_post_unif3 = squeeze(trapz(x1range,trapz(x2range,post_unif,2),1));

marg_post_gauss1 = squeeze(trapz(x2range,trapz(x3range,post_gauss,3),2));%./evidence_tot_gauss;
marg_post_gauss2 = squeeze(trapz(x1range,trapz(x3range,post_gauss,3),1));
marg_post_gauss3 = squeeze(trapz(x1range,trapz(x2range,post_gauss,2),1));

figure(4);clf;

subplot(2,2,1); hold on;
plot(x1range,marg_post_unif1,'-b','LineWidth',2)
plot(x1range,marg_post_gauss1,'-.r','LineWidth',2)
xline(param_star(1),'--k','LineWidth',1);

subplot(2,2,2); hold on;
plot(x2range,marg_post_unif2,'-b','LineWidth',2)
plot(x2range,marg_post_gauss2,'-.r','LineWidth',2)
xline(param_star(2),'--k','LineWidth',1);

subplot(2,2,3); hold on;
plot(x3range,marg_post_unif3,'-b','LineWidth',2)
plot(x3range,marg_post_gauss3,'-.r','LineWidth',2)
xline(param_star(3),'--k','LineWidth',1);
legend('Uniform','Gaussian')
%%
figure(5);clf
subplot(2,2,1); hold on;
contourf(x1range,x2range,squeeze(trapz(x3range,post_unif ,3)),20, 'LineColor', 'none');
plot(param_star(1),param_star(2),'*r','LineWidth',1);
axis tight
title('Gaussian posterior 2D Contours')
subplot(2,2,2); hold on;
contourf(x2range,x3range,squeeze(trapz(x1range,post_unif ,1)),20, 'LineColor', 'none');
plot(param_star(2),param_star(3),'*r','LineWidth',1);

subplot(2,2,3); hold on;
contourf(x1range,x3range,squeeze(trapz(x2range,post_unif ,2)),20, 'LineColor', 'none');
plot(param_star(1),param_star(3),'*r','LineWidth',1);
colorbar;



figure(6);clf
subplot(2,2,1); hold on;
contourf(x1range,x2range,squeeze(trapz(x3range,post_gauss ,3)),20, 'LineColor', 'none');
plot(param_star(1),param_star(2),'*r','LineWidth',1);
axis tight
title('Gaussian posterior 2D Contours')
subplot(2,2,2); hold on;
contourf(x2range,x3range,squeeze(trapz(x1range,post_gauss ,1)),20, 'LineColor', 'none');
plot(param_star(2),param_star(3),'*r','LineWidth',1);

subplot(2,2,3); hold on;
contourf(x1range,x3range,squeeze(trapz(x2range,post_gauss ,2)),20, 'LineColor', 'none');
plot(param_star(1),param_star(3),'*r','LineWidth',1);
colorbar;
