% cubic_local_sensitivity.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% clear; clc; close all;

xspace = linspace(-5,5, 20);
% xspace = linspace(-0.5,0.5, 20);

cubic_model = @(q) q(1) + q(2).*xspace ...
            + q(3).*xspace.^2 + q(4).*xspace.^3;

param_star = [5 0.1 0.05 0.01];
% param_star = [0.1 2 1e-5 1e-5];

num_param = 4;
n_xpts = length(xspace);
%%
figure;
plot(xspace,cubic_model(param_star));

%%
f_base = cubic_model(param_star);
f_sens = zeros(num_param,n_xpts);
h = 1e-16;
for i=1:num_param
    param_step = param_star;
    param_step(i) = param_star(i)+h;
    f_step = cubic_model(param_step);
    f_sens(i,:) = (f_step - f_base)./h;
end

figure;
plot(xspace,f_sens);
legend({'$S_1$','$S_2$','$S_3$','$S_4$'},'Interpreter','latex')
