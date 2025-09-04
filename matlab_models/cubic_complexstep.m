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
f_FD = zeros(num_param,n_xpts);
f_CS = zeros(num_param,n_xpts);

h = 1e-18;
for i=1:num_param
    param_step = param_star;
    param_step(i) = param_star(i)+h;
    param_complex = param_star;
    param_complex(i) = param_star(i)+1i*h;
    f_step = cubic_model(param_step);
    f_FD(i,:) = (f_step-f_base)./h;
    f_step = cubic_model(param_complex);
    f_CS(i,:) = imag(f_step)./h;
end

figure; 
for i=1:4
    subplot(2,2,i); hold on;
    plot(xspace,f_FD(i,:),'r');
    plot(xspace,f_CS(i,:),'--b');
end
legend('FD','CS')