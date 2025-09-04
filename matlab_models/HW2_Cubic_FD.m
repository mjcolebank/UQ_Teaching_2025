% cubic_local_sensitivity.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
clear; clc; close all;

% INSERT SPACE HERE
% xspace = ;

% Function
cubic_model = @(q) q(1) + q(2).*xspace ...
            + q(3).*xspace.^2 + q(4).*xspace.^3;

% INSERT PARAMETER VALUE HERE
% param_star = ;

num_param = 4;
n_xpts = length(xspace);
%% Plot function in output space
figure(1); clf;
plot(xspace,cubic_model(param_star));

%% Loop over the parameter values and perturb one element at a time
f_base = cubic_model(param_star); % Reference solutuion
f_sens = zeros(num_param,n_xpts); % Declare storage
h = 1e-6; % STEP SIZE
for i=1:num_param
    param_step = param_star;
    param_step(i) = param_star(i)+h;
    f_step = cubic_model(param_step);
    f_sens(i,:) = (f_step - f_base)./h;
end

