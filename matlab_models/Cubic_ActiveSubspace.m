% OLS_Cubic.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines a cubic model and then calibrates to data under the assumption of
% normal, gaussian iid errors
clear; clc; close all;
%%
% Region where parameters are identifiable.
xspace = linspace(-5,5,20 )';
% Region where parameters are NOT identifiable.
% xspace = linspace(-0.01,0.01,20 )';
% xspace = linspace(-40,40,20 )';

cubic_model = @(q,x) q(1) + q(2).*x ...
    + q(3).*x.^2 + q(4).*x.^3 ;
param_star = [5; 0.1; 0.05; 0.01];
A = [ones(length(xspace),1) xspace xspace.^2 xspace.^3];

% cubic_model = @(q,x) q(1) + q(2).*x ...
    % + q(3).*x.^2 + q(4).*x.^3 + q(5).*x.^3;
% param_star = [5; 0.1; 0.05; 0.05; 0.01];
% A = [ones(length(xspace),1) xspace xspace.^2 xspace.^3 xspace.^3];


% We can also write this as a linear model
[U,S,V] = svd(A);

figure(1); clf; 
for i=1:length(param_star)
    r = i;
    subplot(2,3,i); hold on;
    plot(xspace,A*param_star,'b');

    A_red = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    plot(xspace,A_red*param_star,'--r');
end

figure; plot(diag(S(1:length(param_star),1:length(param_star))),'x','LineWidth',3);
