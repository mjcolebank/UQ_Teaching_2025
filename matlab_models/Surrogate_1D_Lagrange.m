% Surrogate_1D_Lagrange.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: April, 2025
%
clear; clc; close all;
%% Try to build a surrogate for the following function
% using Lagrange polynomials
f_HF = @(q) (6.*q-2).^2 .* sin(12.*q - 4);
q_space = linspace(0,1,100);
% q_data = 0.1:0.1:0.9;
q_data = [0 0.1:0.1:0.9 1];

y_data = f_HF(q_data);

%% Now build a lagrange polynomial
M = length(y_data);
lag = @(x) 0;
for i=1:M
    f = make_lagrange(q_data,i,M);
    lag = @(x) lag(x) + y_data(i).*f(x);
end
f_LF = lag;
%% Evaluate the LF and HF models
figure(1);clf; hold on;
plot(q_space,f_HF(q_space),'k','LineWidth',3);
plot(q_space,f_LF(q_space),'--r','LineWidth',3);
plot(q_data,y_data,'ko');
grid on; set(gca,'FontSize',20);
legend('HF','LF','Data')
%%
function f = make_lagrange(q_points,m_curr,M)
%% Make a recursive function call for Lagrange Polynomials
%
% q_points = vector of q values from your sampling
% m_curr   = number of values in the q_vector
% q        = current parameter value you are interested in

ids = 1:M;
ids = ids(ids~=m_curr);
q_eval = q_points(ids);
q_m    = q_points(m_curr);
f = @(q,q_eval) 1;
for i=1:M-1
    f = @(q,q_eval) f(q,q_eval) .* (q - q_eval(i))./(q_m - q_eval(i));
end
f = @(q) f(q,q_eval);
end