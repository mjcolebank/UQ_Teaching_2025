% HW3_Saltelli
% ORIGINAL AUTHOR: Ralph Smith, NC State University
% Editor: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
% Uses the Saltelli algorithm to construct global sensitivity indices
%%
clear; clc; close all;
%% Initial conditions for the SIR model
% Region where parameters are identifiable.
%% STUDENTS: define a function, f, that returns the cubic function
% evaluated at x and q, where q are the parameters
% Region where parameters are identifiable.
f = [];

% STUDENTS: Assign upper and lower bounds
UB = [];
LB = [];
num_par = length(UB);

%% Assume that all parameters are uniform
parameter_names = {'$q_1$','$q_2$','$q_3$','$q_4$'};

M = 10000; % Number of samples for Sobol (think parameter values, NOT function evaluations)
A = zeros(M,4);
B = zeros(M,4);

for i=1:M
    A(i,:) = unifrnd(LB,UB);
    B(i,:) = unifrnd(LB,UB);
end

%% Now generate C, which is the matrix A with one column replaced by B
C_cell = cell(4,1); % I use a cell here, but you could also define four matricies, as below in the NOTE
for i=1:4
    C_cell{i} = A; % Assign to A
    C_cell{i}(:,i) = B(:,i); % Replace with ith column
end

%% Lastly initalize the output storage structure
% Integral of infected
f_A = zeros(M,1); % A matrix
f_B = zeros(M,1); % B matrix
f_C = zeros(M,4); % This will have four columns for each C submatrix


tic
w = waitbar(0,'Please wait...');
for j=1:M
    
    f_A(j,1) = f(A(j,:));
    f_B(j,1) = f(B(j,:));
    % Loop over C matrix
    for k=1:4
        f_C(j,k) = f(C_cell{k}(j,:));
    end
    waitbar(j/M)
end

toc
close(w)
%%
y_D = [f_A; f_B];
var_Y = var(y_D);
S_i = zeros(4,1); ST_i = zeros(4,1);
for k=1:4
    S_i(k) = sum(f_B.*(f_C(:,k)-f_A))/var_Y./M;
    ST_i(k) = sum((f_A - f_C(:,k)).^2)/var_Y./(2*M);
   
end
figure;
bar([S_i ST_i]);
legend('First-order','Total');
title('Integral QoI');
xticklabels(parameter_names);


%% Generate a density looking at how shifting one column of A affects the outputs
figure(100);clf;
for k=1:4
    subplot(2,2,k);
    hold on;
    [fA,xA] = ksdensity(f_A,'Function','pdf');
    [fC,xC] = ksdensity(f_C(:,k),'Function','pdf');
    plot(xA,fA,'b','LineWidth',2);
    plot(xC,fC,'r','LineWidth',2);
    legend('A','C');
    set(gca,'FontSize',20);
    grid on;
    title(parameter_names{k})
end

