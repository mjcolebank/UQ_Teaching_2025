% Sobol_Saltelli_SIR
% ORIGINAL AUTHOR: Ralph Smith, NC State University
% Editor: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Uses the Saltelli algorithm to construct global sensitivity indices
%%
clear; clc; close all;
SIR_case = 2; % We can look at the effects of different prior distributions

%% Initial conditions for the SIR model
S0 = 900; I0 = 100; R0 = 0;
X0 = [S0; I0; R0];
N = sum(X0);

% Time variables
t_final = 20;
dt = 0.05;
t_data = 0:dt:t_final;

% We will consider two quantities of interest: the integral of the infected
% state and the steady-state equilibrium value of the infected state
QoI_time = length(t_data);
Y0 = [S0; I0; R0];

%% Assume that mu, eta, and gamma are uniform 
M = 1000; % Number of samples for Sobol (think parameter values, NOT function evaluations)

% Upper and lower bounds for three of the parameters
gamma_range = [0 0.1];
r_range     = [0 1.0];
mu_range    = [0 0.5];

%% We consider two different distributions for the parameter
if SIR_case==1
    alpha = 5;
    beta = 9;
else
    alpha = 0.2;
    beta = 15;
end

%% Generate A and B matrices
% Note that these are the parameters for the system
A(:,1) = unifrnd(gamma_range(1),gamma_range(2),M,1);
A(:,2) = betarnd(alpha,beta,M,1);
A(:,3) = unifrnd(r_range(1),r_range(2),M,1);
A(:,4) = unifrnd(mu_range(1),mu_range(2),M,1);

B(:,1) = unifrnd(gamma_range(1),gamma_range(2),M,1);
B(:,2) = betarnd(alpha,beta,M,1);
B(:,3) = unifrnd(r_range(1),r_range(2),M,1);
B(:,4) = unifrnd(mu_range(1),mu_range(2),M,1);


figure;
plot(0:1e-2:1,betapdf(0:1e-2:1,alpha,beta),'LineWidth',3);
title('Beta distribution for k')
%% Now generate C, which is the matrix A with one column replaced by B
C_cell = cell(4,1); % I use a cell here, but you could also define four matricies, as below in the NOTE
for i=1:4
    C_cell{i} = A; % Assign to A
    C_cell{i}(:,i) = B(:,i); % Replace with ith column
end
%% NOTE: you could also define C in this manner (but this may be cumbersome for large dimensions)
% C1 = A; C2 = A; C3 = A; C4 = A;
% C1(:,1) = B(:,1); C2(:,2) = B(:,2); C3(:,3) = B(:,3); C4(:,4) = B(:,4);

%% Lastly initalize the output storage structure
% Integral of infected
I_int_A = zeros(M,1); % A matrix
I_int_B = zeros(M,1); % B matrix
I_int_C = zeros(M,4); % This will have four columns for each C submatrix

% Steady state infected contribution
I_equi_A = zeros(M,1); % A matrix
I_equi_B = zeros(M,1); % B matrix
I_equi_C = zeros(M,4); % This will have four columns for each C submatrix

%%
ode_options = odeset('RelTol',1e-6);

tic
wait = waitbar(0,'Please wait...');
for j=1:M
    params = A(j,:);
    [~,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    I_int_A(j,1) = sum(dt*Y(:,3));
    I_equi_A(j,1) = Y(QoI_time,3);

    params = B(j,:);
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    I_int_B(j,1) = sum(dt*Y(:,3));
    I_equi_B(j,1) = Y(QoI_time,3);

    % Loop over C matrix
    for k=1:4
        params = C_cell{k}(j,:);
        [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
        I_int_C(j,k) = sum(dt*Y(:,3));
        I_equi_C(j,k) = Y(QoI_time,3);
    end
    waitbar(j/M)
end
toc
%%
y_D_int = [I_int_A; I_int_B];
y_D_equi = [I_equi_A; I_equi_B];
f02_int = ((1/(2*M))^2)*sum(y_D_int)*sum(y_D_int); % Squared expectation
f02_equi = ((1/(2*M))^2)*sum(y_D_equi)*sum(y_D_equi);% Squared expectation

S_i_int = zeros(4,1); ST_i_int = zeros(4,1);
S_i_equi = zeros(4,1); ST_i_equi = zeros(4,1);
for k=1:4
    S_i_int(k) = ((1/M)*(I_int_B'*I_int_C(:,k) - I_int_B'*I_int_A))/((1/(2*M))*(y_D_int'*y_D_int) - f02_int);
    ST_i_int(k) = ((1/(2*M))*(I_int_A'*I_int_A - 2*I_int_A'*I_int_C(:,k) + I_int_C(:,k)'*I_int_C(:,k)))/((1/(2*M))*(y_D_int'*y_D_int) - f02_int);

    S_i_equi(k) = ((1/M)*(I_equi_B'*I_equi_C(:,k) - I_equi_B'*I_equi_A))/((1/(2*M))*(y_D_equi'*y_D_equi) - f02_equi);
    ST_i_equi(k) = ((1/(2*M))*(I_equi_A'*I_equi_A - 2*I_equi_A'*I_equi_C(:,k) + I_equi_C(:,k)'*I_equi_C(:,k)))/((1/(2*M))*(y_D_equi'*y_D_equi) - f02_equi);
end
figure;
bar([S_i_int ST_i_int]);
legend('First-order','Total');
title('Integral QoI');
xticklabels({'\gamma','k','r','\mu'});

figure;
bar([S_i_equi ST_i_equi]);
legend('First-order','Total');
title('Equilibrium QoI');
xticklabels({'\gamma','k','r','\mu'});


%% Generate a density looking at how shifting one column of A affects the outputs
figure(100);clf;
for k=1:4
    subplot(2,2,k);
    hold on;
    [fA,xA] = ksdensity(I_int_A,'Function','pdf');
    [fC,xC] = ksdensity(I_int_C(:,k),'Function','pdf');
    plot(xA,fA,'b','LineWidth',2);
    plot(xC,fC,'r','LineWidth',2);
    legend('A','C');
    set(gca,'FontSize',20);
    grid on;
    xlim([0 1000])
end

