% run_HIVmodel.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the HIV model

clear; clc; close all;
%%
%
%  Set parameters and initial conditions
%
% Note: We split the coefficients between params, which is passed to ode45
%       and coef which is global.  The former comprises the set employed in 
%       Chapters 12 and 13 for Bayesian inference and uncertainty propagation.
%
Names = {'$T_1$','$T_2$','$T^*_1$','$T^*_2$','$V$','$E$'};
  lam1 = 1e+4; % T1 cell target rate
  d1 = .01;    % T1 cell death rate
  epsilon = 0; % T1 treatment efficacy
  k1 = 8.0e-7; %T1 infection rate
  lam2 = 31.98; %T2 cell target rate
  d2 = 0.01; % T2 cell death rate
  f = 0.34; %Treatment efficacy reduction
  k2 = 1e-4; % T2 cell infection rate
  delta = 0.7; %Infected cell death rate
  m1 = 1.0e-5; %T1 clearance rate
  m2 = 1.0e-5; % T2 clearance rate
  NT = 100; %Virions produced per infected cell
  c = 13; % Virus natural death rate
  rho1 = 1; %average virons infecting T1 cell
  rho2 = 1; %average virons infecting T2 cell
  lamE = 1; %Immune effector production rate
  bE = 0.3; % Max birth rate of immune effectors
  Kb = 100; % Birth constant for immune effectors
  dE = 0.25; % Max death rate of immune effectors
  Kd = 500; %Death constant for immune effectors
  deltaE = 0.1; % Natrual death rate of immune effectors

  parameters = [epsilon;k1;lam2;d2;f;NT;c;rho1;rho2;lamE;...
            dE;deltaE;m1;m2;Kd;bE; delta; d1; k2; lam1; Kb];

  Y0 = [0.9e+6; 4000; 1e-1; 1e-1; 1; 12];

%
% Compute and plot the Effector cell population on the interval [0,200] to illustrate 
% the dynamics for the nominal parameter values
%

  t_data = 0:.1:200;
  ode_options = odeset('RelTol',1e-6);
  [t,Y] = ode15s(@HIV_model,t_data,Y0,ode_options,parameters);
  E = Y(:,6);
%%
  figure(1);clf
  for i=1:6
      subplot(2,3,i);
  plot(t_data,Y(:,i),'-','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time (Days)')
  ylabel(Names{i},'Interpreter','latex')
  end
%%
  figure(2)
  plot(t_data,E,'-k','linewidth',3)
  hold on
  plot(t_data(501),E(501),'ob','linewidth',5)
  hold off
  set(gca,'Fontsize',[22]);
  xlabel('Time (Days)')
  ylabel('Immune Effector Response E')

%%
%

% Compute the Effector cell population at time T = 50 for the values in lamE_vec
% and rho1
%

  per = 0.1;
  N_pert = 20;
  lamE_low = lamE-lamE*per;
  lamE_high = lamE+lamE*per;
  lamE_step = (lamE_high-lamE_low)/N_pert;
  lamE_vec = lamE_low:lamE_step:lamE_high;

  rho1_low = rho1-rho1*per;
  rho1_high = rho1+rho1*per;
  rho1_step = (rho1_high-rho1_low)/N_pert;
  rho1_vec = rho1_low:rho1_step:rho1_high;

  tvals = 0:.01:50;
  N = length(tvals);
  for j=1:N_pert+1
    lamE = lamE_vec(j);
    coef = [epsilon,k1,lam2,d2,f,NT,c,rho1,rho2,lamE,dE,deltaE,m1,m2,Kd];
    parameters = [bE; delta; d1; k2; lam1; Kb];
    [t,Y] = ode15s(@HIV_model,tvals,Y0,ode_options,parameters);
    V_lamE(j) = Y(N,5);
    E_lamE(j) = Y(N,6);
  end
  lamE = 1;
  for j=1:N_pert+1
    rho1 = rho1_vec(j);
    coef = [epsilon,k1,lam2,d2,f,NT,c,rho1,rho2,lamE,dE,deltaE,m1,m2,Kd];
    parameters = [bE; delta; d1; k2; lam1; Kb];
    [t,Y] = ode15s(@HIV_model,tvals,Y0,ode_options,parameters);
    V_rho1(j) = Y(N,5);
    E_rho1(j) = Y(N,6);
  end
  
  figure(3)
  plot(rho1_vec,E_rho1,'-k',lamE_vec,E_lamE,'--b','linewidth',3);
  axis([0.9 1.1 34 42])
  set(gca,'Fontsize',[23]);
  xlabel('Parameter Value')
  ylabel('Immune Effector Response E')
  legend('\rho_1','\lambda_E','location','SouthEast')

  figure(4)
  plot(rho1_vec,V_rho1,'-k',lamE_vec,V_lamE,'--b','linewidth',3);
  axis([0.9 1.1 1.525e+5 1.55e+5])
  set(gca,'Fontsize',[23]);
  xlabel('Parameter Value')
  ylabel('Free Virus V')
  legend('\rho_1','\lambda_E','location','SouthEast')


%
%  End hiv_model
%
