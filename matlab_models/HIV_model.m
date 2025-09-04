% SIR_model.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Inputs:
% t: time (used in ODE solver)
% y: current states of the system
% params: parameters of the system
% N: the total population considered in the model
%
% Outputs:
% dy: the vector containing RHS equations
% Modified from original script by Ralph Smith at North Carolina State
% University

  function dy = HIV_model(t,y,params)

  T1     = y(1);
  T2     = y(2);
  T1_inf = y(3);
  T2_inf = y(4);
  V      = y(5);
  E      = y(6);

  epsilon = params(1);
  k1     = params(2);
  lam2 = params(3);
  d2 = params(4);
  f = params(5);
  NT = params(6);
  c = params(7);
  rho1 = params(8);
  rho2 = params(9);
  lamE = params(10);
  dE = params(11);
  deltaE = params(12);
  m1 = params(13);
  m2 = params(14);
  Kd = params(15);

  bE = params(16);
  delta = params(17);
  d1 = params(18);
  k2 = params(19);
  lam1 = params(20);
  Kb = params(21);

 

  dy = [lam1 - d1*T1 - (1 - epsilon)*k1*V*T1;
        lam2 - d2*T2 - (1 - f*epsilon)*k2*V*T2;
        (1 - epsilon)*k1*V*T1 - delta*T1_inf - m1*E*T1_inf;
        (1 - f*epsilon)*k2*V*T2 - delta*T2_inf - m2*E*T2_inf;
        NT*delta*(T1_inf+T2_inf) - c*V - ((1 - epsilon)*rho1*k1*T1 + (1 - f*epsilon)*rho2*k2*T2)*V;
        lamE + (bE*(T1_inf+T2_inf)./(T1_inf+T2_inf+Kb))*E - (dE*(T1_inf+T2_inf)./(T1_inf+T2_inf+Kd))*E - deltaE*E];
