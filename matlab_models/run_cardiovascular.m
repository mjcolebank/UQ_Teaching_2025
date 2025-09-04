% run_cardiovascular.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the 2-component cardiovascular model.
close all; clear; clc;
%%
% Now define the parameters for the system
P_LA     = 5; % left atrial pressure (mmHg)
P_SysCap = 20; % systemic capillary pressure (mmHg)
Rmv      = 5e-3; % resistance in the mitral valve (mmHg s/micro l)
Rav      = 1e-2; % resistance in the aortic valve (mmHg s/micro l)
Rart     = 1.2;  % resistance of the systemic arteries/arterioles (mmHg s/ml)
Emax     = 1.3;  % End systolic elastance (mmHg/ml)
Emin     = 0.03;  % Nonlinear elastance term (1/ml) 
T_peak   = 0.25; % peak elastance (s)
T_relax  = 0.55;  % end of systole (s)
T        = 1.0;  % total cycle length
Cao      = 1.1;  % aortic compliance (ml s / mmHg)
Vlv_d    = 5;    % volume not ejected in the heart

% Stack all the parameters into a vector
param = [P_LA,P_SysCap,Rmv,Rav,Rart,Emax,Emin,T_peak,T_relax,T,Vlv_d,Cao];

% Volumes are in microliters
Vlv_init = 200; % LV Volume
Vao_init = 380;  % Aortic Volume
% Starting and ending time
tstart = 0;
tend   = 50.*T; % 30 cycles
dt = 1e-3;
tspace = tstart:dt:tend;
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Necessary for the ODE solver


Initial_Conditions = [Vlv_init; Vao_init]; % Our initial conditions
% Now solve our system of ODEs
y = ode45(@cardiovascular_model,[tstart, tend],Initial_Conditions,options,param);


% We solved the model for 30 heartbeats; we want to extract the last two
% for plotting
tplot = linspace(29*T-dt,50*T,1000);
yout = deval(y,tplot);
tplot = tplot-tplot(1);


% Extract the solutions to the two differential equations
Vlv = yout(1,:);
Vao = yout(2,:);
% param = [P_LA,P_SysCap,Rmv,Rav,Rart,Emax,Emin,T_peak,T_relax,T,Vlv_d,Cao];
P_LA     = param(1);
P_SysCap = param(2);
Rmv      = param(3);
Rav      = param(4);
Rart     = param(5);
Cao      = param(12);
% Now, recompute the pressures and flows
plv = LinearElastance(Vlv,tplot,[Emax,Emin,Vlv_d,T_peak,T_relax,T]);
pao = Vao./Cao;

% Use the 'max' operator to keep positive flows for valves
QMV = max((P_LA-plv)./Rmv,0);
QAV = max((plv-pao)./Rav,0);
Qart_sys = (pao-P_SysCap)./Rart;

%% Plotting
% Pressure
figure(1); clf; hold on;
plot(tplot,plv,':c','LineWidth',2);
plot(tplot,pao,':m','LineWidth',2);
yline([120 80],'--k')
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
grid on;
set(gca,'FontSize',20)

%% PV loop
figure(2); hold on;
plot(Vlv,plv,'r','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
grid on;
set(gca,'FontSize',20);

