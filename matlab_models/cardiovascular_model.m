% FUNCTION LV_Ao.m
% This is an introductory function example for the course on cardiovascular
% physiology: theory and practice.
%
% Author: Mitchel J. Colebank
% Date:5/17/2022


function dY = cardiovascular_model(t,y,param)
% Differential equations
Vlv = y(1);
Vao = y(2);

% Unpack parameters
P_LeftAtr = param(1);
P_SysCap  = param(2);
Rmv       = param(3);
Rav       = param(4);
Rsysart   = param(5);
Emax      = param(6);
Emin      = param(7);
T_peak    = param(8);
T_relax   = param(9);
T         = param(10);
Vlv_d     = param(11);
Cao       = param(12);

% Get LV and Ao pressure
plv = LinearElastance(Vlv,t,[Emax,Emin,Vlv_d,T_peak,T_relax,T]);
pao = Vao./Cao;

% Flows through valves and systemic arteries
qmv = max(((P_LeftAtr-plv)./Rmv),0);
qav = max(((plv-pao)./Rav),0);


qcap = (pao-P_SysCap)./Rsysart;

% Update differential equations
dVlv = qmv - qav;
dVsa = qav - qcap;
dY = [dVlv; dVsa];


end


