% HW3_Morris_Spring
% ORIGINAL AUTHOR: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: February, 2025
%
% Uses the trajectory algorithm to identify Morris' indices
%%
clear; clc; close all;
%%
% Initial conditions
z0 = 2; v0 = -2;
IC = [z0;v0];
tend = 20;
tspace = linspace(0,tend,50);

% The spring model has three parameters
C = 0.66;        % damper friction, kg/s
K = 2;        % spring resistance, kg/s^2
F0 = 1;     % forcing term, kg cm / s^2
omegaF = 0.8; % frequency of forcing term, rad/s

% Identify upper and lower bounds for [gamma, k, r, mu]
parameter_names = {'$C$','$K$','$F_0$','$\omega_F$'};
UB = [3 3 5 1];
LB = [0.01 0.01 0.01 0.01];
UB(4) = 0.1;
num_par = length(UB);

f = @(q) min(call_spring_model(q,IC,tspace));
%%
R =200;           % Number of samples we want
p = num_par;   % Number of parameters
l = 60;                % Number of levels
delta = l./(2*(l-1));  % Step Size
upper = UB(:)';
lower = LB(:)';
d = zeros(R,p); % Store the elementary effects
%% Try to use the randomization algorithm
% Note that all parameters are scaled to be in the range 0,1 and then
% rescaled in the model evaluation.
A = zeros(p+1,p);
for i=1:p
    A(i+1:p+1,i) = ones((p-(i-1)),1);
end
X = zeros((p+1).*R,p.*R);
%% 
F_storage = cell(p+1,R);
qstar = unifrnd(0,1,R,p);
Jp = ones(p+1,p);
J1 = ones(p+1,1);
P = eye(p,p);
UL_MAT = eye(p).*(upper-lower);
func_evals = 1;
for i=1:R
    qcurr = qstar(i,:);
    pm1 = rand(p,1);
    Dstar = eye(p).*(pm1 > 0.5) - eye(p).*(pm1 <= 0.5);
    [who,where] = sort(rand(p,1));
    Pstar = P(where,:);
    Astar = J1*qcurr + (delta./2).*(( (2.*A - Jp)*Dstar + Jp))*Pstar;
    C = J1*(lower) + Astar*UL_MAT;
    fpast = f(C(1,:));
    F_storage{1,i} = fpast; % Store the entire solution in case we want to look at time series
    for j=1:p
        fstep = f(C(j+1,:));
        % Calculate the elementary effect with the QoI
        par_sign = sign(C(j+1,where(j)) - C(j,where(j))); % Determine whether this is a forward or backward difference
        d(i,where(j)) =  par_sign.*(fstep - fpast)./delta; % Elementary effect
        fpast = fstep;
        F_storage{where(j)+1,i} = fpast;
    end
end

mu = mean(d,1);
mu_star= mean(abs(d),1);
sigma = sqrt(sum((d-mu_star).^2,1)./(R-1));
%%
figure(1);clf;hold on;
for i=1:num_par
    plot(mu_star(i),sigma(i),'k.','MarkerSize',20)
    text(mu_star(i).*1.01,sigma(i),parameter_names{i},'FontSize',16,'Interpreter','latex');
end

rank = sqrt(mu_star.^2 + sigma.^2);
figure;
bar(rank);
xticklabels(parameter_names);
%%

function y = call_spring_model(q,IC,tspace)
[~,y] = ode45(@spring_model_4par,tspace,IC,[],q);
y = y(:,1)';
end

function dy = spring_model_4par(t,y,params)

% Unpack parameters
C      = params(1);
K      = params(2);
F0     = params(3);
omegaF = params(4);

% Redefine state variables for convenience
z = y(1);
v = y(2);


% RHS equations
accel_eq = (- C.*v - K.*z) + F0.*cos(omegaF.*t);
dy = [v;
      accel_eq];

end