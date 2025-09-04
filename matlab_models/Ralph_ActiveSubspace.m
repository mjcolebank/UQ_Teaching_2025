%   
%                  SIR_sensitivities.m
%

  clear all

  global N 

%
% Compute the active subspace for the integrated 4-parameter SIR model.  The
% response is taken to be y(\theta) = int_0^5 R(t,\theta). We employ the complex
% step approximation to compute the sensitivities.
%


%
% Set parameters and initial conditions
%

  S0 = 900;
  I0 = 100;
  R0 = 0;
  N = 1000;
  M = 5e+3;

  prob=2;

  tf = 20;
  dt = 0.05;
  t_vals = 0:dt:tf;
  Nt = length(t_vals);
  w = dt*ones(size(t_vals));
  w(1) = dt/2;
  w(Nt) = dt/2; % Missing before

  Y0 = [S0; I0; R0];
  I_delta = 0.5;
  I_gamma = 0.1;
  I_r = 1.2;

  if prob==1
    alpha = 5;
    beta = 9;
  else
    alpha = 0.2;
    beta = 15;
  end

  gamma = 0.2;
  k = 0.1;
  r = 0.6;
  delta = 0.15;
  params = [gamma k r delta];

  rand_params(1,:) = I_gamma*rand(M,1);
  rand_params(2,:) = betarnd(alpha,beta,1,M);
  rand_params(3,:) = I_r*rand(M,1);
  rand_params(4,:) = I_delta*rand(M,1);

  ode_options = odeset('RelTol',1e-8);

%
% Compute the local sensitivity matrix for the Parameter Subset Selection 
% Algorithm 8.13
%

  h = 1e-16;
  gamma_complex = complex(gamma,h);
  params = [gamma_complex k r delta];
  [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
  R_gamma = imag(Y(:,3))/h;

  k_complex = complex(k,h);
  params = [gamma k_complex r delta];
  [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
  R_k = imag(Y(:,3))/h;

  r_complex = complex(r,h);
  params = [gamma k r_complex delta];
  [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
  R_r = imag(Y(:,3))/h;

  delta_complex = complex(delta,h);
  params = [gamma k r delta_complex];
  [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
  R_delta = imag(Y(:,3))/h;

  Sens = [R_gamma R_k R_r R_delta];
  Fisher = Sens'*Sens;

  [V_f,D_f] = eig(Fisher);

  Sens_red = [R_k R_r R_delta];
  Fisher_red = Sens_red'*Sens_red;
  [V_r,D_r] = eig(Fisher_red);

%
% Compute the sensitivities using complex-step derivative approximations.  This
% requires p=3 solutions of the ODE system.
%

  C = zeros(4,4);
  DB = zeros(4,1);

  wait = waitbar(0,'Please wait...');
  for j=1:M
    gamma = rand_params(1,j);
    k = rand_params(2,j);
    r = rand_params(3,j);
    delta = rand_params(4,j);  

    gamma_complex = complex(gamma,h);
    params = [gamma_complex k r delta];
    [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
    R_gamma = imag(Y(:,3))/h;
    RG(:,j) = R_gamma;
    y_gamma = w*R_gamma;

    k_complex = complex(k,h);
    params = [gamma k_complex r delta];
    [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
    R_k = imag(Y(:,3))/h;
    RK(:,j) = R_k;
    y_k = w*R_k;

    r_complex = complex(r,h);
    params = [gamma k r_complex delta];
    [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
    R_r = imag(Y(:,3))/h;
    RR(:,j) = R_r;
    y_r= w*R_r;

    delta_complex = complex(delta,h);
    params = [gamma k r delta_complex];
    [t,Y] = ode45(@SIR_rhs,t_vals,Y0,ode_options,params);
    R_delta = imag(Y(:,3))/h;
    RD(:,j) = R_delta;
    y_delta = w*R_delta;
    
    grad = [y_gamma; y_k; y_r; y_delta];
    grad2 = [y_gamma^2; y_k^2; y_r^2; y_delta^2];    
%    grad2 = [abs(y_gamma); abs(y_k); abs(y_r); abs(y_delta)];

    G(:,j) =  (1/sqrt(M))*grad;
    C = C + (1/M)*grad*grad';
    DB = DB + (1/M)*grad2;
    waitbar(j/M)
  end

%
% Compute the eigenvalues and singular values of the matrices C and G.
%

  tol = 1e-8;
  [U,S,V] = svd(G,'econ');
  eigs = diag(S).^2;

  [U2,S2] = eig(C); 
  eigs2 = diag(S2);

%
% Compute the activity scores.
% 
% Case i: Beta(2,7) The eigenvalues are within 2 orders of magnitude so we 
%                   add the components of all 4
% Case ii: Beta(0.2,15) We keep only the first component
%

  if prob==1
    act_score(1) = eigs(1)*(U(1,1)^2);
    act_score(2) = eigs(1)*(U(2,1)^2);
    act_score(3) = eigs(1)*(U(3,1)^2);
    act_score(4) = eigs(1)*(U(4,1)^2);
    act_score

    act_score2(1) = eigs2(1)*(U2(1,1)^2);
    act_score2(2) = eigs2(1)*(U2(2,1)^2);
    act_score2(3) = eigs2(1)*(U2(3,1)^2);
    act_score2(4) = eigs2(1)*(U2(4,1)^2);
    act_score2

  else
    act_score(1) = eigs(1)*(U(1,1)^2) + eigs(2)*(U(1,2)^2);
    act_score(2) = eigs(1)*(U(2,1)^2) + eigs(2)*(U(2,2)^2);
    act_score(3) = eigs(1)*(U(3,1)^2) + eigs(2)*(U(3,2)^2);
    act_score(4) = eigs(1)*(U(4,1)^2) + eigs(2)*(U(4,2)^2);
    act_score

    act_score2(1) = eigs2(1)*(U2(1,1)^2) + eigs2(2)*(U2(1,2)^2);
    act_score2(2) = eigs2(1)*(U2(2,1)^2) + eigs2(2)*(U2(2,2)^2);
    act_score2(3) = eigs2(1)*(U2(3,1)^2) + eigs2(2)*(U2(3,2)^2);
    act_score2(4) = eigs2(1)*(U2(4,1)^2) + eigs2(2)*(U2(4,2)^2);
    act_score2

  end

  DB
%
% Plot the solutions.
%

  figure(1)
  subplot(2,2,1)
  plot(t,R_gamma,'linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_\gamma')

  subplot(2,2,2)
  plot(t,R_k,'linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_k')

  subplot(2,2,3)
  plot(t,R_r,'linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_r')

  subplot(2,2,4)
  plot(t,R_delta,'linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_\delta')
 
  
  figure(2)
  semilogy(eigs/eigs(1),'bo','MarkerSize',10,'LineWidth',4)
  set(gca,'Fontsize',[22]);
  xlabel('Index')
  ylabel('Normalized Eigenvalues')

  figure(3)
  subplot(2,2,1)
  plot(t,RG')
  subplot(2,2,2)
  plot(t,RK')
  subplot(2,2,3)
  plot(t,RR)
  subplot(2,2,4)
  plot(t,RD)

  figure(4)
  plot(t,RG','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_\gamma')

  figure(5)
  plot(t,RK','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_k')

  figure(6)
  plot(t,RR','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_r')

  figure(7)
  plot(t,RD','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Time')
  ylabel('R_\delta')





%
%          SIR_rhs
%
  function dy = SIR_rhs(t,y,params);

  global N 

  gamma = params(1);
  k = params(2);
  r = params(3);
  delta = params(4);

  S = y(1);
  I = y(2);
  R = y(3);

%
%
  dy = [delta*(N-S) - gamma*k*I*S;
        gamma*k*I*S - (r+delta)*I;
        r*I - delta*R];

%
  end