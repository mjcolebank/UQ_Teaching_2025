
  clear all

  global N

  example = input('Example = ');

  S0 = 900;
  I0 = 100;
  R0 = 0;
  N = 1000;

  tf = 20;
  dt = 0.05;
  t_data = 0:dt:tf;
  ind = length(t_data);
  Y0 = [S0; I0; R0];
  I_mu = 0.5;
  I_eta = 0.1;
  I_gamma = 1.2;

%  M = 1000;
  M = 5000;
  if example==1
    alpha = 5;
    beta = 9;
  else
    alpha = 0.2;
    beta = 15;
  end
  %A = int*rand(M,4);
  %B = int*rand(M,4);

  A(:,1) = I_eta*rand(M,1);
  A(:,2) = betarnd(alpha,beta,M,1);
  A(:,3) = I_gamma*rand(M,1);
  A(:,4) = I_mu*rand(M,1);
  B(:,1) = I_eta*rand(M,1);
  B(:,2) = betarnd(alpha,beta,M,1);
  B(:,3) = I_gamma*rand(M,1);
  B(:,4) = I_mu*rand(M,1);
  C1 = A;
  C2 = A;
  C3 = A;
  C4 = A;
  C1(:,1) = B(:,1);
  C2(:,2) = B(:,2);
  C3(:,3) = B(:,3);
  C4(:,4) = B(:,4);

  factor = 1e-3;
  %factor = 1;
  ode_options = odeset('RelTol',1e-6);

  tic
  wait = waitbar(0,'Please wait...');
  for j=1:M
    params = A(j,:);
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    y_A(j,1) = factor*sum(dt*Y(:,3));
    y_A_end(j,1) = Y(ind,3);

    if example==1
      params(2) = 0.4;
    else
      params(1) = 0.05;
      params(3) = 0.6;
      params(4) = 0.25;
    end
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    y_A_end_fixed(j,1) = Y(ind,3);
    
    params = B(j,:);
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    y_B(j,1) = factor*sum(dt*Y(:,3));

    params = C1(j,:);
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    y_C1(j,1) = factor*sum(dt*Y(:,3));

    params = C2(j,:);
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    y_C2(j,1) = factor*sum(dt*Y(:,3));

    params = C3(j,:);
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    y_C3(j,1) = factor*sum(dt*Y(:,3));

    params = C4(j,:);
    [t,Y] = ode45(@SIR_model,t_data,Y0,ode_options,params,N);
    y_C4(j,1) = factor*sum(dt*Y(:,3));

    waitbar(j/M)
  end
  toc

  y_D = [y_A; y_B];
  f02 = ((1/(2*M))^2)*sum(y_D)*sum(y_D);
  S1 = ((1/M)*(y_B'*y_C1 - y_B'*y_A))/((1/(2*M))*y_D'*y_D - f02);
  S2 = ((1/M)*(y_B'*y_C2 - y_B'*y_A))/((1/(2*M))*y_D'*y_D - f02);
  S3 = ((1/M)*(y_B'*y_C3 - y_B'*y_A))/((1/(2*M))*y_D'*y_D - f02);
  S4 = ((1/M)*(y_B'*y_C4 - y_B'*y_A))/((1/(2*M))*y_D'*y_D - f02);

  ST1 = ((1/(2*M))*(y_A'*y_A - 2*y_A'*y_C1 + y_C1'*y_C1))/((1/(2*M))*y_D'*y_D - f02);
  ST2 = ((1/(2*M))*(y_A'*y_A - 2*y_A'*y_C2 + y_C2'*y_C2))/((1/(2*M))*y_D'*y_D - f02);
  ST3 = ((1/(2*M))*(y_A'*y_A - 2*y_A'*y_C3 + y_C3'*y_C3))/((1/(2*M))*y_D'*y_D - f02);
  ST4 = ((1/(2*M))*(y_A'*y_A - 2*y_A'*y_C4 + y_C4'*y_C4))/((1/(2*M))*y_D'*y_D - f02);

  S = [S1 S2 S3 S4]
  ST = [ST1 ST2 ST3 ST4]
%%
  x_dens = 0:10:1000;
%  [dens,x_i] = ksdensity(y_A_end,x_dens);
%  [dens_fixed,x_i] = ksdensity(y_A_end_fixed,x_dens);
  [~,dens,x_dens,~] = kde(y_A_end,2^4,min(y_A_end)-5,max(y_A_end)+5);
  [~,dens_fixed,x_dens,~] = kde(y_A_end_fixed,2^4,min(y_A_end)-5,max(y_A_end)+5);
  
%
%  Compute the energy statistics to determine if y_A_end and y_A_end_fixed are sampled
%  from the same distribution.
%

  n1 = 1000;
  n2 = 1000;
  W = [y_A_end; y_A_end_fixed];
  M_stat = 499;
  alpha1 = 0.05;
  alpha2 = 0.01;
  ind1 = (M_stat+1)*(1-alpha1);
  ind2 = (M_stat+1)*(1-alpha2);
  
  [test_stat] = test_statistic(n1,n2,y_A_end,y_A_end_fixed)

  for k=1:M_stat
    %A_stat = datasample(W,n1,'Replace',false);
    %B_stat = datasample(W,n2,'Replace',false);
    Wnew = datasample(W,n1+n2,'Replace',false);
    A_stat = Wnew(1:n1);
    B_stat = Wnew(n1+1:n1+n2);
    energy_stat(k) = test_statistic(n1,n2,A_stat,B_stat);
  end

  energy_stat_sort = sort(energy_stat);
  for k=1:M_stat
    if k==ind1
      crit1 = energy_stat_sort(k)
    elseif k==ind2
      crit2 = energy_stat_sort(k)
    end
  end


%
%  Plot the results
%

  Sval = Y(:,1);
  Ival = Y(:,2);
  Rval = Y(:,3); 

  figure(1)
  plot(t,Sval,'b',t,Ival,'k--',t,Rval,'r-.','linewidth',3)
  axis([0 20 0 1000])
  set(gca,'Fontsize',[22]);
  legend('Susceptible','Infected','Recovered','Location','Northeast')
  xlabel('Time (Days)') 
  ylabel('Number of Individuals')

  figure(2)
  plot(x_dens,dens,'k',x_dens,dens_fixed,'k--','linewidth',3)
%  axis([0 1000 0 0.012])
  axis([0 1000 0 2.2e-3]) 
  set(gca,'Fontsize',[22]);
  xlabel('Number of Recovered Individuals')
  ylabel('PDF')
  if example==1
    legend(' Random \eta, k, \gamma, \mu',' Random \eta, \gamma, \mu; Fixed k=0.4','Location','South')
%    print -depsc density5_9
  else
    legend(' Random \eta, k, \gamma, \mu',' Random k; Fixed \eta=0.05, \gamma=0.6, \mu=0.25','Location','Northeast')
   %print -depsc density2_15
  end

  x = 0:.01:1;
  y1 = betapdf(x,5,9);
  y2 = betapdf(x,0.2,15);
  figure(3)
  plot(x,y1,'k',x,y2,'--b','linewidth',3)
  set(gca,'Fontsize',[22]);
  xlabel('Normalized Number of Interactions')
  ylabel('PDF')
  legend('Beta(5,9)','Beta(0.2,15)','Location','Northeast')
%  print -depsc beta_plot


