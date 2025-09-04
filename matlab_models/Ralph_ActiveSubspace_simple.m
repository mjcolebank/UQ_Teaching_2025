%
%  Compute the active subspace for the function
%
%    y = exp(0.7q_1 + 0.3q_2)
%

  clear all
%%
  c1 = 0.7;
  c2 = 0.3;

  M = 1e+3;

  C = zeros(2,2);
  for j=1:M
    q1 = rand(1,1);
    q2 = rand(1,1);
    grad_f = [c1*exp(c1*q1 + c2*q2); c2*exp(c1*q1 + c2*q2)];
    C = C + grad_f*grad_f';
  end
  C = (1/M)*C
%%
  q1 = rand(1,M);
  q2 = rand(1,M);
  G = 1/sqrt(M)*[c1*exp(c1*q1 + c2*q2); c2*exp(c1*q1 + c2*q2)];
%%
  C_anal = [c1/(4*c2) 1/4;1/4 c2/(4*c1)]*(exp(2*c1)-1)*(exp(2*c2)-1)

  [V1,D] = eig(C_anal);
  [U,S,V] = svd(G);
%%
  U1 = -U(:,1)';%[0.9191 0.3939];
  for j=1:100
    q1 = rand(1,1);
    q2 = rand(1,1);
    q = [q1; q2];
    xi(j) = U1*q;
    f(j) = exp(0.7*q1 + 0.3*q2);
  end

  for j=1:20
    q1 = rand(1,1);
    q2 = rand(1,1);
    q = [q1; q2];
    xi_test(j) = U1*q;
    f_test(j) = exp(0.7*q1 + 0.3*q2);
  end
  uppx = max(max(xi),max(xi_test));
  uppy = max(max(f),max(f_test));

  lowx = min(min(xi),min(xi_test));
  lowy = min(min(f),min(f_test));
%%
  
%
% Fit a 3rd order polynomial
%

  p = polyfit(xi,f,3);
  xvals = linspace(lowx,uppx);
  fit = polyval(p,xvals);
  fit_test = polyval(p,xi_test);

  err = max(abs(f_test - fit_test))
%%
%
% plot the curves
%

  figure(1); clf;
  subplot(1,2,1); hold on;
  plot(xvals,fit,'-k','LineWidth',3)
  % xlim([0 limx]);
  % ylim([1 limy])
  hold on
  plot(xi,f,'rx','MarkerSize',8,'LineWidth',2)
  plot(xi_test,f_test,'bo','MarkerSize',8,'LineWidth',2)
  hold off
  set(gca,'Fontsize',[22]);
  legend('Response Surface','Training Points','Test Points','Location','NorthWest')
  xlabel('\zeta')
  ylabel('y') 
  grid on;

  subplot(1,2,2);hold on;
  plot(xi,f-polyval(p,xi),'rx');
  plot(xi_test,f_test-polyval(p,xi_test),'bo');
    legend('Training Points','Test Points','Location','NorthWest')
  set(gca,'Fontsize',[22]);
grid on;
%  print -depsc active_subspace_ex1 
  
