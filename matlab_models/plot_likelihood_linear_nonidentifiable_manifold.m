%
%   Plot the liklihood contours and manifold h(\theta_{sub}) = 4 - \theta-1 \theta_2 = 0 
%   for the structurally nonidentifiable parameters in Example 2.6.
%
% ORIGINALLY WRITTEN BY R. SMITH

  clear all
  close all

%

  dt1 = 0.1;
  dt2 = 0.1;

  theta1 = 1:dt1:3;
  theta2 = 1:dt2:3;
  nt1 = length(theta1);
  nt2 = length(theta2);
  theta1_val = 2;
  theta2_val = 2;
  [Theta1,Theta2] = meshgrid(theta1',theta2'); 
  T1 = 1:dt1:4;
  T1 = 1:dt1:3;

%  t = 0:1:10;
  t = 0:.1:1;
  nt = length(t);
  sigma = 0.1;
  sigma2 = sigma^2;
  ve = sigma*randn(1,nt);

  A = theta1'*theta2;

  ss1 = zeros(size(A));
  for j = 1:nt
    y = ones(nt1,nt2)*(theta1_val*theta2_val*t(j) + ve(j));
    ss1 = ss1 + (y - A*t(j)).^2; 
  end
  likelihood = exp(-(1/2*sigma2)*ss1);

  T2 = theta1_val*theta2_val./T1;

%
%  Plot contours 
%

  c = gray;
%  c = c(1:40,:);
  c = c(1:50,:);
  c = flipud(c);
  figure(1)
  contour(Theta1,Theta2,likelihood,10,'linewidth',2)  
  %colormap(winter)
  colormap(c);
  axis('equal') 
  axis([1 3 1 3])
  hold on
  plot(T1,T2,'--b','linewidth',3)
  hold off
  set(gca,'Fontsize',[20]);
  xlabel('\theta_1')
  ylabel('\theta_2')

  figure(2)
  plot(T1,T2,'-k','linewidth',3)
  axis('equal')
  axis([1 3 1 3])
  set(gca,'Fontsize',[20]);
  xlabel('\theta_1')
  ylabel('\theta_2')