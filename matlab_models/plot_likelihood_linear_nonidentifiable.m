%
%   Plot the likelihood contours and maximum likelihood value for the practically nonidentifiable
%   model in Example 2.5.
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
  theta1_val = 2.0;
  theta2_val = 2.0;
  [Theta1,Theta2] = meshgrid(theta1,theta2'); 

  %t = 0:1:10;
  eta = 1e-3;
  n = 11;
  %t = 1 - eta*(n-1)/2:eta:1 + eta*(n-1)/2;
  %t = -eta*(n-1)/2:eta:eta*(n-1)/2;
  %t = 1:eta:1+eta*(n-1); 
  t = 0:eta:eta*(n-1);
  nt = length(t);
  sigma = 0.1;
  sigma2 = sigma^2;
  ve = sigma*randn(1,nt);

  ss1 = zeros(size(Theta1));
  for j = 1:nt
    y = ones(nt1,nt2)*(theta1_val*t(j) + theta2_val + ve(j));
    %ss1 = ss1 + (y - (A*t(j) + ones(size(theta1'))*theta2)).^2; 
    ss1 = ss1 + (y - (Theta1*t(j) + Theta2)).^2;
  end
  likelihood = exp(-(1/2*sigma2)*ss1);
  maxval = max(max(likelihood));
 [xval,yval] = find(likelihood==maxval);

%
%  Plot contours 
%

  figure(1)
  c = gray;
  c = c(1:30,:);
  c = flipud(c);
  contour(Theta1,Theta2,likelihood,15,'linewidth',2)  
  colormap(c)
  axis('equal') 
  hold on
  plot(theta1(yval),theta2(xval),'bx','linewidth',8) 
  hold off 
  set(gca,'Fontsize',[24]);
  xlabel('\theta_1')
  ylabel('\theta_2')
  %print -depsc likelihood_ci
