%
%   Plot the likelihood contours for the identifiable parameters in Example 2.4. 
%
% ORIGINALLY WRITTEN BY R. SMITH

  clear all
  close all

%

  dt1 = 0.1;
  dt2 = 0.1;

  theta1 = 0:dt1:4;
  theta2 = 0:dt2:4;
  nt1 = length(theta1);
  nt2 = length(theta2);
  theta1_val = 2.0;
  theta2_val = 2.0;
  [Theta1,Theta2] = meshgrid(theta1',theta2'); 

  t = 0:0.1:1;
  nt = length(t);
  sigma = 0.1;
  sigma2 = sigma^2;
  ve = sigma*randn(1,nt);

  ss1 = zeros(size(Theta1));
  for j = 1:nt
    y = ones(nt1,nt2)*(theta1_val*t(j) + theta2_val + ve(j));
    ss1 = ss1 + (y - (Theta1*t(j) + Theta2)).^2; 
  end
  likelihood = exp(-(1/2*sigma2)*ss1);
  maxval = max(max(likelihood));
 [xval,yval] = find(likelihood==maxval);
%
%  Plot contours 
%

  c = gray;
%  c = c(1:40,:);
  c = c(1:50,:);
  c = flipud(c);
  figure(1)
  contour(Theta1,Theta2,likelihood,15,'linewidth',2)  
  colormap(c)
  axis('equal') 
  hold on
  plot(theta1(yval),theta2(xval),'bx','linewidth',8)
  hold off 
  set(gca,'Fontsize',[24]);
  xlabel('\theta_1')
  ylabel('\theta_2')
  print -depsc likelihood_a

