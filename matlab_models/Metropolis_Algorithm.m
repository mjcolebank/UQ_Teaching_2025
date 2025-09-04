% Metropolis_Algorithm

function chain = Metropolis_Algorithm(f,data,prior_F,theta0,noise_var,M,covar,UB,LB)

n_par = length(theta0);
chain = zeros(n_par,M);
n_y   = length(data);

if isempty(covar)
    covar = eye(n_par).*0.1.*(abs(theta0));
end
lik_func = @(param,data) exp(-sum((f(param)-data).^2) ./ (2.*noise_var)).* (2.*pi.*noise_var).^(-n_y/2);
chain(:,1) = theta0;
lik_old = lik_func(theta0,data);
prior_old = prior_F(theta0);
R = chol(covar);
num_acc = 0;

    prob_acc = unifrnd(0,1,M,1);
for i=2:M
    theta_star = chain(:,i-1)+R'*normrnd(0,1,n_par,1);

    if sum(any(theta_star(:)>UB(:)))>0 || sum(any(theta_star(:)<LB(:)))>0
        lik_star = -inf;
    else
        lik_star = lik_func(theta_star,data);
    end
    prior_star = prior_F(theta_star);

    acc_prob = prior_star.*lik_star./(prior_old.*lik_old);
    if acc_prob>prob_acc(i)
        chain(:,i) = theta_star;
        lik_old = lik_star;
        prior_old = prior_star;
        num_acc = num_acc+1;
    else
        chain(:,i) = chain(:,i-1);
    end
end
disp('Acceptance rate')
disp((num_acc./M) .* 100);
end