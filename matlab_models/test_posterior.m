% Generate synthetic data
n = 10; % Number of observations
theta_true = [5; 3; 2]; % True parameter values
sigma_error = 2; % Standard deviation of measurement errors
x = linspace(0, 1, n)'; % Predictor variable x
z = linspace(1, 2, n)'; % Predictor variable z
y = theta_true(1) + theta_true(2)*x + theta_true(3)*z + sigma_error * randn(n, 1); % Observed data

% Define the prior: Multivariate Gaussian
mu_prior = [0; 0; 0]; % Prior mean for [theta_1, theta_2, theta_3]
sigma_prior = diag([5, 5, 5]); % Prior covariance matrix (diagonal for independence)

% Define grid points for theta_1, theta_2, and theta_3
theta_1 = linspace(-5, 15, 50);
theta_2 = linspace(-5, 10, 50);
theta_3 = linspace(-5, 10, 50);
[Theta1, Theta2, Theta3] = ndgrid(theta_1, theta_2, theta_3);

% Compute the prior
prior = mvnpdf([Theta1(:), Theta2(:), Theta3(:)], mu_prior', sigma_prior);
prior = reshape(prior, size(Theta1));

% Compute the likelihood
likelihood = ones(size(Theta1)); % Initialize likelihood
sigma_likelihood = sigma_error; % Measurement error standard deviation
for i = 1:n
    y_pred = Theta1 + Theta2 * x(i) + Theta3 * z(i); % Predicted value for each point on the grid
    likelihood = likelihood .* normpdf(y(i), y_pred, sigma_likelihood); % Update likelihood
end

% Compute the posterior (unnormalized)
posterior_unnormalized = prior .* likelihood;

% Normalize the posterior
posterior = posterior_unnormalized / trapz(theta_3, trapz(theta_2, trapz(theta_1, posterior_unnormalized, 1), 2));

% Plot marginal distributions of posterior
figure;
subplot(1, 3, 1);
contourf(theta_1, theta_2, squeeze(sum(posterior, 3)), 20, 'LineColor', 'none');
title('Posterior Marginal (\theta_1, \theta_2)');
xlabel('\theta_1');
ylabel('\theta_2');
colorbar;

subplot(1, 3, 2);
contourf(theta_1, theta_3, squeeze(sum(posterior, 2)), 20, 'LineColor', 'none');
title('Posterior Marginal (\theta_1, \theta_3)');
xlabel('\theta_1');
ylabel('\theta_3');
colorbar;

subplot(1, 3, 3);
contourf(theta_2, theta_3, squeeze(sum(posterior, 1)), 20, 'LineColor', 'none');
title('Posterior Marginal (\theta_2, \theta_3)');
xlabel('\theta_2');
ylabel('\theta_3');
colorbar;