% Generate synthetic data from a sine function
x = linspace(0, 10, 100)'; % Input data points
y = sin(x) + 0.1 * randn(size(x)); % Noisy sine function output

% Divide data into training and test sets
train_idx = randperm(length(x), 10); % Select 80 points for training
test_idx = setdiff(1:length(x), train_idx); % Remaining points for testing
x_train = x(train_idx);
y_train = y(train_idx);
x_test = x(test_idx);
y_test = y(test_idx);



% Hyperparameters for Gaussian Process
length_scale = 2.0;
variance = 1.0;
noise_variance = 0.1;

% Compute covariance matrices
K_train = squared_exponential_kernel(x_train, x_train, length_scale, variance) + ...
          noise_variance * eye(length(x_train)); % Add noise variance
K_train_test = squared_exponential_kernel(x_train, x_test, length_scale, variance);
K_test = squared_exponential_kernel(x_test, x_test, length_scale, variance);

% GP Prediction
L = chol(K_train, 'lower'); % Cholesky decomposition
alpha = L' \ (L \ y_train); % Solve for alpha
y_pred = K_train_test' * alpha;

% Compute predictive variance
v = L \ K_train_test;
y_pred_var = diag(K_test) - sum(v.^2, 1)';

% Visualize results
figure;
hold on;
% plot(x, y, 'k.', 'DisplayName', 'Original data');
plot(x_train, y_train, 'b*', 'DisplayName', 'Training data');
plot(x_test, y_test, 'r*', 'DisplayName', 'Testing data');
plot(x_test, y_pred, 'g-', 'DisplayName', 'GP Predictions');
fill([x_test; flipud(x_test)], ...
    [y_pred - 2 * sqrt(y_pred_var); flipud(y_pred + 2 * sqrt(y_pred_var))], ...
    'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Prediction Uncertainty');
legend('show');
xlabel('x');
ylabel('y');
title('Gaussian Process Regression on Sine Function');
hold off;

% Define squared_exponential_kernel as a nested function above

%%
% Kernel function (squared exponential)
function K = squared_exponential_kernel(x1, x2, length_scale, variance)
    K = variance * exp(-(pdist2(x1, x2).^2) / (2 * length_scale^2));
end