% Proper Orthogonal Decomposition (POD)

% Generate or load the data matrix. Each column represents a snapshot.
% For demonstration, we'll create synthetic data.
n = 100; % number of spatial points
m = 50;  % number of snapshots
data = rand(n, m); % synthetic data

% Compute the mean of the snapshots and subtract it to center the data
mean_snapshot = mean(data, 2); % mean across snapshots
data_centered = data - mean_snapshot;

% Compute the covariance matrix
cov_matrix = data_centered' * data_centered / (m - 1);

% Perform Singular Value Decomposition (SVD)
[U, S, V] = svd(data_centered, 'econ');

% POD modes are given by the left singular vectors
POD_modes = U;

% Singular values from S represent the energy contribution of each mode
singular_values = diag(S);

% Energy contribution of each mode
energy_fraction = singular_values.^2 / sum(singular_values.^2);

% Plot the energy contribution of modes
figure;
plot(energy_fraction, 'o-');
xlabel('Mode Number');
ylabel('Fraction of Energy');
title('Energy Contribution of POD Modes');
grid on;

% Optional: Reconstruct data using the first few dominant modes
n_modes = 5; % number of modes to keep
reconstructed_data = mean_snapshot + POD_modes(:, 1:n_modes) * S(1:n_modes, 1:n_modes) * V(:, 1:n_modes)';

% Plot an example of original vs reconstructed snapshot
snapshot_idx = 1; % index of a snapshot to compare
figure;
plot(data(:, snapshot_idx), 'r', 'DisplayName', 'Original Snapshot');
hold on;
plot(reconstructed_data(:, snapshot_idx), 'b--', 'DisplayName', 'Reconstructed Snapshot');
xlabel('Spatial Point Index');
ylabel('Amplitude');
title(['Snapshot Comparison for Mode Number ', num2str(n_modes)]);
legend;
grid on;