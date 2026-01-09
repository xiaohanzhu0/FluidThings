clear
Nx1 = 40;
Nx2 = 40;
N = Nx1*Nx2;
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);

[x1, x2] = meshgrid(s1, s2);
[x1_exact, x2_exact]  = meshgrid(s1, s2);

M = zeros(2, 2, Nx2, Nx1);
M(1,1,:,:) = 40000 * (1 + 15 * x1).^(-2);
M(2,2,:,:) = 40000 * (1 + 15 * x2).^(-2);


M_tilde = grade_structured_mesh(x1, x2, M, 2, 50, 1e-4);
%%
clear
% --- Example Usage ---
% Define a simple structured grid
[X, Y] = meshgrid(0:0.1:1, 0:0.1:1);
[rows, cols] = size(X);

% Define an initial metric (e.g., identity matrix everywhere)
initial_metrics = repmat(eye(2), [1, 1, rows, cols]);

% Make the metric anisotropic in a certain region (e.g., right half of the domain)
for i = 1:rows
    for j = 1:cols
        if X(i,j) > 0.5
            initial_metrics(:, :, i, j) = [10 0; 0 1]; % Stretch elements horizontally
        end
    end
end

% --- Set iterative parameters ---
gradation_coefficient = 1.1;
max_iterations = 50; % Set a limit to prevent infinite loops
tolerance = 1e-4;    % Stop when the change is very small

% Apply the grading algorithm
graded_metrics = grade_structured_mesh(X, Y, initial_metrics, gradation_coefficient, max_iterations, tolerance);

% The 'graded_metrics' variable now holds the final, converged metric field.
fprintf('Final grading complete. The output variable "graded_metrics" is a %s array.\n', class(graded_metrics));
disp(['Size: ', num2str(size(graded_metrics))]);



%%
function graded_metrics = grade_structured_mesh(X, Y, initial_metrics, gradation_coefficient, max_iterations, tolerance)
    % grade_structured_mesh: Applies an iterative anisotropic mesh grading algorithm to a structured grid.
    %
    % Inputs:
    %   X, Y:                  MxN matrices of vertex x and y coordinates.
    %   initial_metrics:       2x2xMxN matrix of initial metric tensors at each vertex.
    %   gradation_coefficient: Scalar controlling the rate of size change.
    %   max_iterations:        Maximum number of iterations to perform.
    %   tolerance:             Convergence tolerance. The iteration stops when the change in metrics is below this value.
    %
    % Outputs:
    %   graded_metrics:        2x2xMxN matrix of the resulting graded metrics at each vertex.
    omega = 1;
    [rows, cols] = size(X);
    graded_metrics = initial_metrics;
    
    % Define the offsets for the 4-connectivity neighbors (up, down, left, right)
    di = [-1, 1, 0, 0];
    dj = [0, 0, -1, 1];

    fprintf('Starting iterative grading...\n');
    
    % --- Outer loop for iterative propagation ---
    for iter = 1:max_iterations
        metrics_before_iter = graded_metrics;
        total_change = 0;

        % Iterate through each vertex in the grid
        for i = 1:rows
            for j = 1:cols
                M_old = graded_metrics(:, :, i, j);
                v_current_coords = [X(i, j), Y(i, j)];
                M_current = graded_metrics(:, :, i, j);

                % Iterate through the neighbors of the current vertex
                for k = 1:4
                    ni = i + di(k); % Neighbor i-index
                    nj = j + dj(k); % Neighbor j-index

                    % Check if the neighbor is within the grid boundaries
                    if ni >= 1 && ni <= rows && nj >= 1 && nj <= cols
                        v_neighbor_coords = [X(ni, nj), Y(ni, nj)];
                        % The key is to compare against the metric field from the *previous* full iteration
                        M_neighbor = metrics_before_iter(:, :, ni, nj);

                        % Propagate the neighbor's metric to the current vertex
                        dist = norm(v_current_coords - v_neighbor_coords);
                        propagated_M = M_neighbor / (1 + gradation_coefficient * dist)^2;

                        % Intersect the current metric with the propagated one
                        M_current = metric_intersection(M_current, propagated_M);
                    end
                end
                graded_metrics(:, :, i, j) = (1-omega)*M_old + omega*M_current;
            end
        end

        % --- Check for convergence ---
        % Calculate the total change using the Frobenius norm of the difference
        change_matrix = graded_metrics - metrics_before_iter;
        for i = 1:rows
            for j = 1:cols
                total_change = total_change + norm(change_matrix(:,:,i,j), 'fro');
            end
        end
        
        fprintf('Iteration %d: Total Change = %e\n', iter, total_change);

        if total_change < tolerance
            fprintf('Convergence reached after %d iterations.\n', iter);
            break;
        end
        
        if iter == max_iterations
            fprintf('Maximum number of iterations reached.\n');
        end
    end
end

function M_intersect = metric_intersection(M1, M2)
    % metric_intersection: Computes the intersection of two metric tensors.
    % This finds the smallest ellipsoid that contains the two input ellipsoids.
    %
    % Inputs:
    %   M1, M2: 2x2 symmetric positive definite matrices (metric tensors).
    %
    % Outputs:
    %   M_intersect: The resulting intersected metric.

    % The intersection is found by simultaneous reduction.
    % We solve the generalized eigenvalue problem M2*v = lambda*M1*v
    A_star = M1 \ M2;
    [P, ~] = eig(A_star);

    % In the basis P, both matrices become diagonal
    lambda = diag(P' * M1 * P);
    mu = diag(P' * M2 * P);

    % The intersection corresponds to taking the max of the eigenvalues
    D_intersect = diag(min(lambda, mu));
    
    % Transform back to the original basis
    % Note: for symmetric M1, inv(P') is equivalent to P
    M_intersect = inv(P') * D_intersect * inv(P);
end

function M_interp = log_euclidean_interpolation(M1, M2, t)
    % log_euclidean_interpolation: Interpolates between two metric tensors.
    % This method ensures the interpolated metric remains symmetric positive definite.
    %
    % Inputs:
    %   M1, M2: 2x2 symmetric positive definite matrices.
    %   t:      Interpolation parameter (0 <= t <= 1).
    %
    % Outputs:
    %   M_interp: The interpolated metric.

    log_M1 = logm(M1);
    log_M2 = logm(M2);
    log_M_interp = (1 - t) * log_M1 + t * log_M2;
    M_interp = expm(log_M_interp);
end

function len = calculate_edge_length(v1_coords, v2_coords, M1, M2)
    % calculate_edge_length: Computes the length of an edge in the Riemannian space defined by the metrics.
    %
    % Inputs:
    %   v1_coords, v2_coords: 1x2 vectors of the two vertex coordinates.
    %   M1, M2:               Metric tensors (2x2) at the respective vertices.
    %
    % Outputs:
    %   len: The length of the edge.

    % Using a 2-point Gaussian quadrature for accurate integration along the edge
    gauss_points = [0.211324865405187; 0.788675134594813];
    weights = [0.5; 0.5];
    len = 0;
    edge_vec = v2_coords - v1_coords;

    for i = 1:length(gauss_points)
        t = gauss_points(i);
        % Interpolate the metric at the Gaussian point
        M_interp = log_euclidean_interpolation(M1, M2, t);
        % Calculate the length element at that point
        integrand = sqrt(edge_vec * M_interp * edge_vec');
        len = len + weights(i) * integrand;
    end
    % The length is the integral over the path, which is approximated by the quadrature
    len = len * norm(edge_vec);
end
