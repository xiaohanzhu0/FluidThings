clear
Nx1 = 10;
Nx2 = 10;
N = Nx1*Nx2;
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);

[x1, x2] = meshgrid(s1, s2);
[x1_exact, x2_exact]  = meshgrid(s1, s2);

M = zeros(2, 2, Nx2, Nx1);
M(1,1,:,:) = 40000 * (1 + 10*15 * x1).^(-2);
M(2,2,:,:) = 40000 * (1 + 10*15 * x2).^(-2);

M = reshape(M, [2,2,N]);
X = [x1(:),x2(:)];
M_tilde = quadratic_gradation_2d(X, M, 0.1);
max(abs(M_tilde-M),[],'all')


function new_metrics = quadratic_gradation_2d(vertices, metrics, h_gradation)
    % Performs anisotropic mesh gradation using the fully quadratic algorithm.
    %
    % Args:
    %     vertices (matrix): An (N, 2) matrix of vertex coordinates.
    %     metrics (3D array): An (2, 2, N) array of metric tensors.
    %     h_gradation (double): The gradation factor (typically > 1.0).
    %
    % Returns:
    %     new_metrics (3D array): An (2, 2, N) array of the new, graded metrics.
    
    num_vertices = size(vertices, 1);
    new_metrics = metrics; % Pre-allocate with the initial metrics
    
    fprintf('Starting quadratic gradation for %d vertices...\n', num_vertices);

    % For each vertex v_i, compute its new graded metric.
    for i = 1:num_vertices
        if mod(i, 10) == 0
            fprintf('  Processing vertex %d/%d\n', i, num_vertices);
        end
        
        current_metric = metrics(:, :, i);
        
        % Intersect with the grown metric from every other vertex v_j.
        for j = 1:num_vertices
            % A vertex does not influence itself.
            if i == j
                continue;
            end

            vec_ij = vertices(i, :)' - vertices(j, :)'; % Column vector
            metric_j = metrics(:, :, j);

            % Calculate the length of the vector from j to i in j's metric.
            length_in_j = calculate_length(vec_ij, metric_j);

            % The size of the grown metric is h * l.
            grown_size_at_i = h_gradation * length_in_j;
            
            if grown_size_at_i < 1e-9 % Avoid division by zero
                continue;
            end

            % The grown metric is isotropic, with size `grown_size_at_i`.
            % M_grad = (1 / size^2) * Identity
            grown_metric_val = 1.0 / (grown_size_at_i ^ 2);
            grown_metric = eye(2) * grown_metric_val;
            
            % Intersect the current metric with the grown one.
            current_metric = intersect_metrics(current_metric, grown_metric);
        end
        new_metrics(:, :, i) = current_metric;
    end
    
    fprintf('Gradation complete.\n');
end


% =========================================================================
% Local Helper Functions
% =========================================================================

function m_intersect = intersect_metrics(m1, m2)
    % Computes the intersection of two metric tensors, m1 and m2.
    if rcond(m1) < eps
        warning('M1 is singular or badly conditioned, returning M2.');
        m_intersect = m2;
        return;
    end
    
    m1_inv = inv(m1);
    a_star = m1_inv * m2;

    % The basis P is formed by the eigenvectors of A*
    [p_matrix, ~] = eig(a_star);
    
    if rcond(p_matrix) < eps
       warning('Eigenvector matrix is singular, returning original metric.');
       m_intersect = m1;
       return;
    end
    p_inv = inv(p_matrix);

    % The metrics M1 and M2 are diagonal in the basis P.
    d1_diag = zeros(2, 1);
    d2_diag = zeros(2, 1);
    for k = 1:2
        d1_diag(k) = p_matrix(:, k)' * m1 * p_matrix(:, k);
        d2_diag(k) = p_matrix(:, k)' * m2 * p_matrix(:, k);
    end

    % The diagonal of the intersected metric is the maximum of the two.
    d_intersect_diag = max(d1_diag, d2_diag);
    d_intersect = diag(d_intersect_diag);

    % Transform the intersected metric back to the original basis.
    m_intersect = p_inv' * d_intersect * p_inv;
end

function len = calculate_length(vec, metric)
    % Calculates the length of a vector in the space defined by a metric.
    if all(vec == 0)
        len = 0.0;
        return;
    end
    len = sqrt(vec' * metric * vec);
end

function plot_metrics(ax, vertices, metrics, plot_title, color)
    % Plots a metric field as a collection of ellipses.
    title(ax, plot_title);
    hold(ax, 'on');
    
    for i = 1:size(vertices, 1)
        m = metrics(:, :, i);
        [eigvecs, eigvals_matrix] = eig(m);
        eigvals = diag(eigvals_matrix);
        
        % Ellipse radii are 1/sqrt(eigenvalues)
        width = 2 / sqrt(eigvals(1));
        height = 2 / sqrt(eigvals(2));
        angle = atan2(eigvecs(2, 1), eigvecs(1, 1));

        % Scale ellipses for better visualization
        scale_factor = 0.1;
        
        % MATLAB's rectangle function can draw ellipses
        pos = [vertices(i,1) - width*scale_factor/2, vertices(i,2) - height*scale_factor/2, width*scale_factor, height*scale_factor];
        
        % Create a rotated rectangle (ellipse)
        h = rectangle(ax, 'Position', pos, 'Curvature', [1 1], 'EdgeColor', color, 'LineWidth', 1);
        rotate(h, [0 0 1], rad2deg(angle), vertices(i,:));
    end
    
    axis(ax, 'equal');
    grid(ax, 'on');
    hold(ax, 'off');
end