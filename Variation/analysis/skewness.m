function [Theta, Theta_1, Theta_2, Theta_inf] = skewness(x1, x2)
    % --- Input Validation ---
    if ~isequal(size(x1), size(x2))
        error('Input matrices X and Y must have the same dimensions.');
    end

    [Nx1, Nx2] = size(x1);
    Theta = compute_orthogonality_field(x1, x2);
    Theta_1 = norm(Theta(:), 1) / (Nx1*Nx2);
    Theta_2 = norm(Theta(:), 2) / sqrt(Nx1*Nx2);
    Theta_inf = norm(Theta(:), Inf);

    pcolor(x1,x2,Theta); colorbar
    title("Deviation from Orthogonality")
    xlabel(sprintf('$\\left\\|\\Theta\\right\\|_1=%.3g,\\; \\left\\|\\Theta\\right\\|_2=%.3g,\\; \\left\\|\\Theta\\right\\|_{\\infty}=%.3g$', ...
               Theta_1, Theta_2, Theta_inf), ...
               'Interpreter','latex');
end


function Theta = compute_orthogonality_field(X, Y)
% computeOrthogonality Calculates the deviation from orthogonality for a
% structured 2D mesh, including boundary points.
%
% Syntax:
%   Theta = computeOrthogonality(X, Y)
%
% Description:
%   This function computes the deviation from orthogonality for a structured
%   grid defined by the coordinate matrices X and Y. The deviation, Theta,
%   is calculated at each node as Theta = |90 - Gamma|, where Gamma is the
%   angle in degrees between the two intersecting grid lines at that node.
%
%   This version computes the deviation for all nodes:
%   - Interior nodes use a central difference scheme for higher accuracy.
%   - Boundary nodes (edges and corners) use a one-sided (forward or
%     backward) difference scheme.
%
% Inputs:
%   X - A 2D array representing the x-coordinates of the grid nodes.
%   Y - A 2D array representing the y-coordinates of the grid nodes.
%       X and Y must have the same dimensions.
%
% Outputs:
%   Theta - A 2D array of the same size as X and Y, containing the
%           orthogonality deviation in degrees for every node.
%
%   % Expected output will show a non-zero Theta at the perturbed corner.


[num_j, num_i] = size(X);

% --- Initialization ---
% Initialize the output matrix. All elements will be computed.
Theta = zeros(num_j, num_i);

% --- Helper function for angle calculation ---
% This anonymous function avoids repeating the angle formula.
% It takes the x and y components of two vectors and returns the
% deviation from 90 degrees.
calculate_angle = @(v1x, v1y, v2x, v2y) ...
    abs(90 - rad2deg(acos((v1x.*v2x + v1y.*v2y) ./ (sqrt(v1x.^2 + v1y.^2) .* sqrt(v2x.^2 + v2y.^2) + eps))));

% --- Interior Nodes (using central differences) ---
% For a node (j,i), vector v1 is from (j,i-1) to (j,i+1) and v2 is from (j-1,i) to (j+1,i).
if num_j > 2 && num_i > 2
    v1_x_interior = X(2:end-1, 3:end)   - X(2:end-1, 1:end-2);
    v1_y_interior = Y(2:end-1, 3:end)   - Y(2:end-1, 1:end-2);
    v2_x_interior = X(3:end,   2:end-1) - X(1:end-2, 2:end-1);
    v2_y_interior = Y(3:end,   2:end-1) - Y(1:end-2, 2:end-1);
    Theta(2:end-1, 2:end-1) = calculate_angle(v1_x_interior, v1_y_interior, v2_x_interior, v2_y_interior);
end

% --- Boundary Nodes ---

% -- Edges (excluding corners) --

% Top Edge (j=1, i=2:num_i-1)
if num_j > 1 && num_i > 2
    v1_x = X(1, 3:end)   - X(1, 1:end-2); % Central diff for i
    v1_y = Y(1, 3:end)   - Y(1, 1:end-2);
    v2_x = X(2, 2:end-1) - X(1, 2:end-1); % Forward diff for j
    v2_y = Y(2, 2:end-1) - Y(1, 2:end-1);
    Theta(1, 2:end-1) = calculate_angle(v1_x, v1_y, v2_x, v2_y);
end

% Bottom Edge (j=num_j, i=2:num_i-1)
if num_j > 1 && num_i > 2
    v1_x = X(num_j, 3:end)   - X(num_j, 1:end-2); % Central diff for i
    v1_y = Y(num_j, 3:end)   - Y(num_j, 1:end-2);
    v2_x = X(num_j, 2:end-1) - X(num_j-1, 2:end-1); % Backward diff for j
    v2_y = Y(num_j, 2:end-1) - Y(num_j-1, 2:end-1);
    Theta(num_j, 2:end-1) = calculate_angle(v1_x, v1_y, v2_x, v2_y);
end

% Left Edge (j=2:num_j-1, i=1)
if num_j > 2 && num_i > 1
    v1_x = X(2:end-1, 2) - X(2:end-1, 1); % Forward diff for i
    v1_y = Y(2:end-1, 2) - Y(2:end-1, 1);
    v2_x = X(3:end, 1)   - X(1:end-2, 1); % Central diff for j
    v2_y = Y(3:end, 1)   - Y(1:end-2, 1);
    Theta(2:end-1, 1) = calculate_angle(v1_x, v1_y, v2_x, v2_y);
end

% Right Edge (j=2:num_j-1, i=num_i)
if num_j > 2 && num_i > 1
    v1_x = X(2:end-1, num_i) - X(2:end-1, num_i-1); % Backward diff for i
    v1_y = Y(2:end-1, num_i) - Y(2:end-1, num_i-1);
    v2_x = X(3:end, num_i)   - X(1:end-2, num_i);   % Central diff for j
    v2_y = Y(3:end, num_i)   - Y(1:end-2, num_i);
    Theta(2:end-1, num_i) = calculate_angle(v1_x, v1_y, v2_x, v2_y);
end

% -- Corners --
if num_j > 1 && num_i > 1
    % Top-Left Corner (j=1, i=1)
    v1_x = X(1, 2) - X(1, 1); % Forward diff i
    v1_y = Y(1, 2) - Y(1, 1);
    v2_x = X(2, 1) - X(1, 1); % Forward diff j
    v2_y = Y(2, 1) - Y(1, 1);
    Theta(1, 1) = calculate_angle(v1_x, v1_y, v2_x, v2_y);

    % Top-Right Corner (j=1, i=num_i)
    v1_x = X(1, num_i) - X(1, num_i-1); % Backward diff i
    v1_y = Y(1, num_i) - Y(1, num_i-1);
    v2_x = X(2, num_i) - X(1, num_i);   % Forward diff j
    v2_y = Y(2, num_i) - Y(1, num_i);
    Theta(1, num_i) = calculate_angle(v1_x, v1_y, v2_x, v2_y);

    % Bottom-Left Corner (j=num_j, i=1)
    v1_x = X(num_j, 2) - X(num_j, 1);   % Forward diff i
    v1_y = Y(num_j, 2) - Y(num_j, 1);
    v2_x = X(num_j, 1) - X(num_j-1, 1); % Backward diff j
    v2_y = Y(num_j, 1) - Y(num_j-1, 1);
    Theta(num_j, 1) = calculate_angle(v1_x, v1_y, v2_x, v2_y);

    % Bottom-Right Corner (j=num_j, i=num_i)
    v1_x = X(num_j, num_i) - X(num_j, num_i-1); % Backward diff i
    v1_y = Y(num_j, num_i) - Y(num_j, num_i-1);
    v2_x = X(num_j, num_i) - X(num_j-1, num_i); % Backward diff j
    v2_y = Y(num_j, num_i) - Y(num_j-1, num_i);
    Theta(num_j, num_i) = calculate_angle(v1_x, v1_y, v2_x, v2_y);
end

end

