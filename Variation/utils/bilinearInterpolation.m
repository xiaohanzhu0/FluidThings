function M = bilinearInterpolation(M_sampled, x1_sampled, x2_sampled, x1, x2)
%BILINEARINTERPOLATION Interpolate a 2D array using bilinear interpolation.
%
%   M = bilinearInterpolation(M_sampled, x1_sampled, x2_sampled, x1, x2)
%
%   Inputs:
%       M_sampled   - 2D array of sampled data values.
%       x1_sampled  - 2D array of x-coordinates for each element in M_sampled.
%       x2_sampled  - 2D array of y-coordinates for each element in M_sampled.
%       x1          - Target x-coordinates for interpolation (vector or 2D array).
%       x2          - Target y-coordinates for interpolation (vector or 2D array).
%
%   Output:
%       M           - 2D array of interpolated values at the points defined by (x1, x2).
%
%   This function uses interp2 with the 'linear' method (i.e. bilinear interpolation).
%   If x1 and x2 are provided as vectors, a grid is created using meshgrid.

    % Check if target grid points are provided as vectors.
    if isvector(x1) && isvector(x2)
        [X1_target, X2_target] = meshgrid(x1, x2);
    else
        % Assume x1 and x2 are already 2D arrays.
        X1_target = x1;
        X2_target = x2;
    end

    % Perform bilinear interpolation using interp2.
    M = interp2(x1_sampled, x2_sampled, M_sampled, X1_target, X2_target, 'linear');
end
