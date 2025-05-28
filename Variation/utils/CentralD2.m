function d2 = CentralD2(f, direction, params)
% compute_second_derivative computes the second derivative of a 2D array.
%
%   d2 = compute_second_derivative(f, direction, grid_size) returns the 
%   second derivative of the 2D array f in the specified direction.
%   The 'direction' flag can be:
%       - 'horizontal' for the second derivative along columns (x-direction)
%       - 'vertical' for the second derivative along rows (y-direction)
%   The grid_size is the spacing between grid points.
%
%   The derivative is computed using the central difference formula:
%       f''(x) ≈ (f(x+grid_size) - 2*f(x) + f(x-grid_size)) / grid_size^2
%
%
% Example:
%   A = magic(5);
%   d2x = compute_second_derivative(A, 'horizontal', 1);
%   d2y = compute_second_derivative(A, 'vertical', 1);

    [nrows, ncols] = size(f);
    d2 = zeros(nrows, ncols);
    
    % Ensure the input array has sufficient size.
    if nrows < 3 || ncols < 3
        error('Input array must have at least 3 rows and 3 columns.');
    end

    switch lower(direction)
        case 'horizontal'
            % Compute horizontal second derivative:
            % Operate on the interior rows (2:end-1) and compute along columns
            d2(2:end-1,2:end-1) = ( f(2:end-1, 3:end) - 2 * f(2:end-1, 2:end-1) + f(2:end-1, 1:end-2) ) / params.dxi^2;
            
        case 'vertical'
            % Compute vertical second derivative:
            % Operate on the interior columns (2:end-1) and compute along rows
            d2(2:end-1,2:end-1) = ( f(3:end, 2:end-1) - 2 * f(2:end-1, 2:end-1) + f(1:end-2, 2:end-1) ) / params.deta^2;

        case 'mixed'
            % Compute mixed second derivative:
            d2 = (f(3:end, 3:end)-f(3:end, 1:end-2)-f(1:end-2, 3:end)+f(1:end-2,1:end-2))...
                 / (4*params.dxi*params.deta);
            
        otherwise
            error('Invalid direction. Use ''horizontal'' or ''vertical''.');
    end
end