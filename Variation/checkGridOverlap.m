function [isValid, badCells] = checkGridOverlap(x, y)
%CHECKGRIDOVERLAP   Verify a 2D structured grid has no overlapping or inverted cells (vectorized).
%   [isValid, badCells] = checkGridOverlap(x, y)
%   Inputs:
%     x, y   - (ni x nj) arrays of node coordinates.
%   Outputs:
%     isValid  - true if all cells have consistent orientation
%     badCells - (k x 2) array of [i,j] indices of invalid cells

[ni, nj] = size(x);
if ~isequal(size(y), [ni, nj])
    error('x and y must have the same dimensions');
end

% Compute local edge vectors for all cells
v1x = x(2:ni, 1:nj-1) - x(1:ni-1, 1:nj-1);
v1y = y(2:ni, 1:nj-1) - y(1:ni-1, 1:nj-1);
v2x = x(1:ni-1, 2:nj) - x(1:ni-1, 1:nj-1);
v2y = y(1:ni-1, 2:nj) - y(1:ni-1, 1:nj-1);

% Compute Jacobian determinant for each cell
detJ = v1x .* v2y - v1y .* v2x;

% Establish expected sign from first cell (fallback to +1)
orientationSign = sign(detJ(1,1));
if orientationSign == 0
    orientationSign = 1;
end

% Identify invalid cells where sign flips
badMask = (sign(detJ) ~= orientationSign);

% Results
isValid = ~any(badMask, 'all');
badCells = [];
if ~isValid
    [I, J] = find(badMask);
    badCells = [I, J];
end

% Throw error if requested and invalid
if nargout < 2 && ~isValid
    error('Grid has invalid cells at indices: %s', mat2str(badCells));
end
end
