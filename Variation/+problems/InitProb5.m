function [T1, T2, M_samp, Mfun] = InitProb5(cf)
    if cf.new_airfoil
        [data, ~, ~, ~, ~] = readTurtleFields(cf.metric_datapath);
        [x_cv, ~, ~, ~, ~, ~, ~, ~, ~, ~] = readTurtleGrid(cf.airfoil_datapath);
        M_samp.x_metric = [x_cv{1,9}(:,:,1,1), x_cv{1,8}(:,:,1,1);
        x_cv{1,10}(:,:,1,1),x_cv{1,1}(:,:,1,1);
        x_cv{1,11}(:,:,1,1),x_cv{1,2}(:,:,1,1);
        x_cv{1,12}(:,:,1,1),x_cv{1,3}(:,:,1,1);
        x_cv{1,13}(:,:,1,1),x_cv{1,4}(:,:,1,1);
        x_cv{1,14}(:,:,1,1),x_cv{1,5}(:,:,1,1);
        x_cv{1,15}(:,:,1,1),x_cv{1,6}(:,:,1,1);
        x_cv{1,16}(:,:,1,1),x_cv{1,7}(:,:,1,1)];

        M_samp.y_metric = [x_cv{1,9}(:,:,1,2), x_cv{1,8}(:,:,1,2);
        x_cv{1,10}(:,:,1,2),x_cv{1,1}(:,:,1,2);
        x_cv{1,11}(:,:,1,2),x_cv{1,2}(:,:,1,2);
        x_cv{1,12}(:,:,1,2),x_cv{1,3}(:,:,1,2);
        x_cv{1,13}(:,:,1,2),x_cv{1,4}(:,:,1,2);
        x_cv{1,14}(:,:,1,2),x_cv{1,5}(:,:,1,2);
        x_cv{1,15}(:,:,1,2),x_cv{1,6}(:,:,1,2);
        x_cv{1,16}(:,:,1,2),x_cv{1,7}(:,:,1,2)];

        aux1 = cat(1,data{1,1}{1,9}(:,:,1,1:2,1:2), data{1,1}{1,10}(:,:,1,1:2,1:2), ...
            data{1,1}{1,11}(:,:,1,1:2,1:2), data{1,1}{1,12}(:,:,1,1:2,1:2), ...
            data{1,1}{1,13}(:,:,1,1:2,1:2), data{1,1}{1,14}(:,:,1,1:2,1:2), ...
            data{1,1}{1,15}(:,:,1,1:2,1:2), data{1,1}{1,16}(:,:,1,1:2,1:2));
        aux2 = cat(1,data{1,1}{1,8}(:,:,1,1:2,1:2), data{1,1}{1,1}(:,:,1,1:2,1:2), ...
            data{1,1}{1,2}(:,:,1,1:2,1:2), data{1,1}{1,3}(:,:,1,1:2,1:2), ...
            data{1,1}{1,4}(:,:,1,1:2,1:2), data{1,1}{1,5}(:,:,1,1:2,1:2), ...
            data{1,1}{1,6}(:,:,1,1:2,1:2), data{1,1}{1,7}(:,:,1,1:2,1:2));
        aux = cat(2,aux1,aux2);
        dim = size(aux);
        M = reshape(aux, dim(1), dim(2), dim(4), dim(5));
        M_samp.metric(:,:,1) = M(:,:,1,1);
        M_samp.metric(:,:,2) = M(:,:,1,2);
        M_samp.metric(:,:,3) = M(:,:,2,2);
    else
        addpath('~/Files/data/Mesh_generation/Airfoil');
        % Processing sampled metric tensor point
        M_samp = load('metric_airfoil_bodyfitted.mat');
    end

    h = gaussian1d(50, 100);
    
    for i=1:5
        M_samp.metric(:,:,1) = conv2(h,h,M_samp.metric(:,:,1),'same');
        M_samp.metric(:,:,2) = conv2(h,h,M_samp.metric(:,:,2),'same');
        M_samp.metric(:,:,3) = conv2(h,h,M_samp.metric(:,:,3),'same');
    end


    M_samp.F11 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
    M_samp.F12 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
    M_samp.F22 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));
    
    Mfun = @(x1,x2) Prob5Metric(x1, x2, M_samp);
    function M = Prob5Metric(x1,x2,M_samp)
        M.M11 = M_samp.F11(x1,x2);
        M.M12 = M_samp.F12(x1,x2);
        M.M22 = M_samp.F22(x1,x2);
    end

    gd = Hyperbolic(cf.Nx1, cf.Nx2, cf.alpha, cf.append_trail);
    T1 = gd.x';
    T2 = gd.y';

end

%% Apply TFI
function [X, Y] = transfiniteInterp(b1,b2,b3,b4, bi)
%TRANSFINITEINTERP  Structured grid via transfinite (Coons) interpolation
%  [X,Y] = TRANSFINITEINTERP(b1,b2,b3,b4)
%  [X,Y] = TRANSFINITEINTERP(b1,b2,b3,b4,bi)
%
%  b1: 2×Nx1 bottom (ξ∈[0,1])
%  b2: 2×Nx2 right  (η∈[0,1])
%  b3: 2×Nx1 top    (ξ∈[0,1])
%  b4: 2×Nx2 left   (η∈[0,1])
%  bi: 2×Nx2 interior layer (optional; η∈[0,1])
%
%  If bi is provided, we split at ξ = 0.5 (mid‐column) into:
%    • Left patch:   (b1L, bi, b3L, b4)
%    • Right patch:  (b1R, b2, b3R, bi)
%  and then stitch them, dropping the duplicate seam column.

  narginchk(4,5);
  [~, Nx1] = size(b1);
  [~, Nx2] = size(b2);
  assert(all(size(b3)==[2 Nx1]), 'b3 must be 2×Nx1');
  assert(all(size(b4)==[2 Nx2]), 'b4 must be 2×Nx2');
  if nargin==5
    assert(all(size(bi)==[2 Nx2]), 'bi must be 2×Nx2');
  end

  %— BASE CASE: no interior curve → standard Coons patch
  if nargin<5
    %— parametric coords
    xi  = linspace(0,1,Nx1);
    eta = linspace(0,1,Nx2);
    [XI, ETA] = meshgrid(xi, eta);    % size Nx2×Nx1

    %— unpack
    b1x = b1(1,:);  b1y = b1(2,:);
    b3x = b3(1,:);  b3y = b3(2,:);
    b4x = b4(1,:)'; b4y = b4(2,:)';
    b2x = b2(1,:)'; b2y = b2(2,:)';

    %— corners
    Pbl = b1(:,1);   Pbr = b1(:,end);
    Ptl = b3(:,1);   Ptr = b3(:,end);

    %— blends
    P1x = (1-ETA).*b1x(ones(Nx2,1),:) + ETA.*b3x(ones(Nx2,1),:);
    P1y = (1-ETA).*b1y(ones(Nx2,1),:) + ETA.*b3y(ones(Nx2,1),:);
    P2x = (1-XI).*b4x + XI.*b2x;
    P2y = (1-XI).*b4y + XI.*b2y;

    %— bilinear corner correction
    Cx = (1-XI).*(1-ETA)*Pbl(1) + XI.*(1-ETA)*Pbr(1) ...
       + XI.*ETA*Ptr(1)     + (1-XI).*ETA*Ptl(1);
    Cy = (1-XI).*(1-ETA)*Pbl(2) + XI.*(1-ETA)*Pbr(2) ...
       + XI.*ETA*Ptr(2)     + (1-XI).*ETA*Ptl(2);

    %— assemble
    X = P1x + P2x - Cx;
    Y = P1y + P2y - Cy;

    return
  end

  %— RECURSIVE CASE: interior curve supplied → two patches
  midcol = ceil(Nx1/2);

  % split bottom/top into left/right chunks
  b1L = b1(:,     1:midcol);    b3L = b3(:,     1:midcol);
  b1R = b1(:, midcol: end);     b3R = b3(:, midcol: end);

  % left patch   uses (b1L,  bi, b3L, b4)
  [Xl, Yl] = transfiniteInterp(b1L, bi, b3L, b4);

  % right patch  uses (b1R, b2, b3R, bi)
  [Xr, Yr] = transfiniteInterp(b1R, b2, b3R, bi);

  % stitch seams (drop duplicate interior column at Xr(:,1)==bi)
  X = [ Xl(:,1:midcol),  Xr(:,2:end) ];
  Y = [ Yl(:,1:midcol),  Yr(:,2:end) ];
end


function h = gaussian1d(sigma, ksize)
%GAUSSIAN1D  Return a normalized 1-D Gaussian kernel.
%   h = GAUSSIAN1D(sigma, ksize) returns a row vector of length ksize
%   containing samples of a Gaussian with standard deviation sigma.
%
%   Example:
%     h = gaussian1d(2, 11);
%     plot(h, '-o'); axis tight; title('1D Gaussian Kernel');

    % ensure kernel size is a positive integer
    validateattributes(ksize, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(sigma, {'numeric'}, {'scalar','positive'});

    % define sample positions centered at zero
    half = (ksize-1)/2;
    if mod(ksize,2)==0
        % even length—center between the two middle samples
        x = -half:(half-1);
    else
        % odd length—symmetric around zero
        x = -floor(half):floor(half);
    end

    % compute the Gaussian
    h = exp( - (x.^2) / (2*sigma^2) );

    % normalize so sum(h) == 1
    h = h / sum(h);
end

function M_graded = grade_metric_field_iterative(M, h, max_iter, tol)
%GRADE_METRIC_FIELD_ITERATIVE Grades a metric field using an efficient iterative method.
%
%   M_graded = GRADE_METRIC_FIELD_ITERATIVE(M, h, max_iter, tol) takes a
%   N*M*2*2 metric field M, a gradation factor h, a maximum number of
%   iterations, and a convergence tolerance.
%
%   This function implements a local, iterative approach based on the
%   methods described in "Size gradation control of anisotropic meshes"
%   by F. Alauzet. Instead of a global comparison, each metric is updated
%   based on its immediate neighbors. This process is repeated until the
%   metric field converges, which is significantly faster.
%
%   Inputs:
%       M - A N*M*2*2 array, where M(i,j,:,:) is a 2x2 SPD metric tensor.
%       h - A scalar gradation factor, typically > 1.
%       max_iter - The maximum number of iterations to perform (e.g., 100).
%       tol - Convergence tolerance. The iteration stops when the maximum
%             change in any metric is below this value (e.g., 1e-6).
%
%   Output:
%       M_graded - A N*M*2*2 array representing the graded metric field.

% Input validation
if nargin < 4
    tol = 1e-6; % Default tolerance
end
if nargin < 3
    max_iter = 100; % Default max iterations
end
if ndims(M) ~= 4 || size(M, 3) ~= 2 || size(M, 4) ~= 2
    error('Input metric field M must be a N*M*2*2 array.');
end
if ~isscalar(h) || h <= 1
    error('Gradation factor h must be a scalar greater than 1.');
end

[N, M_dim, ~, ~] = size(M);

% Initialize the graded metric field with the original field
M_graded = M;

fprintf('Starting iterative gradation...\n');

% --- Main Iteration Loop ---
for iter = 1:max_iter
    max_change = 0;
    M_old = M_graded; % Store the state at the start of the iteration (Jacobi-style)

    % Loop over each point in the metric field
    for i = 1:N
        for j = 1:M_dim
            % Start with the original (ungraded) metric at this point.
            % This ensures the final metric is always a refinement of the original.
            M_current = squeeze(M(i, j, :, :));

            % --- Neighborhood Loop ---
            % Iterate over the 8-connected neighbors (and the point itself)
            for p = max(1, i-1):min(N, i+1)
                for q = max(1, M_dim-1):min(M_dim, j+1)
                    
                    % Get the neighbor's metric from the *previous* iteration
                    M_neighbor = squeeze(M_old(p, q, :, :));
                    
                    % Euclidean distance to the neighbor (1 for cardinal, sqrt(2) for diagonal)
                    dist = sqrt((i-p)^2 + (j-q)^2);
                    
                    % If dist is 0, it's the point itself. The intersection
                    % would be identity, so we can skip it.
                    if dist == 0
                        continue;
                    end
                    
                    % Propagate the neighbor's metric to the current point
                    scaling_factor = h^(2*dist);
                    M_modified = M_neighbor / scaling_factor;
                    
                    % The singularity check is less critical here since dist is small,
                    % but it's good practice to keep it.
                    if rcond(M_modified) < eps
                        continue;
                    end

                    % Intersect the current metric with the propagated neighbor metric
                    [V, D] = eig(M_current, M_modified);
                    D_intersect = max(eye(2), D);
                    M_intersect = inv(V') * D_intersect * inv(V);
                    
                    % Update the current metric with the result of the intersection
                    M_current = M_intersect;
                end
            end
            
            % After intersecting with all neighbors, update the graded field
            M_graded(i, j, :, :) = M_current;
            
            % Track the maximum change for the convergence check
            change = norm(squeeze(M_graded(i, j, :, :)) - squeeze(M_old(i, j, :, :)), 'fro');
            if change > max_change
                max_change = change;
            end
        end
    end
    
    fprintf('Iteration %d: Max change = %e\n', iter, max_change);
    
    % --- Convergence Check ---
    if max_change < tol
        fprintf('Convergence reached after %d iterations.\n', iter);
        break;
    end

    figure(30)
    f = pcolor(M_graded(:,:,1,1)); 
    set(f, 'EdgeColor', 'none'); caxis([1e3,5e4]);
end

if iter == max_iter && max_change >= tol
    fprintf('Warning: Maximum iterations reached without full convergence.\n');
end

end
