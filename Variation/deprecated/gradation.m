clear
beta = 1.1;
nbr = [-1 0;   % left
        1 0;   % right
        0 -1;  % down
        0  1;
        -1 -1;
        1 -1;
        -1 1;
        1 1]; % up
problem = "one blk";

if problem == "entire domain"
    cf.new_airfoil = 1;
    cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
    cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
    cf.Nx1 = Nx1; cf.Nx2 = Nx2;
    cf.alpha = 1.005;
    cf.append_trail = 0;
    [x1, x2, M_samp, Mfun] = problems.InitProb5(cf);
    x1 = M_samp.x_metric;
    x2 = M_samp.y_metric;
elseif problem == "one blk"
    %cf.block_idx = 8;
    cf.block_idx = 5;
    cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
    cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
    [x1, x2, M11_fun,M12_fun,M22_fun] = problems.InitProb9(cf);
    M_samp.metric(:,:,1) = M11_fun(x1,x2);
    M_samp.metric(:,:,2) = M12_fun(x1,x2);
    M_samp.metric(:,:,3) = M22_fun(x1,x2);
    M_samp.x_metric = x1;
    M_samp.y_metric = x2;
    x1 = M_samp.x_metric;
    x2 = M_samp.y_metric;
elseif problem == "cartesian"
    Nx1 = 101; Nx2 = 101;
    x1 = linspace(0,1,Nx1);
    x2 = linspace(0,1,Nx2);
    [x1, x2] = ndgrid(x1,x2);
    x = cat(1,reshape(x1,1,Nx1,Nx2), reshape(x2,1,Nx1,Nx2));
    M_in = zeros(2,2,Nx1,Nx2);
    M_in(1,1,:,:) = 40000*(1+15*x1).^(-2);
    M_in(2,2,:,:) = 40000*(1+15*x2).^(-2);
end
%%

[Nx1, Nx2] = size(x1);
x = cat(1,reshape(x1,1,Nx1,Nx2), reshape(x2,1,Nx1,Nx2));
M = zeros(2,2,Nx1,Nx2);
M(1,1,:,:) = M_samp.metric(:,:,1);
M(1,2,:,:) = M_samp.metric(:,:,2);
M(2,1,:,:) = M_samp.metric(:,:,2);
M(2,2,:,:) = M_samp.metric(:,:,3);
x = x/100;
M_in = M;
M_out = zeros(2,2,Nx1,Nx2);

%%
for iter=1:5
for i=1:Nx1
    for j=1:Nx2
        p = x(:,i,j);
        M_p = M(:,:,i,j);

        for k = 1:4
            ii = i + nbr(k,1);
            jj = j + nbr(k,2);
    
            if ii < 1 || ii > Nx1 || jj < 1 || jj > Nx2
                continue
            end

            q = x(:,ii,jj);
            M_q = M(:,:,ii,jj);
            pq = q-p;
            %pq_eu = norm(pq,2);

            eta = (1 + sqrt(pq'*M_p*pq)*log(beta))^(-2);
            M_xp = eta*M_p;

            M_inter = inter(M_xp,M_q);
            M_q = M_inter;
            M(:,:,ii,jj) = M(:,:,ii,jj) + 0.5*(real(M_q)-M(:,:,ii,jj));
        end
    end
end
end
%%
figure
surf(x1,x2,reshape(M(1,1,:,:),Nx1,Nx2))
figure
surf(x1,x2,reshape(M_in(1,1,:,:),Nx1,Nx2))
%%
figure
plotMetricEllipses(M, x);
figure
plotMetricEllipses(M_in, x);
%%
S1 = linspace(0,1,Nx1);
S2 = linspace(0,1,Nx2);
[S1,S2] = ndgrid(S1,S2);
S = cat(1,reshape(S1,1,Nx1,Nx2), reshape(S2,1,Nx1,Nx2));

figure
plotMetricEllipses(M, S);
figure
plotMetricEllipses(M_in, S);

%%
function M_inter = inter(M_xp, M_q)
    %N = M_q \ M_xp;
    N = M_xp \ M_q;
    [V,D] = eig(N);
    V(:,1) = V(:,1) / norm(V(:,1));
    V(:,2) = V(:,2) / norm(V(:,2));
    lam_xp = V'*M_xp*V;
    lam_q = V'*M_q*V;
    lam1_inter = max(lam_xp(1,1), lam_q(1,1));
    lam2_inter = max(lam_xp(2,2), lam_q(2,2));

    M_inter = inv(V')*diag([lam1_inter, lam2_inter])*inv(V);
end


function h = plotMetricEllipses(M, X, varargin)
%PLOTMETRICELLIPSES Plot a cloud of 2D ellipses from a 2×2 metric field.
%   h = plotMetricEllipses(M, X) draws ellipses centered at positions X
%   whose principal axes align with the eigenvectors of the 2×2 SPD metric
%   tensors M, with axis lengths proportional to the eigenvalues.
%
%   INPUTS
%     M : [2, 2, Nx1, Nx2] array of symmetric positive-(semi)definite metrics
%     X : [2, Nx1, Nx2] array of ellipse centers (same grid as M)
%
%   PARAMETER-VALUE OPTIONS (all optional)
%     'Step'        : [s1 s2] stride for subsampling the grid (default auto)
%     'Npts'        : number of points per ellipse (default 64)
%     'Scale'       : scalar scale factor for axis lengths (default auto)
%     'Normalize'   : true|false. If true, scales so max axis ≈ 0.4*min spacing (default true)
%     'Color'       : line/edge color (default [0 0 0])
%     'LineWidth'   : line width (default 1.0)
%     'Filled'      : true|false to fill ellipses (default false)
%     'FaceAlpha'   : face transparency if Filled=true (default 0.15)
%
%   OUTPUT
%     h : array of line/patch handles for the drawn ellipses
%
%   EXAMPLE
%     Nx1=30; Nx2=20;
%     [x1,x2] = ndgrid(linspace(0,1,Nx1), linspace(0,1,Nx2));
%     X = cat(1, reshape(x1,1,Nx1,Nx2), reshape(x2,1,Nx1,Nx2));
%     % Anisotropic metric varying with x1:
%     lam1 = 50 + 30*sin(2*pi*x1); lam2 = 5 + 2*cos(2*pi*x2);
%     theta = pi/6 * ones(size(x1));
%     c=cos(theta); s=sin(theta);
%     M11 = c.^2.*lam1 + s.^2.*lam2;
%     M22 = s.^2.*lam1 + c.^2.*lam2;
%     M12 = c.*s.*(lam1 - lam2);
%     M = zeros(2,2,Nx1,Nx2);
%     M(1,1,:,:) = M11; M(2,2,:,:) = M22; M(1,2,:,:) = M12; M(2,1,:,:) = M12;
%     figure; axis equal; hold on;
%     plotMetricEllipses(M, X, 'Step', [2 2], 'Color', [0 0.3 0.8]);
%
%   NOTE
%     - If some metrics are not SPD (e.g., tiny negative eigenvalues from
%       roundoff), they are symmetrized and negative eigs are clamped to 0.
%     - Default scaling uses local grid spacing so ellipses are legible.

% -------------------------- parse inputs ---------------------------------
assert(ndims(M) == 4 && size(M,1)==2 && size(M,2)==2, 'M must be [2,2,Nx1,Nx2].');
assert(ndims(X) == 3 && size(X,1)==2, 'X must be [2,Nx1,Nx2].');
[Nx1, Nx2] = deal(size(M,3), size(M,4));
assert(isequal([size(X,2), size(X,3)], [Nx1, Nx2]), 'X grid must match M.');

p = inputParser;
p.addParameter('Step', [], @(v) (isnumeric(v) && (numel(v)==2 || isempty(v))));
p.addParameter('Npts', 64, @(v) isnumeric(v) && isscalar(v) && v>=16);
p.addParameter('Scale', [], @(v) isempty(v) || (isscalar(v) && v>0));
p.addParameter('Normalize', true, @(v) islogical(v) && isscalar(v));
p.addParameter('Color', [0 0 0], @(v) isnumeric(v) && numel(v)==3);
p.addParameter('LineWidth', 1.0, @(v) isnumeric(v) && isscalar(v) && v>0);
p.addParameter('Filled', false, @(v) islogical(v) && isscalar(v));
p.addParameter('FaceAlpha', 0.15, @(v) isscalar(v) && v>=0 && v<=1);
p.parse(varargin{:});
opt = p.Results;

% default stride based on grid size
if isempty(opt.Step)
    s = max(1, round(min(Nx1, Nx2)/30)); % heuristic density
    opt.Step = [s s];
else
    opt.Step = [opt.Step(:).' ];
end
s1 = opt.Step(1); s2 = opt.Step(2);

% ------------------------ prepare geometry --------------------------------
% Compute a heuristic local spacing to auto-scale ellipses
x1 = squeeze(X(1,:,:)); x2 = squeeze(X(2,:,:));
dx1 = []; dx2 = [];
if Nx1 > 1
    dx1 = vecnorm([diff(x1,1,1), zeros(Nx1-1,1)] + 1i*[diff(x2,1,1), zeros(Nx1-1,1)], 2, 2);
end
if Nx2 > 1
    dx2 = vecnorm([diff(x1,1,2); zeros(1,Nx2-1)] + 1i*[diff(x2,1,2); zeros(1,Nx2-1)], 2, 1);
end
gridLen = [];
if ~isempty(dx1), gridLen = [gridLen; dx1(:)]; end
if ~isempty(dx2), gridLen = [gridLen; dx2(:)]; end
if isempty(gridLen), gridLen = 1; end
hmin = prctile(gridLen, 10); % robust local spacing

% Precompute max eigenvalue for normalization if needed
lamMax = 1;
if opt.Normalize
    lamMax = 0;
    for i = 1:s1:Nx1
        for j = 1:s2:Nx2
            A = M(:,:,i,j);
            A = 0.5*(A + A.');         % symmetrize
            [V,D] = eig(A);
            lam = max(real(diag(D)), 0); % clamp
            lamMax = max(lamMax, max(lam));
        end
    end
    if lamMax <= 0, lamMax = 1; end
end

% Default scale: make the largest axis around 0.4*hmin
if isempty(opt.Scale)
    opt.Scale = 0.4 * hmin / lamMax;
end

% ellipse parameterization
theta = linspace(0, 2*pi, opt.Npts);
C = [cos(theta); sin(theta)];

% --------------------------- plotting -------------------------------------
holdState = ishold;
hold on;
axis equal;

h = gobjects(0);
for i = 1:s1:Nx1
    for j = 1:s2:Nx2
        xc = [X(1,i,j); X(2,i,j)];
        A = M(:,:,i,j);
        if any(~isfinite(A), 'all') || any(~isfinite(xc))
            continue;
        end
        A = 0.5*(A + A.');         % ensure symmetric
        [V,D] = eig(A);
        lam = max(real(diag(D)), 0);  % clamp negatives to zero

        % sort by descending eigenvalue for consistent major/minor axes
        [lam, idx] = sort(lam, 'descend');
        V = V(:, idx);

        % axis lengths proportional to eigenvalues
        axesLen = opt.Scale * lam(:).';
        if all(axesLen <= eps)
            continue;
        end

        R = V * diag(axesLen);
        P = xc + R * C;            % 2×Npts

        if opt.Filled
            h(end+1) = patch(P(1,:), P(2,:), opt.Color, ...
                             'FaceAlpha', opt.FaceAlpha, ...
                             'EdgeColor', opt.Color, ...
                             'LineWidth', opt.LineWidth); %#ok<AGROW>
        else
            h(end+1) = plot(P(1,:), P(2,:), '-', ...
                            'Color', opt.Color, ...
                            'LineWidth', opt.LineWidth); %#ok<AGROW>
        end
    end
end

if ~holdState, hold off; end
end
