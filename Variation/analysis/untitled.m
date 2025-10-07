x1 = linspace(0,1,40);
x2 = linspace(0,1,40);
[x1,x2] = meshgrid(x1,x2);

Mfun = GetM(1);
M = Mfun(x1,x2);
[eigvals, eigvecs] = eig2x2_metric(M);
h = plotMetricEigenvectors(eigvals, eigvecs, x1, x2, ColorMode='colormap', BaseScale=0.02);


function [eigvals, eigvecs] = eig2x2_metric(M)
    [Nx1,Nx2] = size(M.M11);
    eigvals = zeros(Nx1,Nx2,2);
    eigvecs = zeros(Nx1,Nx2,2,2);
    for i = 1:Nx1
        for j = 1:Nx2
            Mij = [M.M11(i,j), M.M12(i,j); M.M12(i,j), M.M22(i,j)];
            [V,D] = eig(Mij);
            eigvals(i,j,:) = diag(D);
            eigvecs(i,j,:,:) = V;
        end
    end
end


function h = plotMetricEigenvectors(lambdas, eigvecs, X, Y, varargin)
% plotMetricEigenvectors  Plot eigenvectors of a 2x2 SPD metric field
%
% h = plotMetricEigenvectors(lambdas, eigvecs, X, Y, Name,Value,...)
%
% Inputs
%   lambdas : N1 x N2 x 2    eigenvalues (lambda_1, lambda_2)
%   eigvecs : N1 x N2 x 2 x 2 eigenvectors: (:,:,1,k)=vx, (:,:,2,k)=vy, k=1,2
%   X, Y    : N1 x N2 grids
%
% Name-Value options
%   'Which'      : 1, 2, or 'both' (default 'both')
%   'Stride'     : positive integer (default 1)
%   'Centered'   : logical (default false)  % used in ColorMode="length"
%   'Scale'      : scalar (default 1)       % multiplies sqrt(lambda) in length mode
%   'ColorMode'  : 'length' (default) or 'colormap'
%   'BaseScale'  : scalar (default 1)       % constant length when ColorMode='colormap'
%   'CMap'       : Mx3 colormap (default parula(256))
%   'CLim'       : [vmin vmax] for coloring on sqrt(lambda); default = auto from data
%   'Color1'     : RGB for eigvec 1 in length mode (default [0 0.4470 0.7410])
%   'Color2'     : RGB for eigvec 2 in length mode (default [0.8500 0.3250 0.0980])
%   'LineWidth'  : line width (default 1.2)
%
% Output
%   h : struct of graphics handles

% ---- basic checks ----
szL = size(lambdas);
if numel(szL) ~= 3 || szL(3) ~= 2, error('lambdas must be N1 x N2 x 2.'); end
N1 = szL(1); N2 = szL(2);
if ~isequal(size(X), [N1 N2]) || ~isequal(size(Y), [N1 N2])
    error('X and Y must be N1 x N2 to match lambdas.');
end
szV = size(eigvecs);
if numel(szV) ~= 4 || ~isequal(szV(1:2), [N1 N2]) || szV(3) ~= 2 || szV(4) ~= 2
    error('eigvecs must be N1 x N2 x 2 x 2 (components x eigenvector index).');
end

% ---- options ----
p = inputParser;
addParameter(p, 'Which', 'both');
addParameter(p, 'Stride', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
addParameter(p, 'Centered', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'Scale', 1, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'ColorMode', 'length', @(s)ischar(s)||isstring(s));
addParameter(p, 'BaseScale', 1, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'CMap', parula(256));
addParameter(p, 'CLim', [], @(x)isnumeric(x)&&numel(x)==2);
addParameter(p, 'Color1', [0 0.4470 0.7410]);
addParameter(p, 'Color2', [0.8500 0.3250 0.0980]);
addParameter(p, 'LineWidth', 1.2, @(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p, varargin{:});
opt = p.Results;

% which eigenvectors
switch lower(string(opt.Which))
    case "1", whichList = 1;
    case "2", whichList = 2;
    case "both", whichList = [1 2];
    otherwise, error('Which must be 1, 2, or ''both''.');
end

% stride
s = opt.Stride; I = 1:s:N1; J = 1:s:N2;
Xs = X(I,J); Ys = Y(I,J);

% prepare output
h = struct('q1',[], 'q2',[], 'l1',[], 'l2',[], 'cmaps',[]);

holdState = ishold; hold on;

% ---------- COLOR MODE: length (previous behavior) ----------
if strcmpi(opt.ColorMode, 'length')
    for k = whichList
        vx = eigvecs(I,J,1,k); vy = eigvecs(I,J,2,k);
        L  = sqrt(lambdas(I,J,k)) * opt.Scale;
        U = L .* vx; V = L .* vy;

        if ~opt.Centered
            switch k
                case 1
                    h.q1 = quiver(Xs, Ys, U, V, 0, 'Color', opt.Color1, ...
                                  'LineWidth', opt.LineWidth, 'MaxHeadSize', 0.75);
                case 2
                    h.q2 = quiver(Xs, Ys, U, V, 0, 'Color', opt.Color2, ...
                                  'LineWidth', opt.LineWidth, 'MaxHeadSize', 0.75);
            end
        else
            X1 = Xs - 0.5*U;  Y1 = Ys - 0.5*V;
            X2 = Xs + 0.5*U;  Y2 = Ys + 0.5*V;
            npts = numel(X1);
            Xlines = [X1(:) X2(:) nan(npts,1)].';  % NaN-separated
            Ylines = [Y1(:) Y2(:) nan(npts,1)].';
            switch k
                case 1
                    h.l1 = plot(Xlines(:), Ylines(:), '-', 'Color', opt.Color1, 'LineWidth', opt.LineWidth);
                case 2
                    h.l2 = plot(Xlines(:), Ylines(:), '-', 'Color', opt.Color2, 'LineWidth', opt.LineWidth);
            end
        end
    end

% ---------- COLOR MODE: colormap (length fixed, color encodes sqrt(lambda)) ----------
elseif strcmpi(opt.ColorMode, 'colormap')
    % auto CLim on sqrt(lambda) if not provided
    Lvals = []; % collect values used for color
    for k = whichList
        aux = sqrt(lambdas(I,J,k));
        Lvals = [Lvals; aux(:)]; %#ok<AGROW>
    end
    if isempty(opt.CLim)
        vmin = min(Lvals); vmax = max(Lvals);
        if vmin==vmax, vmin = vmin-1e-12; vmax = vmax+1e-12; end
        clim = [vmin vmax];
    else
        clim = opt.CLim;
    end
    cmap = opt.CMap;
    nC = size(cmap,1);

    for k = whichList
        vx = eigvecs(I,J,1,k); vy = eigvecs(I,J,2,k);

        % Normalize eigenvectors to unit length, then apply constant BaseScale
        nrm = hypot(vx, vy); nrm(nrm==0) = 1;  % protect zeros
        vx = vx ./ nrm; vy = vy ./ nrm;

        U = opt.BaseScale * vx; V = opt.BaseScale * vy;

        % Build centered segments for consistent appearance
        X1 = Xs - 0.5*U;  Y1 = Ys - 0.5*V;
        X2 = Xs + 0.5*U;  Y2 = Ys + 0.5*V;

        % Color by sqrt(lambda_k)
        val = sqrt(lambdas(I,J,k));
        t = (val - clim(1)) ./ (clim(2)-clim(1));
        t = max(0,min(1,t));                % clamp
        idx = 1 + floor(t*(nC-1));          % 1..nC
        % Draw each segment with its own RGB (loop is OK with stride)
        nSeg = numel(X1);
        Lh = gobjects(nSeg,1);
        for n = 1:nSeg
            Lh(n) = line([X1(n), X2(n)], [Y1(n), Y2(n)], ...
                         'Color', cmap(idx(n),:), 'LineWidth', opt.LineWidth);
        end
        if k==1, h.l1 = Lh; else, h.l2 = Lh; end
    end

    % add a colorbar that reflects sqrt(lambda)
    caxis(clim); colormap(cmap);
    cb = colorbar; cb.Label.String = '\surd\lambda';
    h.cbar = cb;

else
    error('ColorMode must be ''length'' or ''colormap''.');
end

if ~holdState, hold off; end
axis equal tight; box on
xlabel('x'); ylabel('y');

% title + legend
if strcmpi(opt.ColorMode,'length')
    ttl = 'Eigenvectors scaled by \surd\lambda';
else
    ttl = 'Eigenvectors colored by \surd\lambda (constant length)';
end
title(ttl);

legendEntries = {};
if any(whichList==1) && strcmpi(opt.ColorMode,'length')
    legendEntries{end+1} = 'eigvec 1 (\surd\lambda_1 length)'; %#ok<AGROW>
elseif any(whichList==1)
    legendEntries{end+1} = 'eigvec 1 (colored by \surd\lambda_1)';
end
if any(whichList==2) && strcmpi(opt.ColorMode,'length')
    legendEntries{end+1} = 'eigvec 2 (\surd\lambda_2 length)';
elseif any(whichList==2)
    legendEntries{end+1} = 'eigvec 2 (colored by \surd\lambda_2)';
end
if ~isempty(legendEntries) && strcmpi(opt.ColorMode,'length')
    legend(legendEntries, 'Location', 'bestoutside');
end
end
