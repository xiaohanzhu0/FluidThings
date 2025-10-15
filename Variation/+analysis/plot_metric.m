function [h1,h2,h3] = plot_metric(x1, x2, M)
    [eigvals, eigvecs] = eig2x2_metric(M);
    figure; colorbar;
    h1 = plotMetricEigenvectors(eigvals, eigvecs, x1, x2, ColorMode='colormap', BaseScale=0.02, ColorBy='sqrtlambda', CMap=winter(256));
    figure; colorbar;
    h2 = plotMetricEigenvectors(eigvals, eigvecs, x1, x2, ColorMode='colormap', BaseScale=0.02, ColorBy='sqrt_ratio', CMap=copper(256));
    figure; colorbar;
    h3 = plotMetricEigenvectors(eigvals, eigvecs, x1, x2, ColorMode='colormap', BaseScale=0.02, ColorBy='axis_deviation', CMap=autumn(256));
end

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

            %gap = abs(D(1,1) - D(2,2));
            %scale = max(1, max(abs(D),[],'all'));
            %if gap <= 1e-4 * scale
            %    % Nearly repeated eigenvalues: use basis-aligned directions
            %    eigvecs(i,j,:,:) = eye(2);
            %end

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
%   'ColorBy'    : 'sqrtlambda' (default), 'sqrt_ratio', 'axis_deviation'  % only when ColorMode='colormap'
%   'BaseScale'  : scalar (default 1)       % constant length in colormap mode
%   'CMap'       : Mx3 colormap (default parula(256))
%   'CLim'       : [vmin vmax] for coloring scalar; default auto from data
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
addParameter(p, 'ColorBy', 'sqrtlambda', @(s)ischar(s)||isstring(s));  % only for 'colormap'
addParameter(p, 'BaseScale', 1, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'CMap', winter(256));
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
h = struct('q1',[], 'q2',[], 'l1',[], 'l2',[], 'cbar',[]);

holdState = ishold; hold on;

% ---------- COLOR MODE: length ----------
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

% ---------- COLOR MODE: colormap ----------
elseif strcmpi(opt.ColorMode, 'colormap')
    % constant length segments; color carries the scalar info
    cmap = opt.CMap; nC = size(cmap,1);

    % Precompute anisotropy (sqrt_ratio) which is shared across k
    lam1 = lambdas(I,J,1);
    lam2 = lambdas(I,J,2);
    lam_min = min(lam1, lam2);
    lam_max = max(lam1, lam2);
    % protect against exact zeros (shouldn't happen for SPD, but be safe)
    lam_min(lam_min==0) = eps;
    sqrt_ratio_val = sqrt(lam_max ./ lam_min);

    % Gather all scalar values (depending on ColorBy) to set CLim automatically
    allVals = [];

    % First pass: compute scalar fields for each k as needed
    % We'll store handles per k after drawing
    scalarPerK = cell(1,2);

    for k = whichList
        switch lower(string(opt.ColorBy))
            case "sqrtlambda"
                scalarPerK{k} = sqrt(lambdas(I,J,k));
            case "sqrt_ratio"
                scalarPerK{k} = sqrt_ratio_val; % same for both eigenvectors at a point
            case "axis_deviation"
                vx = eigvecs(I,J,1,k); vy = eigvecs(I,J,2,k);
                % Normalize direction (protect zeros)
                nrm = hypot(vx, vy); nrm(nrm==0) = 1;
                vx = vx ./ nrm; vy = vy ./ nrm;
                theta = atan2(vy, vx);           % [-pi, pi]
                % Deviation from nearest axis (period pi/2):
                th_mod = mod(theta, pi/2);       % [0, pi/2)
                dev = min(th_mod, (pi/2) - th_mod); % [0, pi/4]
                scalarPerK{k} = dev / (pi/4);    % normalize to [0,1]
            otherwise
                error('ColorBy must be ''sqrtlambda'', ''sqrt_ratio'', or ''axis_deviation''.');
        end
        allVals = [allVals; scalarPerK{k}(:)]; %#ok<AGROW>
    end

    % Determine CLim
    if isempty(opt.CLim)
        vmin = min(allVals); vmax = max(allVals);
        if vmin==vmax, vmin = vmin-1e-12; vmax = vmax+1e-12; end
        clim = [vmin vmax];
    else
        clim = opt.CLim;
    end

    % Draw per k
    for k = whichList
        vx = eigvecs(I,J,1,k); vy = eigvecs(I,J,2,k);
        % Normalize direction then set constant length
        nrm = hypot(vx, vy); nrm(nrm==0) = 1;
        vx = vx ./ nrm; vy = vy ./ nrm;

        U = opt.BaseScale * vx; V = opt.BaseScale * vy;

        % Centered segments look best with per-segment colors
        X1 = Xs - 0.5*U;  Y1 = Ys - 0.5*V;
        X2 = Xs + 0.5*U;  Y2 = Ys + 0.5*V;

        val = scalarPerK{k};
        % map values to colormap indices
        t = (val - clim(1)) ./ (clim(2)-clim(1));
        t = max(0,min(1,t));
        idx = 1 + floor(t*(nC-1));  % 1..nC

        nSeg = numel(X1);
        Lh = gobjects(nSeg,1);
        for n = 1:nSeg
            Lh(n) = line([X1(n), X2(n)], [Y1(n), Y2(n)], ...
                         'Color', cmap(idx(n),:), 'LineWidth', opt.LineWidth);
        end
        if k==1, h.l1 = Lh; else, h.l2 = Lh; end
    end

    % Colorbar label depending on ColorBy
    caxis(clim); colormap(cmap);
    cb = colorbar;
    switch lower(string(opt.ColorBy))
        case "sqrtlambda"
            cb.Label.String = '\surd\lambda_k';
        case "sqrt_ratio"
            cb.Label.String = '\surd(\lambda_{max}/\lambda_{min})  (anisotropy)';
        case "axis_deviation"
            cb.Label.String = 'axis deviation (0=axis, 1=45^\circ)';
    end
    h.cbar = cb;

else
    error('ColorMode must be ''length'' or ''colormap''.');
end

if ~holdState, hold off; end
axis equal tight; box on
%xlabel('x'); ylabel('y');

% title + legend
if strcmpi(opt.ColorMode,'length')
    ttl = 'Eigenvectors scaled by \surd\lambda';
else
    switch lower(string(opt.ColorBy))
        case "sqrtlambda",   ttl = 'Eigenvectors colored by \surd\lambda (constant length)';
        case "sqrt_ratio",   ttl = 'Eigenvectors colored by anisotropy \surd(\lambda_{max}/\lambda_{min})';
        case "axis_deviation", ttl = 'Eigenvectors colored by axis deviation (0=axis, 1=45^\circ)';
    end
end
title(ttl);

% Only show legend in length mode (colors are uniform there)
if strcmpi(opt.ColorMode,'length')
    legendEntries = {};
    if any(whichList==1), legendEntries{end+1} = 'eigvec 1'; end %#ok<AGROW>
    if any(whichList==2), legendEntries{end+1} = 'eigvec 2'; end %#ok<AGROW>
    if ~isempty(legendEntries)
        legend(legendEntries, 'Location', 'bestoutside');
    end
end
end
