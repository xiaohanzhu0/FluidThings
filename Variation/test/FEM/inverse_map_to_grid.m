function [Xs, Ys, S1g, S2g, info] = inverse_map_to_grid(p, t, s1, s2, N_s1, N_s2, varargin)
%INVERSE_MAP_TO_GRID  Interpolate inverse map x(s) from nodal s(x) on a tri mesh.
%
%  [Xs, Ys, S1g, S2g, info] = inverse_map_to_grid(p, t, s1, s2, N_s1, N_s2, ...)
%
% Inputs
%   p     : (nP x 2) node coordinates [x y] of the physical FEM mesh
%   t     : (nT x 3) triangles (1-based). (Not required by default path)
%   s1    : (nP x 1) nodal values of s_1(x)
%   s2    : (nP x 1) nodal values of s_2(x)
%   N_s1  : number of grid points along s_1 (columns)
%   N_s2  : number of grid points along s_2 (rows)
%
% Name-Value options
%   'S1Range' : [s1_min s1_max] (default [0 1])
%   'S2Range' : [s2_min s2_max] (default [0 1])
%   'Method'  : 'natural' (default) | 'linear' | 'nearest'
%   'Fill'    : 'nearest' (default) | 'none'
%               How to fill NaNs outside convex hull of (s1,s2).
%
% Outputs
%   Xs, Ys : (N_s2 x N_s1) arrays s.t. Xs(i,j),Ys(i,j) = x(s1_j, s2_i)
%   S1g,S2g: the (N_s2 x N_s1) grids of s1 and s2 used for evaluation
%   info   : diagnostics (nan counts, duplicates, etc.)
%
% Notes
%   • This assumes the mapping is one-to-one (no folds). If your mapping
%     folds, some (s1,s2) targets map to multiple x; this method will pick
%     one branch implicitly via Delaunay in (s1,s2).
%   • Grid orientation: consistent with meshgrid. Columns index s1, rows s2.

% --------- parse options ----------
opts.S1Range = [0, 1];
opts.S2Range = [0, 1];
opts.Method  = 'natural';   % 'natural' = natural neighbor (good quality)
opts.Fill    = 'nearest';   % fill NaNs by nearest in (s1,s2)
opts = parseopts(opts, varargin{:});

% --------- sanity checks ----------
nP = size(p,1);
assert(isvector(s1) && isvector(s2) && numel(s1)==nP && numel(s2)==nP, ...
    's1,s2 must be nP-by-1 vectors aligned with p.');
s1 = s1(:); s2 = s2(:);
x  = p(:,1); y = p(:,2);

% --------- build target (s1,s2) grid ----------
s1v = linspace(opts.S1Range(1), opts.S1Range(2), N_s1);
s2v = linspace(opts.S2Range(1), opts.S2Range(2), N_s2);
[S1g, S2g] = meshgrid(s1v, s2v);   % size N_s2 x N_s1

% --------- handle duplicate (s1,s2) nodes (rare, e.g., shared corners) ----------
S = [s1, s2];
[Su, ia, ic] = unique(S, 'rows', 'stable');
xu = accumarray(ic, x, [], @mean);
yu = accumarray(ic, y, [], @mean);
numDup = size(S,1) - size(Su,1);

% --------- primary interpolants in (s1,s2) space ----------
Fx = scatteredInterpolant(Su(:,1), Su(:,2), xu, opts.Method, 'none');
Fy = scatteredInterpolant(Su(:,1), Su(:,2), yu, opts.Method, 'none');

Xs = Fx(S1g, S2g);
Ys = Fy(S1g, S2g);

% --------- fill NaNs near the hull if requested ----------
nanMask = isnan(Xs) | isnan(Ys);
if any(nanMask,'all') && ~strcmpi(opts.Fill,'none')
    FxN = scatteredInterpolant(Su(:,1), Su(:,2), xu, 'nearest', 'nearest');
    FyN = scatteredInterpolant(Su(:,1), Su(:,2), yu, 'nearest', 'nearest');
    Xs(nanMask) = FxN(S1g(nanMask), S2g(nanMask));
    Ys(nanMask) = FyN(S1g(nanMask), S2g(nanMask));
end

info = struct();
info.duplicatesCollapsed = numDup;
info.nanBeforeFill = nnz(isnan(Fx(S1g,S2g)) | isnan(Fy(S1g,S2g)));
info.nanAfterFill  = nnz(isnan(Xs) | isnan(Ys));
info.method = opts.Method;
info.fill   = opts.Fill;
info.S1Range = opts.S1Range;
info.S2Range = opts.S2Range;

end

% ================= helper: parseopts =================
function opts = parseopts(opts, varargin)
    if mod(numel(varargin),2)~=0
        error('Name-value pairs expected.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k}; val = varargin{k+1};
        if ~isfield(opts, name)
            error('Unknown option: %s', name);
        end
        opts.(name) = val;
    end
end
