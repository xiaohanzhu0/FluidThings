function [p, t, model] = triangulate_from_edgeArcs(edgeArcs, varargin)
%TRIANGULATE_FROM_EDGEARCS  Build PDE geometry from boundary arcs and mesh it.
%  [p,t,model] = triangulate_from_edgeArcs(edgeArcs, 'Hmax', hmax, ...)
%
%  INPUT (required)
%    edgeArcs : struct describing the outer boundary as four ordered arcs
%               and (optionally) holes. Fields:
%       .g1, .g2, .g3, .g4 : (Ni x 2) arrays of [x y] points along each arc
%                            in traversal order around the domain boundary.
%                            Endpoints of consecutive arcs must coincide.
%       .holes (optional)   : cell array {H1, H2, ...}, each Hj is (Nj x 2)
%                             point list for one hole boundary (clockwise or
%                             counter-clockwise; will be handled).
%  NAME-VALUE OPTIONS (any subset)
%       'Hmax'            : target maximum element size (default: auto)
%       'Hmin'            : minimum element size (default: [])
%       'Hgrad'           : mesh growth rate (default: 1.5)
%       'GeometricOrder'  : 'linear' (default) or 'quadratic'
%
%  OUTPUT
%    p     : (nP x 2) node coordinates
%    t     : (nT x 3) triangle connectivity (1-based indices into rows of p)
%    model : PDEModel object (contains geometry + mesh if you need it)
%
%  NOTES
%   • This function constructs a 2-D analytic geometry via DECSG from your
%     polylines, then calls generateMesh (PDE Toolbox) to triangulate.
%   • If you have more than four outer arcs, concatenate them beforehand or
%     adapt the code where noted.
%   • After meshing, you can obtain PDE edge IDs via pdegplot(model,'EdgeLabels','on').
%
%  Example:
%    EA.g1 = [x1 y1; ...]; EA.g2 = [...]; EA.g3 = [...]; EA.g4 = [...];
%    EA.holes = {H1_pts, H2_pts};   % optional
%    [p,t] = triangulate_from_edgeArcs(EA,'Hmax',0.02);
%
%  Xiaohan-ready: returns (p,t) just like you asked.

% ---------- parse options ----------
Hmax = []; Hmin = []; Hgrad = 1.5; GeometricOrder = 'linear';
for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'hmax',  Hmax = varargin{k+1};
        case 'hmin',  Hmin = varargin{k+1};
        case 'hgrad', Hgrad = varargin{k+1};
        case 'geometricorder', GeometricOrder = varargin{k+1};
        otherwise, error('Unknown option: %s', varargin{k});
    end
end

% ---------- build outer polygon from arcs ----------
assert(all(isfield(edgeArcs, {'g1','g2','g3','g4'})), ...
    'edgeArcs must contain g1,g2,g3,g4 as polylines.');

P1 = edgeArcs.g1; P2 = edgeArcs.g2; P3 = edgeArcs.g3; P4 = edgeArcs.g4;

% drop duplicate endpoints when concatenating
outer = [ P1; P2(2:end,:); P3(2:end,:); P4(2:end,:) ];
outer = remove_duplicate_vertices(outer, 1e-12);

% ensure consistent orientation (outer boundary CCW is conventional)
if signed_area(outer) < 0
    outer = flipud(outer);
end

% ---------- handle holes (if any) ----------
holeList = {};
if isfield(edgeArcs,'holes') && ~isempty(edgeArcs.holes)
    holeList = edgeArcs.holes(:)';
    % make hole loops clockwise (recommended when subtracting)
    for i = 1:numel(holeList)
        if signed_area(holeList{i}) > 0
            holeList{i} = flipud(holeList{i});
        end
        holeList{i} = remove_duplicate_vertices(holeList{i}, 1e-12);
    end
end

% ---------- build DECSG geometry description ----------
polyCells = [{outer}, holeList];
maxNV = max(cellfun(@(Q) size(Q,1), polyCells));
M = 2 + 2*maxNV;               % column height for one polygon in DECSG gd
cols = zeros(M, numel(polyCells)); cols(:) = NaN;

for j = 1:numel(polyCells)
    Q = polyCells{j}; N = size(Q,1);
    col = NaN(M,1);
    col(1) = 2;                % polygon code for DECSG
    col(2) = N;                % number of vertices
    col(3:2+N) = Q(:,1);       % x's (pad with NaN automatically)
    col(3+maxNV:2+maxNV+N) = Q(:,2); % y's after padding block
    cols(:,j) = col;
end

% names & set formula: outer P1 minus all holes Hk
names = ['P1'; char(arrayfun(@(k) sprintf('H%d',k), 1:numel(holeList), 'uni', false))];
ns = names';                    % DECSG expects names as columns
sf = 'P1';
for k = 1:numel(holeList)
    sf = sprintf('%s - H%d', sf, k);
end

[dl, ~] = decsg(cols, sf, ns);

% ---------- build PDE model & generate mesh ----------
model = createpde(1);
geometryFromEdges(model, dl);

meshArgs = {'GeometricOrder', GeometricOrder, 'Hgrad', Hgrad};
if ~isempty(Hmax), meshArgs = [meshArgs, {'Hmax', Hmax}]; end
if ~isempty(Hmin), meshArgs = [meshArgs, {'Hmin', Hmin}]; end
%meshArgs = [meshArgs, {'Hedge', {40:79,0.01}}];
%meshArgs = [meshArgs, {'Hedge', {[1,2],0.001}}];

msh = generateMesh(model, meshArgs{:}); %#ok<NASGU>

% ---------- extract p,t (transpose to row-wise form) ----------
p = model.Mesh.Nodes';                 % nP x 2
T = model.Mesh.Elements;               % 3 x nT (1-based)
if size(T,1) > 3
    % quadratic elements present
    T = T(1:3,:);                      % keep vertex connectivity only
end
t = T';                                % nT x 3

end

% ================= helper functions =================
function A = signed_area(P)
% signed polygon area (positive for CCW)
    x = P(:,1); y = P(:,2);
    x2 = [x(2:end); x(1)]; y2 = [y(2:end); y(1)];
    A = 0.5 * sum(x.*y2 - x2.*y);
end

function Q = remove_duplicate_vertices(P, tol)
% remove immediate duplicates and collapse coincident joins
    if nargin < 2, tol = 0; end
    keep = true(size(P,1),1);
    for i = 2:size(P,1)
        if norm(P(i,:) - P(i-1,:)) <= tol
            keep(i) = false;
        end
    end
    % also check last-first
    if norm(P(end,:) - P(1,:)) <= tol
        keep(end) = false; % not needed by DECSG, polygon closes implicitly
    end
    Q = P(keep,:);
end
