%% smooth_metric_admm_curvilinear_edgeLipschitz.m
% This script smooths a sampled 2D SPD metric field M(x,y) on a *curvilinear
% structured grid* using an ADMM solver with *physical Lipschitz per edge*
% constraints in Log–Euclidean space.
%
% INPUTS (all 2D arrays of the same size Ny-by-Nx):
%   X   : x-position of grid nodes
%   Y   : y-position of grid nodes
%   M11 : metric tensor component (1,1) at nodes
%   M12 : metric tensor component (1,2) = (2,1) at nodes
%   M22 : metric tensor component (2,2) at nodes
%
% OUTPUTS:
%   M11_s, M12_s, M22_s : smoothed metric components at nodes
%
% MODEL (edge Lipschitz in physical space):
%   Let V = log(M) (matrix log), represented as a 3-vector per node:
%       V = [V11; V12; V22]
%   Solve:
%       minimize   1/2 sum_nodes ||V - U||_W^2      (U = log(M0))
%       subject to ||V_q - V_p||_2 <= beta * ell_pq for each physical edge (p,q)
%   where ell_pq is the physical edge length between nodes p and q.
%
% ADMM splitting (per edge e=(p,q)):
%   z_e  ≈ (V_q - V_p),   with constraint ||z_e|| <= beta*ell_e
%   y_e  is the scaled dual variable.
%
% -------------------------------------------------------------------------
% IMPORTANT:
% - This assumes the input metric is SPD at every node.
% - You may want to clamp/fix small eigenvalues before log if your data is noisy.
% -------------------------------------------------------------------------

function metric_s = metric_gradation2(X, Y, metric, params)
%% ---------------------- USER PARAMETERS ---------------------------------
beta    = params.beta;    % gradation / Lipschitz strength (units: "per unit physical length")
rho     = 5.0;    % ADMM penalty parameter (tune for speed; see notes below)
maxIter = params.grade_iterations;    % ADMM iterations
tolPri  = 1e-4;   % primal residual tolerance
tolDual = 1e-4;   % dual residual tolerance

usePCG  = true;   % true: use PCG; false: use direct sparse backslash

%% ---------------------- LOAD / CHECK INPUTS -----------------------------
% You should have X,Y,M11,M12,M22 already in your workspace.
% Uncomment the following lines if you want to load from a MAT file.
% load('metric_data.mat','X','Y','M11','M12','M22');
M11 = metric(:,:,1);
M12 = metric(:,:,2);
M22 = metric(:,:,3);

[Ny,Nx] = size(X);
assert(all(size(Y)   == [Ny,Nx]), 'Y must be Ny-by-Nx');
assert(all(size(M11) == [Ny,Nx]), 'M11 must be Ny-by-Nx');
assert(all(size(M12) == [Ny,Nx]), 'M12 must be Ny-by-Nx');
assert(all(size(M22) == [Ny,Nx]), 'M22 must be Ny-by-Nx');

N = Ny*Nx;  % total number of nodes

%% ---------------------- STEP 1: LOG OF INPUT METRIC ---------------------
% Work in Log–Euclidean space: U = log(M0).
% Represent each symmetric 2x2 matrix by a 3-vector [a; b; c] for [a b; b c].
U = zeros(Ny, Nx, 3);  % U(:,:,1)=U11, U(:,:,2)=U12, U(:,:,3)=U22

for i = 1:Ny
    for j = 1:Nx
        U(i,j,:) = logSPD2(M11(i,j), M12(i,j), M22(i,j));
    end
end

% Initialize V with U (a natural warm start).
V = U;

%% ---------------------- STEP 2: EDGE GEOMETRY (LENGTHS) -----------------
% We use 4-neighbor edges on the structured index grid:
%   - horizontal edges: (i,j) <-> (i,j+1) for j=1..Nx-1
%   - vertical edges:   (i,j) <-> (i+1,j) for i=1..Ny-1
%
% Physical edge lengths come from X,Y:
ellH = hypot( X(:,2:Nx) - X(:,1:Nx-1),  Y(:,2:Nx) - Y(:,1:Nx-1) );  % Ny-by-(Nx-1)
ellV = hypot( X(2:Ny,:) - X(1:Ny-1,:),  Y(2:Ny,:) - Y(1:Ny-1,:) );  % (Ny-1)-by-Nx

% Radii for Lipschitz constraints: ||Vq - Vp|| <= beta * ell
rH = beta * ellH;
rV = beta * ellV;

%% ---------------------- STEP 3: ADMM VARIABLES (z, y) -------------------
% For each edge we store a 3-vector:
%   zH(i,j,:) corresponds to edge (i,j) -> (i,j+1)
%   zV(i,j,:) corresponds to edge (i,j) -> (i+1,j)
zH = zeros(Ny,   Nx-1, 3);
zV = zeros(Ny-1, Nx,   3);

% Scaled dual variables y (same shapes as z).
yH = zeros(Ny,   Nx-1, 3);
yV = zeros(Ny-1, Nx,   3);

% A reasonable initialization is to set z to the projected initial differences.
[dVH, dVV] = edgeDifferences(V);    % raw edge differences from current V
zH = projectEdgeBall(dVH, rH);      % per-edge projection onto ball radius rH
zV = projectEdgeBall(dVV, rV);

% Duals start at zero.
% yH,yV already zero.

%% ---------------------- STEP 4: BUILD SPARSE LINEAR SYSTEM --------------
% The V-update solves (component-wise):
%   V^{k+1} = argmin_V  1/2||V-U||_W^2 + (rho/2)||D V - z + y||^2
%
% This yields a sparse SPD system (per component):
%   (Wcomp * I + rho * L) vcomp = Wcomp * ucomp + rho * D^T (z - y)
%
% where L = D^T D is the (unweighted) graph Laplacian for 4-neighbor grid edges.
%
% Objective weight W for symmetric 2x2 matrices:
%   ||V-U||_W^2 = (V11-U11)^2 + 2*(V12-U12)^2 + (V22-U22)^2
% so:
wComp = [1, 2, 1];   % per-component data fidelity weight

% Build a single Laplacian L on the node graph induced by the structured edges.
% We'll assemble L as a sparse matrix in the flattened node ordering.
% Flattening index:
%   idx(i,j) = i + (j-1)*Ny   (column-major like MATLAB)
L = buildStructuredLaplacian(Ny, Nx);

% Build the three SPD matrices A_c = wComp(c)*I + rho*L.
A = cell(1,3);
for c = 1:3
    A{c} = spdiags(wComp(c)*ones(N,1), 0, N, N) + rho * L;
end

% If using PCG, build an incomplete Cholesky preconditioner for each component.
if usePCG
    setup.type = 'ict';
    setup.droptol = 1e-3;
    P = cell(1,3);
    for c = 1:3
        % ichol expects SPD; A{c} is SPD for rho>0 and wComp(c)>0
        P{c} = ichol(A{c}, setup);
    end
end

%% ---------------------- STEP 5: ADMM ITERATIONS -------------------------
fprintf('ADMM start: Ny=%d, Nx=%d, nodes=%d\n', Ny, Nx, N);
fprintf('%5s  %12s  %12s\n', 'iter', 'r_pri', 'r_dual');

% Store previous z for dual residual.
zH_prev = zH;
zV_prev = zV;

for k = 1:maxIter

    % ---------------------- (1) V-UPDATE --------------------------------
    % Compute rhs = W*u + rho * D^T (z - y), component-wise.
    %
    % First form edge fields (z - y) for each component.
    dH = zH - yH;   % Ny-by-(Nx-1)-by-3
    dV = zV - yV;   % (Ny-1)-by-Nx-by-3

    % Compute node-wise divergence: div = D^T(d) (same shape as V).
    div = edgeDivergence(dH, dV, Ny, Nx);  % Ny-by-Nx-by-3

    % Solve three independent SPD systems (one per component).
    for c = 1:3
        rhs = wComp(c) * U(:,:,c) + rho * div(:,:,c);
        rhs = rhs(:);  % flatten to N-by-1

        if usePCG
            % Preconditioned Conjugate Gradient solve.
            % Tight tolerance isn't usually needed because ADMM tolerances dominate.
            [vsol,flag] = pcg(A{c}, rhs, 1e-8, 200, P{c}, P{c}');
            if flag ~= 0
                warning('PCG did not fully converge for component %d at iter %d (flag=%d).', c, k, flag);
            end
        else
            % Direct sparse solve (can be expensive for large grids).
            vsol = A{c} \ rhs;
        end

        V(:,:,c) = reshape(vsol, Ny, Nx);
    end

    % ---------------------- (2) z-UPDATE (local projections) -------------
    % Compute current edge differences DV.
    [dVH, dVV] = edgeDifferences(V);  % shapes match zH and zV

    % Each edge solves:
    %   z = argmin_{||z||<=r} (rho/2)|| (DV + y) - z ||^2
    % i.e. z = projection of (DV + y) onto the l2 ball radius r.
    zH = projectEdgeBall(dVH + yH, rH);
    zV = projectEdgeBall(dVV + yV, rV);

    % ---------------------- (3) y-UPDATE (scaled dual) -------------------
    % y := y + (DV - z)
    yH = yH + (dVH - zH);
    yV = yV + (dVV - zV);

    % ---------------------- (4) STOPPING CRITERIA ------------------------
    % Primal residual: r_pri = ||DV - z||
    rPri = sqrt( sum((dVH - zH).^2,'all') + sum((dVV - zV).^2,'all') );

    % Dual residual (scaled form):
    % r_dual = rho * ||D^T (z - z_prev)||
    dzH = zH - zH_prev;
    dzV = zV - zV_prev;
    div_dz = edgeDivergence(dzH, dzV, Ny, Nx);  % D^T(dz)
    rDual = rho * sqrt( sum(div_dz.^2,'all') );

    if mod(k,10)==0 || k==1
        fprintf('%5d  %12.4e  %12.4e\n', k, rPri, rDual);
    end

    if (rPri <= tolPri) && (rDual <= tolDual)
        fprintf('Converged at iter %d: r_pri=%.3e, r_dual=%.3e\n', k, rPri, rDual);
        break;
    end

    % Update z_prev for next dual residual.
    zH_prev = zH;
    zV_prev = zV;

end

%% ---------------------- STEP 6: EXP BACK TO SPD METRIC ------------------
% Convert V = log(M) back to M via matrix exponential at each node.
M11_s = zeros(Ny,Nx);
M12_s = zeros(Ny,Nx);
M22_s = zeros(Ny,Nx);

for i = 1:Ny
    for j = 1:Nx
        [a,b,c] = expSPD2(V(i,j,1), V(i,j,2), V(i,j,3));
        M11_s(i,j) = a;
        M12_s(i,j) = b;
        M22_s(i,j) = c;
    end
end
metric_s = zeros(size(metric));
metric_s(:,:,1) = M11_s;
metric_s(:,:,2) = M12_s;
metric_s(:,:,3) = M22_s;
% You now have M11_s, M12_s, M22_s in the workspace.
fprintf('Done. Smoothed metric returned as M11_s, M12_s, M22_s.\n');
end

%% ====================== HELPER FUNCTIONS =================================
function L = buildStructuredLaplacian(Ny, Nx)
% Build the unweighted graph Laplacian for a Ny-by-Nx grid with 4-neighbor edges.
% Node indexing (column-major):
%   id(i,j) = i + (j-1)*Ny, with i=1..Ny, j=1..Nx.

N = Ny*Nx;

% We'll assemble using triplets (I,J,V) for sparse().
% Each undirected edge contributes:
%   L(p,p) += 1, L(q,q) += 1, L(p,q) -= 1, L(q,p) -= 1
I = [];
J = [];
V = [];

% Horizontal edges: (i,j) -- (i,j+1)
for j = 1:(Nx-1)
    for i = 1:Ny
        p = i + (j-1)*Ny;
        q = i + (j)*Ny;

        I = [I; p; q; p; q];
        J = [J; p; q; q; p];
        V = [V; 1; 1; -1; -1];
    end
end

% Vertical edges: (i,j) -- (i+1,j)
for j = 1:Nx
    for i = 1:(Ny-1)
        p = i     + (j-1)*Ny;
        q = (i+1) + (j-1)*Ny;

        I = [I; p; q; p; q];
        J = [J; p; q; q; p];
        V = [V; 1; 1; -1; -1];
    end
end

L = sparse(I, J, V, N, N);
end

function [dH, dV] = edgeDifferences(V)
% Compute forward differences along structured edges.
% V is Ny-by-Nx-by-3.
% dH(i,j,:) = V(i,j+1,:) - V(i,j,:) for j=1..Nx-1
% dV(i,j,:) = V(i+1,j,:) - V(i,j,:) for i=1..Ny-1
dH = V(:,2:end,:) - V(:,1:end-1,:);
dV = V(2:end,:,:) - V(1:end-1,:,:);
end

function div = edgeDivergence(dH, dV, Ny, Nx)
% Compute D^T(d) where dH and dV are edge fields:
%   dH: Ny-by-(Nx-1)-by-3 for horizontal edges (i,j)->(i,j+1)
%   dV: (Ny-1)-by-Nx-by-3 for vertical edges (i,j)->(i+1,j)
%
% For each horizontal edge (p->q):
%   add -d to node p, +d to node q
% For each vertical edge similarly.
div = zeros(Ny, Nx, 3);

% Horizontal contributions
div(:,1:end-1,:) = div(:,1:end-1,:) - dH;
div(:,2:end,  :) = div(:,2:end,  :) + dH;

% Vertical contributions
div(1:end-1,:,:) = div(1:end-1,:,:) - dV;
div(2:end,  :,:) = div(2:end,  :,:) + dV;
end

function z = projectEdgeBall(w, r)
% Project per-edge vectors w onto an l2 ball of radius r (elementwise).
% w: ...-by-3 array, r: ... array (same leading dimensions).
% Output z has same size as w.
%
% Projection:
%   if ||w|| <= r: z=w
%   else          z=(r/||w||)*w

% Compute norms across the 3 components.
nw = sqrt( sum(w.^2, 3) );     % leading dims (Ny-by-(Nx-1) or (Ny-1)-by-Nx)

% Avoid division by zero.
eps0 = 1e-14;

% Scaling factor s = min(1, r/nw)
s = ones(size(nw));
mask = (nw > r);
s(mask) = r(mask) ./ max(nw(mask), eps0);

% Apply scaling to each component.
z = w;
z(:,:,1) = s .* w(:,:,1);
z(:,:,2) = s .* w(:,:,2);
z(:,:,3) = s .* w(:,:,3);
end

function v = logSPD2(a,b,c)
% Compute vecsym(log([a b; b c])) for a 2x2 SPD matrix.
% Uses eigen-decomposition (stable for SPD).
M = [a b; b c];
% Symmetrize defensively
M = 0.5*(M+M');
[R,D] = eig(M);
lam = diag(D);

% If your input might be nearly singular/noisy, clamp eigenvalues here:
lam_min = 1e-14;
lam = max(lam, lam_min);

LogM = R * diag(log(lam)) * R';
v = [LogM(1,1), LogM(1,2), LogM(2,2)];
end

function [a,b,c] = expSPD2(v11,v12,v22)
% Compute exp([v11 v12; v12 v22]) and return components (a,b,c).
V = [v11 v12; v12 v22];
V = 0.5*(V+V');
[R,D] = eig(V);
lam = diag(D);
M = R * diag(exp(lam)) * R';
a = M(1,1);
b = M(1,2);
c = M(2,2);
end
