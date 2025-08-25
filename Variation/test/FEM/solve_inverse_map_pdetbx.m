function [s1, s2, model1, model2] = solve_inverse_map_pdetbx(p, t, edgeArcs, Mfun, bc, Delta)
%SOLVE_INVERSE_MAP_PDETBX  Solve ∇·(Δ^2 M^{-1}(x)∇s)=0 with PDE Toolbox.
%
% Inputs
%   p        : (nP x 2) node coordinates [x y]
%   t        : (nT x 3) triangles (1-based indices into rows of p)
%   edgeArcs : struct with edge IDs for boundary arcs:
%              .g1, .g2, .g3, .g4  (vectors of edge IDs)
%              (Dirichlet: s1=0 on g1, s1=1 on g3; s2=0 on g2, s2=1 on g4.
%               All other edges are natural Neumann (zero flux).)
%              Optional ramped BCs:
%              .s1_u_on_g1, .s1_u_on_g3, .s2_u_on_g2, .s2_u_on_g4
%              each is either a scalar or a function handle @(loc,~)->values
%   Mfun     : function handle @(x,y) -> 2x2 SPD metric tensor M(x)
%   Delta    : [Delta1, Delta2] (defaults to [1,1] if omitted)
%
% Outputs
%   s1, s2   : (nP x 1) nodal solutions for s1 and s2 on the PDE mesh
%   model1/2 : PDEModel objects (contain Mesh, solution, etc.)
%
% Notes
%   • To see/choose edge IDs: pdegplot(model,'EdgeLabels','on').
%   • If you pass a mesh (p,t), we import it as the geometry; PDE Toolbox
%     will reuse it as the analysis mesh via geometryFromMesh.
%   • For strongly varying M, this remains linear—c(x) varies in space only.

    if nargin < 6 || isempty(Delta), Delta = [1,1]; end

    % --- s1: build model from your mesh and set coefficients
    model1 = createpde(1);
    geometryFromMesh(model1, p', t');           % <- uses your triangulation
    % (No need to call generateMesh; geometryFromMesh supplies a usable mesh.)
    c1 = @(loc,state) c_from_M(loc, Delta(1), Mfun);  % anisotropic c(x)

    specifyCoefficients(model1, 'm',0, 'd',0, 'c',c1, 'a',0, 'f',0);

    % Dirichlet on edges: s1=0 on g1, s1=1 on g3 (or use provided handles)
    u_g1 = getU(edgeArcs, 's1_u_on_g1', 0);
    u_g3 = getU(edgeArcs, 's1_u_on_g3', 1);
    %applyBoundaryCondition(model1,'dirichlet','Edge',edgeArcs.g1,'u',u_g1);
    %applyBoundaryCondition(model1,'dirichlet','Edge',edgeArcs.g3,'u',u_g3);
    applyBoundaryCondition(model1,'dirichlet','Edge',bc.x0,'u',u_g1);
    applyBoundaryCondition(model1,'dirichlet','Edge',bc.x1,'u',u_g3);
    % Other edges default to natural Neumann (zero flux): do nothing.

    % Solve and read nodal solution
    R1 = solvepde(model1);
    s1 = R1.NodalSolution;

    % --- s2: repeat with different edges and scaling
    model2 = createpde(1);
    geometryFromMesh(model2, p', t');
    c2 = @(loc,state) c_from_M(loc, Delta(2), Mfun);
    specifyCoefficients(model2, 'm',0, 'd',0, 'c',c2, 'a',0, 'f',0);

    u_g2 = getU(edgeArcs, 's2_u_on_g2', 0);
    u_g4 = getU(edgeArcs, 's2_u_on_g4', 1);
    %applyBoundaryCondition(model2,'dirichlet','Edge',edgeArcs.g2,'u',u_g2);
    %applyBoundaryCondition(model2,'dirichlet','Edge',edgeArcs.g4,'u',u_g4);
    applyBoundaryCondition(model2,'dirichlet','Edge',bc.y0,'u',u_g2);
    applyBoundaryCondition(model2,'dirichlet','Edge',bc.y1,'u',u_g4);

    R2 = solvepde(model2);
    s2 = R2.NodalSolution;

    % (Optional) sanity: ensure both models share the same node set
    % assert(isequal(model1.Mesh.Nodes, model2.Mesh.Nodes), 'Meshes differ.');

    % --- helper: vectorized anisotropic coefficient
    function c = c_from_M(location, Delta_, Mfun_)
        x = location.x; y = location.y; n = numel(x);
        c = zeros(4, n);                    % [c11; c12; c21; c22] for each point
        for k = 1:n
            Minv = Mfun_(x(k), y(k)) \ eye(2);  % robust inverse via solve
            C = (Delta_^2) * Minv;
            c(:,k) = [C(1,1); C(1,2); C(2,1); C(2,2)] * 1;%sqrt(det(Mfun_(x(k), y(k))));
        end
    end

    function u = getU(S, field, def)
        if isfield(S, field) && ~isempty(S.(field))
            u = S.(field);       % either scalar or @(loc,~) -> values
        else
            u = def;
        end
    end
end
