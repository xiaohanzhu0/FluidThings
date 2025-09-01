% Structured mesh (could be curvilinear)
Nx = 11; Ny = 11;
[Xs,Ys] = ndgrid(linspace(0,1,Nx), linspace(0,1,Ny));

% (Example) mild curvilinear map: bow outward near top
x = Xs;
y = Ys + 0.15*Xs.*(1-Xs).*Ys.^2;

% Coefficient: anisotropic tensor varying with y
Kfun = @(x,y) [1+0.5*y, 0; 0, 0.2 + 0.8*y.^2];

% Load
ffun = @(x,y) 1.0;

% Dirichlet boundary
gfun = @(x,y) 0.0;

[u,A,b] = fem_q1_structured_quad(x,y,Kfun,ffun,gfun);
%%
% visualize
figure
surf(x, y, u, 'EdgeColor','none'); view(2); colorbar; axis equal tight
title('Q1 FEM on structured curvilinear quads');

figure
pcolor(x, y, u); view(2); colorbar; axis equal tight
title('Q1 FEM on structured curvilinear quads');

%%





function [u,A,b] = fem_q1_structured_quad(x, y, Kfun, ffun, gfun)
% Solve -div(K grad u) = f on a structured quadrilateral mesh
% x,y : (Nx x Ny) nodal coordinates (could be curvilinear)
% Kfun: @(x,y) -> 2x2 SPD tensor or scalar k (isotropic)
% ffun: @(x,y) -> scalar load
% gfun: @(x,y) -> Dirichlet value on boundary nodes
%
% Returns:
%   u : (Nx x Ny) solution
%   A : (N x N) stiffness (sparse)
%   b : (N x 1) load

[Nx,Ny] = size(x);
N  = Nx*Ny;
id = @(i,j) i + (j-1)*Nx;

% --- Gauss data (2x2) ---
gp = [-1/sqrt(3),  1/sqrt(3)];
w  = [1, 1];

% Precompute reference shape and derivatives at 4 Gauss points
% Node order: [1: (−1,−1), 2:(+1,−1), 3:(+1,+1), 4:(−1,+1)] on reference
Ng  = cell(2,2); dNdxi = cell(2,2); dNdeta = cell(2,2);
for a = 1:2
  xi = gp(a);
  for b = 1:2
    eta = gp(b);
    Nloc = 0.25 * [(1-xi)*(1-eta);
                   (1+xi)*(1-eta);
                   (1+xi)*(1+eta);
                   (1-xi)*(1+eta)];
    dN_dxi  = 0.25 * [-(1-eta); (1-eta); (1+eta); -(1+eta)];
    dN_deta = 0.25 * [-(1-xi); -(1+xi); (1+xi);  (1-xi)];
    Ng{a,b}     = Nloc;
    dNdxi{a,b}  = dN_dxi;
    dNdeta{a,b} = dN_deta;
  end
end

% --- Preallocate sparse triplets (16 entries per element) ---
ne = (Nx-1)*(Ny-1);
I = zeros(16*ne,1);
J = zeros(16*ne,1);
V = zeros(16*ne,1);
b = zeros(N,1);

% --- Loop over elements using implicit connectivity ---
ptr = 1;
for j = 1:Ny-1
  for i = 1:Nx-1
    % Element node ids (counter-clockwise)
    nd = [ id(i,  j);
           id(i+1,j);
           id(i+1,j+1);
           id(i,  j+1) ];

    % Element coordinates
    Xe = [ x(i,  j); x(i+1,j); x(i+1,j+1); x(i,  j+1) ];
    Ye = [ y(i,  j); y(i+1,j); y(i+1,j+1); y(i,  j+1) ];

    Ke = zeros(4,4);
    be = zeros(4,1);

    % 2x2 Gauss integration
    for a = 1:2
      for bq = 1:2
        Nloc    = Ng{a,bq};        % 4x1
        dN_dxi  = dNdxi{a,bq};     % 4x1
        dN_deta = dNdeta{a,bq};    % 4x1

        % Jacobian J = [x_xi x_eta; y_xi y_eta]
        x_xi  = Xe.'*dN_dxi;   x_eta  = Xe.'*dN_deta;
        y_xi  = Ye.'*dN_dxi;   y_eta  = Ye.'*dN_deta;
        Jm = [x_xi, x_eta; y_xi, y_eta];
        detJ = Jm(1,1)*Jm(2,2) - Jm(1,2)*Jm(2,1);
        invJ = (1/detJ) * [ Jm(2,2), -Jm(1,2); -Jm(2,1), Jm(1,1) ];

        % Gradients in physical coords: grad N = J^{-T} * [dN/dxi; dN/deta]
        dN_ref = [dN_dxi.'; dN_deta.'];  % 2x4
        gradN  = invJ.' * dN_ref;        % 2x4

        % Physical GP for coefficient/load
        xgp = Nloc.'*Xe;
        ygp = Nloc.'*Ye;

        Kxy = Kfun(xgp, ygp);
        if isscalar(Kxy), Kxy = Kxy*eye(2); end

        wgt = w(a)*w(bq)*detJ;

        Ke = Ke + (gradN.' * Kxy * gradN) * wgt;
        be = be + Nloc * (ffun(xgp, ygp) * wgt);
      end
    end

    % Scatter-add into global sparse triplets
    [RR,CC] = ndgrid(nd, nd);
    rng = ptr:(ptr+15);
    I(rng) = RR(:);
    J(rng) = CC(:);
    V(rng) = Ke(:);
    ptr = ptr + 16;

    b(nd) = b(nd) + be;
  end
end

A = sparse(I,J,V,N,N);

% --- Dirichlet BC on the outer boundary (strong) ---
isB = false(N,1); uB = zeros(N,1);
for j = 1:Ny
  for i = 1:Nx
    onB = (i==1) || (i==Nx) || (j==1) || (j==Ny);
    if onB
      k = id(i,j);
      isB(k) = true;
      uB(k) = gfun(x(i,j), y(i,j));
    end
  end
end

free = ~isB;
b = b - A(:,isB)*uB(isB);
A(isB,:) = 0; A(:,isB) = 0;
A = A + spdiags(isB,0,N,N);     % put 1 on boundary diagonals
b(isB) = uB(isB);

% --- Solve and reshape ---
u = A \ b;
u = reshape(u, Nx, Ny);
end
