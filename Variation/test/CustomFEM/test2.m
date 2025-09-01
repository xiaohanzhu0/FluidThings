clear
% Example mesh (structured, possibly curved)
Nx = 81; Ny = 81;
problem = 2;

if problem == 1
    [xc, yc] = meshgrid(linspace(0,1,Nx), linspace(0,1,Ny));
    amp = 0.15;
    X = xc + amp*sin(pi*yc).* (xc.*(1-xc));
    Y = yc + amp*sin(pi*xc).* (yc.*(1-yc));
    Mfun = @(x,y) [ 2 + x.^2,        0.5*sin(pi*x).*cos(pi*y);
                    0.5*sin(pi*x).*cos(pi*y), 1 + y.^2 ];
elseif problem == 2
    [X, Y] = InitProb8(Nx, Ny);
    Mfun = @(x,y) [1,0;0,1];
end

[S1, S2, A] = fem_Q1_structured_divMgrad(X, Y, Mfun);
%%

% s1=0 on left, 1 on right; s2=0 on bottom, 1 on top
figure; contour(X, Y, S1, 2*Nx, 'k'); axis equal tight; title('s_1'); 
hold on; contour(X, Y, S2, 2*Ny, 'k'); axis equal tight; title('s_2'); 


%%
Nx = 81; Ny = 81;
%if problem == 1
%    [xc, yc] = meshgrid(linspace(0,1,Nx), linspace(0,1,Ny));
%    Xq = xc + amp*sin(pi*yc).* (xc.*(1-xc));          % interior warp in x
%    Yq = yc + amp*sin(pi*xc).* (yc.*(1-yc));          % interior warp in y
%elseif problem ==2
%    [Xq, Yq] = InitProb8(Nx, Ny);
%end

%S1q = griddata(X,Y,S1,Xq,Yq);
%S2q = griddata(X,Y,S2,Xq,Yq);
S1q = linspace(0,1,Nx);
S2q = linspace(0,1,Ny);
[S1q, S2q] = meshgrid(S1q, S2q);

Xq = griddata(S1,S2,X,S1q,S2q,"cubic");
Yq = griddata(S1,S2,Y,S1q,S2q,"cubic");

figure
plot(Xq, Yq, 'k'); hold on; plot(Xq', Yq', 'k');


%%
function [S1, S2, A] = fem_Q1_structured_divMgrad(X, Y, Mfun)
% FEM solver for ∇·(M(x) ∇ s) = 0 with Q1 (bilinear) elements on a structured quad mesh.
% - X, Y: Ny-by-Nx arrays of node coordinates (like meshgrid output; curved boundaries OK)
% - Mfun: function handle @(x,y) returning a 2x2 SPD matrix M(x,y)
% Returns:
% - S1, S2: Ny-by-Nx arrays for the two components with boundary data described below
% - A: assembled global stiffness (symmetric sparse)
%
% Boundary conditions (by component):
%   s1: Dirichlet s1=0 on Left (i=1), s1=1 on Right (i=Nx); natural Neumann on Bottom/Top
%   s2: Dirichlet s2=0 on Bottom (j=1), s2=1 on Top (j=Ny);  natural Neumann on Left/Right

  [Ny, Nx] = size(X);
  N = Nx*Ny;
  idx = @(i,j) (j-1)*Nx + i;           % linear index for node (i,j)

  %----- Assembly over structured cells -----------------------------------
  ne = (Nx-1)*(Ny-1);                  % number of elements
  I = zeros(16*ne,1); J = I; V = I;    % triplets for sparse A
  ptr = 0;

  % 2x2 Gauss points on reference square [-1,1]^2
  gp = [-1, +1]/sqrt(3);  wgt = [1, 1];

  for j = 1:Ny-1
    for i = 1:Nx-1
      % 4 local nodes (counterclockwise): (i,j)->(i+1,j)->(i+1,j+1)->(i,j+1)
      en = [ idx(i,j), idx(i+1,j), idx(i+1,j+1), idx(i,j+1) ];
      Xe = [ X(j,i);   X(j,  i+1);  X(j+1,i+1);  X(j+1,i)   ];
      Ye = [ Y(j,i);   Y(j,  i+1);  Y(j+1,i+1);  Y(j+1,i)   ];

      Ke = zeros(4,4);

      % Gauss integration
      for a = 1:2
        xi = gp(a);
        for b = 1:2
          eta = gp(b);

          % Q1 shape functions and derivatives on reference element
          Nhat  = 0.25 * [(1-xi)*(1-eta); (1+xi)*(1-eta); (1+xi)*(1+eta); (1-xi)*(1+eta)];
          dNdxi = 0.25 * [ -(1-eta); +(1-eta); +(1+eta); -(1+eta) ];
          dNdeta= 0.25 * [ -(1-xi);  -(1+xi);  +(1+xi);  +(1-xi)  ];

          % Isoparametric map x(ξ,η) = Σ N_k Xk, y(ξ,η) = Σ N_k Yk
          xg = Xe.'*Nhat;  yg = Ye.'*Nhat;

          % Jacobian of the mapping J = [x_ξ x_η; y_ξ y_η]
          x_xi  = Xe.'*dNdxi;   x_eta  = Xe.'*dNdeta;
          y_xi  = Ye.'*dNdxi;   y_eta  = Ye.'*dNdeta;
          Jm = [x_xi, x_eta; y_xi, y_eta];
          detJ = abs(x_xi*y_eta - x_eta*y_xi);
          if detJ <= 0
            error('Element (%d,%d) has nonpositive Jacobian (folded or wrong ordering).', i, j);
          end

          % Gradients of shape functions in physical coords: ∇N = J^{-T} ∇̂N
          % Build 2x4 reference gradient matrix, then apply J^{-T}
          Gref = [dNdxi.'; dNdeta.'];         % size 2x4
          G = (Jm.' \ Gref);                  % 2x4, columns are ∇N_k

          % Material tensor at Gauss point
          M = Mfun(xg, yg);                   % 2x2 SPD
          % Element stiffness contribution: K_e += (∇N)^T M (∇N) detJ w
          Ke = Ke + (G.' * (M * G)) * detJ * wgt(a) * wgt(b);
        end
      end

      % Scatter-add into global triplets
      [ii, jj] = ndgrid(en, en);
      nn = 16;
      I(ptr+(1:nn)) = ii(:);
      J(ptr+(1:nn)) = jj(:);
      V(ptr+(1:nn)) = Ke(:);
      ptr = ptr + nn;
    end
  end

  % Global stiffness (symmetrize for numerical robustness)
  A = sparse(I, J, V, N, N);
  A = 0.5*(A + A.');

  %----- Solve for s1: Dirichlet on Left/Right -----------------------------
  left   = (1:Nx:Nx*Ny).';             % i=1, all j
  right  = (Nx:Nx:Nx*Ny).';            % i=Nx, all j
  D1 = [left; right];
  g1 = [zeros(size(left)); ones(size(right))];

  b = zeros(N,1);
  [A1, b1] = impose_dirichlet(A, b, D1, g1);
  s1 = A1 \ b1;
  S1 = reshape(s1, [Nx, Ny]).';        % back to Ny-by-Nx

  %----- Solve for s2: Dirichlet on Bottom/Top -----------------------------
  bottom = (1:Nx).';                   % j=1, all i
  top    = ((Ny-1)*Nx + (1:Nx)).';     % j=Ny, all i
  D2 = [bottom; top];
  g2 = [zeros(size(bottom)); ones(size(top))];

  [A2, b2] = impose_dirichlet(A, b, D2, g2);
  s2 = A2 \ b2;
  S2 = reshape(s2, [Nx, Ny]).';        % Ny-by-Nx
end

%--------------------------- helpers ---------------------------------------
function [Amod, bmod] = impose_dirichlet(A, b, D, g)
% Strongly impose Dirichlet on nodes D with values g
  N = size(A,1);
  % Move known values to RHS first
  bmod = b - A(:, D) * g;
  % Zero Dirichlet rows/cols and put 1s on the diagonal
  Amod = A;
  Amod(:, D) = 0;
  Amod(D, :) = 0;
  Amod = Amod + sparse(D, D, 1, N, N);
  % Set Dirichlet RHS
  bmod(D) = g;
end
