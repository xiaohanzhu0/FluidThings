clear
N         = 50;        % grid size (including boundaries)
lam       = 6;         % Bratu parameter
tolNewton = 1e-8;
tolKrylov = 1e-3;
maxNewton = 20;

% --- Build and Factor Preconditioner (The Discrete Laplacian A) ---
% We will explicitly form the sparse matrix A once to create the preconditioner.
e = ones(N^2, 1);
% spdiags creates a sparse matrix from its diagonals.
% This is the standard 5-point stencil matrix for the negative Laplacian.
A = spdiags([-e -e 4*e -e -e], [-N^2 -1 0 1 N^2], N^2, N^2);

% Now, compute the incomplete Cholesky factorization of A.
% The '0' means we use ichol(0), which has the same sparsity pattern as A.
% This is very fast and memory-efficient.
setup.type = 'ict';
setup.droptol = 1e-3; % Drop tolerance for incomplete factorization
L = ichol(A, setup);

U = solveBratuNewtonKrylov(N, lam, tolNewton, tolKrylov, maxNewton, L);

% plot
x = linspace(0,1,N);
surf(x, x, U, 'EdgeColor','none');
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
title('Bratu Solution via Newton–Krylov');
view(30,45);




function U = solveBratuNewtonKrylov(N, lam, tolNewton, tolKrylov, maxNewton, L)
% solveBratuNewtonKrylov  Matrix‐free Newton–Krylov for Bratu on [0,1]^2
%   U = solveBratuNewtonKrylov(N, lam, tolNewton, tolKrylov, maxNewton)
%   solves -Δu = lam*exp(u) with u=0 on the boundary, on an N×N grid
%   (including boundary). Returns U as an N×N array.

h = 1/(N-1);           % grid spacing
u = zeros(N*N,1);      % initial guess (including boundaries)

for k = 1:maxNewton
    F = residual(u, N, h, lam);
    resNorm = norm(F);
    fprintf('Newton iter %2d: ||F|| = %.3e\n', k, resNorm);
    if resNorm < tolNewton
        break;
    end
    
    % finite‐difference step for J·v
    epsFD = sqrt(eps) * max(1, norm(u)) * 100000;
    
    
    % matrix‐free Jacobian action
    Afun = @(v) jtv(u, v, N, h, lam, epsFD);
    
    restart = 10;                     % GMRES restart length
    [delta, flag] = gmres(Afun, -F, restart, tolKrylov, 100, L, L');
    %if flag ~= 0
    %    error('GMRES failed to converge (flag=%d)', flag);
    %end
    
    u = u + delta;
end

U = reshape(u, N, N);
end


function F = residual(u, N, h, lam)
    % Compute F(u) = -Δ_h(u) - lam*exp(u) on interior; 
    % Dirichlet BCs are zero so we leave boundary residual = 0.
    U = reshape(u, N, N);
    Fmat = zeros(N, N);
    for i = 2:N-1
        for j = 2:N-1
            lap = (U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) ...
                   - 4*U(i,j)) / h^2;
            Fmat(i,j) = -lap - lam * exp(U(i,j));
        end
    end
    F = Fmat(:);
end


function jv = jtv(u, v, N, h, lam, epsFD)
    % Finite‐difference approximation of J(u)*v
    F0 = residual(u,       N, h, lam);
    F1 = residual(u + epsFD*v, N, h, lam);
    jv = (F1 - F0) / epsFD;
end


