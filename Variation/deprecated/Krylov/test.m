clear
% compare_poisson_solvers.m
% Compare direct solve and GMRES performance for 2D Poisson equation using finite differences

% Grid parameters
n = 10;                 % number of interior points per dimension
h = 1/(n+1);
N = n^2;

% Construct 2D Laplacian with 5-point stencil
e = ones(n,1);
T = spdiags([-e 2*e -e], -1:1, n, n);
I = speye(n);
A = (kron(I,T) + kron(T,I)) / h^2;

% Right-hand side (f = 1 everywhere)
b = ones(N,1);

%load('GT01R.mat');
%load('venkat01.mat');
%load('RM07R.mat');

%A = Problem.A; 
%b = Problem.b;


% Direct solver
%tic;
%u_direct = A \ b;
%t_direct = toc;
%fprintf('Direct solve time: %.4f s\n', t_direct);

% GMRES solver
tol = 1e-6;
maxit = 500;
options = struct("type","ilutp","droptol",1e-2);
tic;
[L, U] = ilu(A, options);

[u_gmres, flag, relres, iter, resvec] = gmres(A, b, [], tol, maxit, L, U);
t_gmres = toc;
fprintf('GMRES time: %.4f s\n', t_gmres);
fprintf('GMRES flag: %d, relres: %.2e, iterations: %d\n', flag, relres, iter(2));

% Compute relative difference between solutions
diff_norm = norm(u_direct - u_gmres) / norm(u_direct);
fprintf('Relative difference: %.2e\n', diff_norm);


