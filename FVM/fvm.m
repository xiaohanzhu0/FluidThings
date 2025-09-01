% poisson_fvm_nonuniform.m
% Solve -Lap(phi) = f on [xmin,xmax]×[ymin,ymax] with Dirichlet BC
% using a cell-centered finite‐volume scheme on a non-uniform grid.

clear; clc; close all;

%% 1. Define computational domain and nonuniform grid
xmin = 0; xmax = 1;
ymin = 0; ymax = 1;
Nx = 50; Ny = 40;

% Example non-uniform spacing (stretching toward boundaries)
beta = 1.5;  % stretching factor >1 clusters points near edges
xi = linspace(0,1,Nx+1).^beta;
yj = linspace(0,1,Ny+1).^beta;
xv = xmin + (xmax-xmin)*xi;  % cell faces in x
yv = ymin + (ymax-ymin)*yj;  % cell faces in y

% cell centers
xc = 0.5*(xv(1:end-1)+xv(2:end));
yc = 0.5*(yv(1:end-1)+yv(2:end));
[XX,YY] = meshgrid(xc,yc);    % (Ny × Nx) arrays

%% 2. Define RHS f(x,y) and Dirichlet BC phi_D
f = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);       % e.g. exact solution sin(pi x)sin(pi y)
phiD = @(x,y) sin(pi*x).*sin(pi*y);           % Dirichlet BC (matches exact)

%% 3. Precompute mesh metrics: Δx_w,e and Δy_s,n and cell areas
dx_w = zeros(Ny,Nx); dx_e = zeros(Ny,Nx);
dy_s = zeros(Ny,Nx); dy_n = zeros(Ny,Nx);
V   = zeros(Ny,Nx);

for j = 1:Ny
  for i = 1:Nx
    % distances to west/east faces
    dx_w(j,i) = xc(i) - xv(i);
    dx_e(j,i) = xv(i+1) - xc(i);
    % distances to south/north faces
    dy_s(j,i) = yc(j) - yv(j);
    dy_n(j,i) = yv(j+1) - yc(j);
    % cell area
    V(j,i) = (xv(i+1)-xv(i))*(yv(j+1)-yv(j));
  end
end

%% 4. Assemble sparse linear system A*phi = b
N = Nx*Ny;
I = []; J = []; S = [];    % for sparse matrix
b = zeros(N,1);

idx = @(i,j) (j-1)*Nx + i;  % linear index

for j = 1:Ny
  for i = 1:Nx
    P = idx(i,j);
    
    % initialize b
    bP = f(xc(i),yc(j)) * V(j,i);
    
    % West neighbor
    if i>1
      AW = 1/dx_w(j,i) * (yv(j+1)-yv(j));
      I(end+1)=P; J(end+1)=idx(i-1,j); S(end+1)=-AW;
      I(end+1)=P; J(end+1)=P;        S(end+1)= AW;
    else
      % Dirichlet at west boundary
      phiW = phiD( xv(1), yc(j) );
      AW = 1/dx_w(j,i) * (yv(j+1)-yv(j));
      bP = bP + AW * phiW;
      I(end+1)=P; J(end+1)=P; S(end+1)= AW;
    end
    
    % East neighbor
    if i<Nx
      AE = 1/dx_e(j,i) * (yv(j+1)-yv(j));
      I(end+1)=P; J(end+1)=idx(i+1,j); S(end+1)=-AE;
      I(end+1)=P; J(end+1)=P;        S(end+1)= AE;
    else
      phiE = phiD( xv(end), yc(j) );
      AE = 1/dx_e(j,i) * (yv(j+1)-yv(j));
      bP = bP + AE * phiE;
      I(end+1)=P; J(end+1)=P; S(end+1)= AE;
    end
    
    % South neighbor
    if j>1
      AS = 1/dy_s(j,i) * (xv(i+1)-xv(i));
      I(end+1)=P; J(end+1)=idx(i,j-1); S(end+1)=-AS;
      I(end+1)=P; J(end+1)=P;        S(end+1)= AS;
    else
      phiS = phiD( xc(i), yv(1) );
      AS = 1/dy_s(j,i) * (xv(i+1)-xv(i));
      bP = bP + AS * phiS;
      I(end+1)=P; J(end+1)=P; S(end+1)= AS;
    end
    
    % North neighbor
    if j<Ny
      AN = 1/dy_n(j,i) * (xv(i+1)-xv(i));
      I(end+1)=P; J(end+1)=idx(i,j+1); S(end+1)=-AN;
      I(end+1)=P; J(end+1)=P;        S(end+1)= AN;
    else
      phiN = phiD( xc(i), yv(end) );
      AN = 1/dy_n(j,i) * (xv(i+1)-xv(i));
      bP = bP + AN * phiN;
      I(end+1)=P; J(end+1)=P; S(end+1)= AN;
    end
    
    b(P) = bP;
  end
end

A = sparse(I,J,S,N,N);

%% 5. Solve linear system
phi_vec = A\b;

%% 6. Reshape and plot
Phi = reshape(phi_vec, Nx, Ny)';  % back to Ny×Nx

figure; 
surf(XX,YY,Phi,'EdgeColor','none');
xlabel('x'); ylabel('y'); zlabel('\phi');
title('FVM solution of 2D Poisson on nonuniform grid');
colorbar; view(2);  % top view

%% 7. Compute and display max error (if exact known)
Phi_ex = phiD(XX,YY);
err = max(abs(Phi_ex(:)-Phi(:)));
fprintf('Max pointwise error = %g\n', err);
