%% Transformed-grid plots for u(x,y), v(x,y) on [0,1]^2  --> (u,v)
% Boundary sets used (from our derivations):
% u: u(0,y) = 1-2y (y<1/2), 0 (y>=1/2),  u(1,y) = 1-y (y<1/2), 1/2 (y>=1/2),
%    u(x,0) = 1, u(x,1) = x/2
% v: v(0,y) = 0 (y<1/2), 2y-1 (y>=1/2), v(1,y) = 1/2 (y<1/2), y (y>=1/2),
%    v(x,0) = x/2, v(x,1) = 1

clear; clc;

% --- choose spectral truncation ---
Nmax = 199;   % only odd n contribute; keep this odd-ish and ~O(100-300)

% Build function handles that accept vectors/matrices X,Y
ufun = @(X,Y) u_series_eval(X,Y,Nmax);
vfun = @(X,Y) v_series_eval(X,Y,Nmax);

% === A) Parametric grid (images of x=const and y=const lines) ===
figure; plotTransformedGrid_parametric(ufun, vfun, 21, 21);

% === B) Wireframe (entire rect grid pushed forward) ===
figure; plotTransformedGrid_wireframe(ufun, vfun, 61, 61);

%% --------- helper: parametric grid ----------
function plotTransformedGrid_parametric(ufun, vfun, nx, ny)
if nargin<3, nx=21; ny=21; end
hold on; axis equal; box on; xlabel('u'); ylabel('v');
col = [0.6 0.6 0.6]; lw = 0.6;
tLong = linspace(0,1,600).';

% verticals x = const
for x0 = linspace(0,1,nx)
    X = x0*ones(size(tLong));
    plot(ufun(X,tLong), vfun(X,tLong), '-', 'Color', col, 'LineWidth', lw);
end
% horizontals y = const
for y0 = linspace(0,1,ny)
    Y = y0*ones(size(tLong));
    plot(ufun(tLong,Y), vfun(tLong,Y), '-', 'Color', col, 'LineWidth', lw);
end
% boundaries (thicker)
plotBoundary(ufun, vfun);
end

%% --------- helper: wireframe grid ----------
function plotTransformedGrid_wireframe(ufun, vfun, nx, ny)
if nargin<3, nx=61; ny=61; end
[xg,yg] = meshgrid(linspace(0,1,nx), linspace(0,1,ny));
U = ufun(xg,yg); V = vfun(xg,yg);
surf(U, V, zeros(size(U)), 'FaceColor','none', 'EdgeColor',[0.55 0.55 0.55]);
hold on; view(2); axis equal tight; box on; xlabel('u'); ylabel('v');
plotBoundary(ufun, vfun);
end

%% --------- helper: draw boundary images ----------
function plotBoundary(ufun, vfun)
t = linspace(0,1,1000).';
plot(ufun(0*t,t),   vfun(0*t,t),   'k-', 'LineWidth',1.5); % x=0
plot(ufun(1+0*t,t), vfun(1+0*t,t), 'k-', 'LineWidth',1.5); % x=1
plot(ufun(t,0*t),   vfun(t,0*t),   'k-', 'LineWidth',1.5); % y=0
plot(ufun(t,1+0*t), vfun(t,1+0*t), 'k-', 'LineWidth',1.5); % y=1
end

%% --------- u-series evaluator (vectorized, stable) ----------
function U = u_series_eval(X,Y,Nmax)
% Particular (harmonic) matching u(x,0)=1, u(x,1)=x/2:
U = 1 - Y + (X.*Y)/2;

% Only odd n contribute
n = (1:2:Nmax).';           % Kx1
N = n*pi;                   % Kx1
sn = sin(N/2);              % = +/-1 for odd n
cn = -4*sn./(N.^2);         % c_n
dn = 0.5*cn;                % d_n = c_n/2

% Stable ratios for sinh terms:
% r1 = sinh(N*(1-x))/sinh(N) = exp(-N*x) * (1 - exp(-2N(1-x))) / (1 - exp(-2N))
% r2 = sinh(N*x)/sinh(N)     = exp(-N*(1-x)) * (1 - exp(-2N*x)) / (1 - exp(-2N))
Xv = X(:).';  Yv = Y(:).';
E  = exp(-N);               % Kx1
den = 1 - E.^2;             % Kx1
r1 = exp(-N.*Xv)        .* (1 - exp(-2*N.*(1 - Xv))) ./ den;  % KxM
r2 = exp(-N.*(1 - Xv)) .* (1 - exp(-2*N.*Xv))         ./ den;  % KxM

S  = sin(N.*Yv);            % KxM
add = sum((cn.*r1 + dn.*r2).*S, 1);   % 1xM
U   = U + reshape(add, size(X));
end

%% --------- v-series evaluator (vectorized, stable) ----------
function V = v_series_eval(X,Y,Nmax)
% Particular (harmonic) matching v(x,0)=x/2, v(x,1)=1:
V = Y + X/2 - (X.*Y)/2;

% Coefficients identical pattern (g_R = 1/2 g_L) so same cn,dn as u-case
n = (1:2:Nmax).';
N = n*pi;
sn = sin(N/2);
cn = -4*sn./(N.^2);
dn = 0.5*cn;

Xv = X(:).';  Yv = Y(:).';
E  = exp(-N);
den = 1 - E.^2;
r1 = exp(-N.*Xv)        .* (1 - exp(-2*N.*(1 - Xv))) ./ den;
r2 = exp(-N.*(1 - Xv)) .* (1 - exp(-2*N.*Xv))         ./ den;

S  = sin(N.*Yv);
add = sum((cn.*r1 + dn.*r2).*S, 1);
V   = V + reshape(add, size(X));
end
