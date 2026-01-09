clear
problem = 5;
Nx1 = 201;
Nx2 = 201;
[x1,x2,M11_fun,M12_fun,M22_fun] = problems.Initialization(problem,Nx1,Nx2);
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
axis equal
%%
%{
clear
addpath('/Users/zhuxiaohan/Documents/repo/Fluid/Variation/utils')
addpath('/Users/zhuxiaohan/Documents/repo/Fluid/Variation/test/Harmonic')
problem = 4;
orthongonal_project = 1;

Nx1 = 51;
Nx2 = 101;

[x1,x2,M11_fun,M12_fun,M22_fun] = problems.Initialization(problem,Nx1,Nx2);


boundary_points = extract_boundary_points(x1,x2);

figure
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');

[x1, x2, info] = solve_harmonic(x1, x2, Omega=0.5, ShowPlot=true, ...
                                BoundaryPoints=boundary_points, PauseTime=0.1);
%}
%%
[Nx2,Nx1] = size(x1);
ds1 = 1 / (Nx1-1);
ds2 = 1 / (Nx2-1);
[M11, M12, M22] = metric2meshv2(x1, x2, ds1, ds2);

M = metric_component_assembly(M11, M12, M22);
X = zeros(2, Nx2, Nx1);
X(1,:,:) = reshape(x1, [1 Nx2 Nx1]);
X(2,:,:) = reshape(x2, [1 Nx2 Nx1]);

figure
plot(x1, x2, 'r'); hold on; plot(x1', x2', 'r'); hold on
h = plotMetricEllipses(M, X);

M11_fun = scatteredInterpolant(x1(:),x2(:),M11(:));
M12_fun = scatteredInterpolant(x1(:),x2(:),M12(:));
M22_fun = scatteredInterpolant(x1(:),x2(:),M22(:));
x1_A = x1;
x2_A = x2;
%%
boundary_points = extract_boundary_points(x1,x2);
[x1, x2, info] = solve_harmonic(x1, x2, Omega=0.5, ShowPlot=true, ...
                                BoundaryPoints=boundary_points, PauseTime=0.1);
%%
[Nx2, Nx1] = size(x1);
% Interpolant of X(T)
[T1,T2] = ndgrid(linspace(0,1,Nx1),linspace(0,1,Nx2));
x1_in_t = griddedInterpolant(T1, T2, x1', 'cubic');
x2_in_t = griddedInterpolant(T1, T2, x2', 'cubic');
% Interpolant of J = dX(t)/dT
[dx1dt1_samp, dx2dt2_samp] = DCentral(x1, x2, 1/(Nx1-1), 1/(Nx2-1));
[dx2dt1_samp, dx1dt2_samp] = DCentral(x2, x1, 1/(Nx1-1), 1/(Nx2-1));


[T1,T2] = meshgrid(linspace(0,1,Nx1),linspace(0,1,Nx2));
t1 = linspace(0,1,313);
t2 = linspace(0,1,201);
[t1, t2] = meshgrid(t1, t2);
dx1dt1 = interp2(T1,T2,dx1dt1_samp,t1,t2, 'cubic');
dx1dt2 = interp2(T1,T2,dx1dt2_samp,t1,t2, 'cubic');
dx2dt1 = interp2(T1,T2,dx2dt1_samp,t1,t2, 'cubic');
dx2dt2 = interp2(T1,T2,dx2dt2_samp,t1,t2, 'cubic');

x1_temp = x1_in_t(t1,t2);
x2_temp = x2_in_t(t1,t2);
M11 = M11_fun(x1_temp, x2_temp);
M12 = M12_fun(x1_temp, x2_temp);
M22 = M22_fun(x1_temp, x2_temp);

Mp11 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt1 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt1;
Mp22 = (dx1dt2.*M11+dx2dt2.*M12).*dx1dt2 + (dx1dt2.*M12+dx2dt2.*M22).*dx2dt2;
Mp12 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt2 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt2;
clear x1_temp x2_temp M11 M12 M22

Mp11_fun = griddedInterpolant(t1', t2', Mp11', 'cubic');
Mp12_fun = griddedInterpolant(t1', t2', Mp12', 'cubic');
Mp22_fun = griddedInterpolant(t1', t2', Mp22', 'cubic');
% Solve inverse instead
detJ = Mp11.*Mp22 - Mp12.^2;
Mp_inv11 = Mp22 ./ detJ;
Mp_inv22 = Mp11 ./ detJ;
Mp_inv12 = -Mp12 ./ detJ;

Nt1 = 313; Nt2 = 201;
N = Nt1*Nt2;
s1 = linspace(0,1,Nt1);
s2 = linspace(0,1,Nt2);
[s1, s2] = meshgrid(s1, s2);
t1 = linspace(0,1,Nt1);
t2 = linspace(0,1,Nt2);
[t1, t2] = meshgrid(t1, t2);

max_iter = 20;
Lt_list1 = zeros(1,max_iter);
Lx_list1 = zeros(1,max_iter);
Theta_1_list1 = zeros(1,max_iter);
Theta_2_list1 = zeros(1,max_iter);
Theta_inf_list1 = zeros(1,max_iter);

for i = 1:max_iter
opt.ortho = false;
opt.ortho_.lambda = 0.001;

[ds1dt1, ds2dt2] = DCentral(s1, s2, 1/(Nt1-1), 1/(Nt2-1));
[ds2dt1, ds1dt2] = DCentral(s2, s1, 1/(Nt1-1), 1/(Nt2-1));
J = abs(ds1dt1.*ds2dt2 - ds2dt1.*ds1dt2);

Mp_inv11_fun = griddedInterpolant(t1', t2', Mp_inv11'.*J', 'cubic');
Mp_inv22_fun = griddedInterpolant(t1', t2', Mp_inv22'.*J', 'cubic');
Mp_inv12_fun = griddedInterpolant(t1', t2', Mp_inv12'.*J', 'cubic');

%
[s1_new, s2_new] = AssembleLinearSystemDual(t1, t2, Mp_inv11_fun, Mp_inv12_fun, Mp_inv22_fun);

s1 = s1 + 0.5*(s1_new-s1);
s2 = s2 + 0.5*(s2_new-s2);

s1_samp = linspace(0,1,313);
s2_samp = linspace(0,1,201);
[s1_samp,s2_samp] = meshgrid(s1_samp,s2_samp);
t1_samp = griddata(s1,s2,t1,s1_samp,s2_samp);
t2_samp = griddata(s1,s2,t2,s1_samp,s2_samp);
x1_samp = x1_in_t(t1_samp,t2_samp);
x2_samp = x2_in_t(t1_samp,t2_samp);

figure(8)
plot(s1, s2, 'k'); hold on; plot(s1', s2', 'k'); hold off
figure(9)
plot(t1_samp, t2_samp, 'k'); hold on; plot(t1_samp', t2_samp', 'k'); hold off
figure(10)
plot(x1_samp, x2_samp, 'k'); hold on; plot(x1_samp', x2_samp', 'k'); hold off
axis equal
pause(0.1)

M11_samp = Mp11_fun(t1_samp, t2_samp);
M12_samp = Mp11_fun(t1_samp, t2_samp);
M22_samp = Mp11_fun(t1_samp, t2_samp);
[Lt1, Lt2] = Cost(t1_samp, t2_samp, M11_samp, M12_samp, M22_samp);
Lt_list1(i) = mean(Lt1+Lt2,'all');

M11_samp = M11_fun(x1_samp, x2_samp);
M12_samp = M12_fun(x1_samp, x2_samp);
M22_samp = M22_fun(x1_samp, x2_samp);
[Lx1, Lx2] = Cost(x1_samp, x2_samp, M11_samp, M12_samp, M22_samp);
Lx_list1(i) = mean(Lx1+Lx2,'all');

[Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1_samp, x2_samp);
Theta_1_list1(i) = Theta_1;
Theta_2_list1(i) = Theta_2;
Theta_inf_list1(i) = Theta_inf;
end
figure()
plot(x1_A, x2_A, 'r'); hold on; plot(x1_A', x2_A', 'r'); 
plot(x1_samp, x2_samp, 'k'); hold on; plot(x1_samp', x2_samp', 'k'); hold off
%%
[M11, M12, M22] = metric2meshv2(x1, x2, ds1, ds2);
[M11_samp, M12_samp, M22_samp] = metric2meshv2(x1_samp, x2_samp, ds1, ds2);
figure
surf(x1_A, x2_A,(M11-M11_samp)./M11)
%%
%{
ds1 = 1 / (Nx1-1);
ds2 = 1 / (Nx1-1);
[g11, g12, g22] = mesh2metric(x1, x2, ds1, ds2);
figure
pcolor(x1, x2, g11);
M = metric_component_assembly(g11, g12, g22);
X = zeros(2, Nx2, Nx1);
X(1,:,:) = reshape(x1, [1 Nx2 Nx1]);
X(2,:,:) = reshape(x2, [1 Nx2 Nx1]);

h = plotMetricEllipses(M, X);
%%
J = g11.*g22 - g12.*g12;
g11_inv = g22 ./ J;
g22_inv = g11 ./ J;
g12_inv = -g12 ./ J;

figure
pcolor(x1, x2, g11_inv);
M = metric_component_assembly(g11_inv, g12_inv, g22_inv);
X = zeros(2, Nx2, Nx1);
X(1,:,:) = reshape(x1, [1 Nx2 Nx1]);
X(2,:,:) = reshape(x2, [1 Nx2 Nx1]);

h = plotMetricEllipses(M, X);
%%
s1 = linspace(0,1,Nx1);
s2 = linspace(0,1,Nx2);
[s1, s2] = meshgrid(s1,s2);
x1_samp = linspace(0,1,Nx1);
x2_samp = linspace(0,1,Nx2);
[x1_samp, x2_samp] = meshgrid(x1_samp,x2_samp);

s1_samp = griddata(x1,x2,s1,x1_samp,x2_samp);
s2_samp = griddata(x1,x2,s2,x1_samp,x2_samp);
plot(s1_samp, s2_samp, 'k'); hold on; plot(s1_samp', s2_samp', 'k');
%%
function [g11, g12, g22] = mesh2metric(x1, x2, ds1, ds2)

    [Nx, Ny] = size(x1);
    x1_s1  = zeros(Nx, Ny);
    x1_s2 = zeros(Nx, Ny);
    x2_s1  = zeros(Nx, Ny);
    x2_s2 = zeros(Nx, Ny);

    % Interior: central difference
    x1_s2(2:Nx-1, :) = (x1(3:Nx, :) - x1(1:Nx-2, :)) / (2 * ds2);
    x2_s2(2:Nx-1, :) = (x2(3:Nx, :) - x2(1:Nx-2, :)) / (2 * ds2);

    % Left boundary: forward difference
    x1_s2(1, :) = (x1(2, :) - x1(1, :)) / ds2;
    x2_s2(1, :) = (x2(2, :) - x2(1, :)) / ds2;

    % Right boundary: backward difference
    x1_s2(Nx, :) = (x1(Nx, :) - x1(Nx-1, :)) / ds2;
    x2_s2(Nx, :) = (x2(Nx, :) - x2(Nx-1, :)) / ds2;

    % Interior: central difference
    x1_s1(:, 2:Ny-1) = (x1(:, 3:Ny) - x1(:, 1:Ny-2)) / (2 * ds1);
    x2_s1(:, 2:Ny-1) = (x2(:, 3:Ny) - x2(:, 1:Ny-2)) / (2 * ds1);

    % Bottom boundary: forward difference
    x1_s1(:, 1) = (x1(:, 2) - x1(:, 1)) / ds1;
    x2_s1(:, 1) = (x2(:, 2) - x2(:, 1)) / ds1;

    % Top boundary: backward difference
    x1_s1(:, Ny) = (x1(:, Ny) - x1(:, Ny-1)) / ds1;
    x2_s1(:, Ny) = (x2(:, Ny) - x2(:, Ny-1)) / ds1;

    g11 = x1_s1.^2  + x2_s1.^2;
    g22 = x1_s2.^2 + x2_s2.^2;
    g12 = x1_s1 .* x1_s2 + x2_s1 .* x2_s2;
end
%}
function [M11, M12, M22] = metric2meshv2(x1, x2, ds1, ds2)
%METRIC_FROM_MESH Construct a physical-space metric M(x) from a known mesh.
%
%   [M11, M12, M22] = METRIC_FROM_MESH(x1, x2)
%   [M11, M12, M22] = METRIC_FROM_MESH(x1, x2, dxi, deta)
%
%   Given a mapping x_A : (xi,eta) -> (x1,x2) specified on a structured grid
%   in computational space (xi,eta), this function computes the metric tensor
%   M(x) in physical coordinates that makes this mesh locally metric-uniform:
%
%       M(x_A) = (J_A^{-1})^T * J_A^{-1},
%
%   where J_A = d(x1,x2)/d(xi,eta) is the Jacobian of the mapping.
%
%   Inputs
%   ------
%   x1, x2 : (Nx x Ny) arrays
%       Physical coordinates of the mesh nodes. First index is xi-direction,
%       second index is eta-direction.
%
%   dxi, deta : scalars (optional)
%       Computational grid spacing in xi and eta. Defaults: dxi = 1, deta = 1.
%
%   Outputs
%   -------
%   M11, M12, M22 : (Nx x Ny) arrays
%       Components of the SPD metric tensor M(x) in physical coordinates:
%           M = [M11  M12;
%                M12  M22].
%
%   Notes
%   -----
%   - Derivatives are computed with central differences in the interior and
%     one-sided differences on the boundaries.
%   - No checks are done on det(J); if it is very small or changes sign,
%     the resulting metric may be ill-conditioned or undefined.

    [Nx, Ny] = size(x1);
    x1_s1  = zeros(Nx, Ny);
    x1_s2 = zeros(Nx, Ny);
    x2_s1  = zeros(Nx, Ny);
    x2_s2 = zeros(Nx, Ny);

    % Interior: central difference
    x1_s2(2:Nx-1, :) = (x1(3:Nx, :) - x1(1:Nx-2, :)) / (2 * ds2);
    x2_s2(2:Nx-1, :) = (x2(3:Nx, :) - x2(1:Nx-2, :)) / (2 * ds2);

    % Left boundary: forward difference
    x1_s2(1, :) = (x1(2, :) - x1(1, :)) / ds2;
    x2_s2(1, :) = (x2(2, :) - x2(1, :)) / ds2;

    % Right boundary: backward difference
    x1_s2(Nx, :) = (x1(Nx, :) - x1(Nx-1, :)) / ds2;
    x2_s2(Nx, :) = (x2(Nx, :) - x2(Nx-1, :)) / ds2;

    % Interior: central difference
    x1_s1(:, 2:Ny-1) = (x1(:, 3:Ny) - x1(:, 1:Ny-2)) / (2 * ds1);
    x2_s1(:, 2:Ny-1) = (x2(:, 3:Ny) - x2(:, 1:Ny-2)) / (2 * ds1);

    % Bottom boundary: forward difference
    x1_s1(:, 1) = (x1(:, 2) - x1(:, 1)) / ds1;
    x2_s1(:, 1) = (x2(:, 2) - x2(:, 1)) / ds1;

    % Top boundary: backward difference
    x1_s1(:, Ny) = (x1(:, Ny) - x1(:, Ny-1)) / ds1;
    x2_s1(:, Ny) = (x2(:, Ny) - x2(:, Ny-1)) / ds1;

    %----------------------------------------------------------------------
    % 2. Jacobian inverse J^{-1} at each node
    %
    %   J = [x1_xi   x1_eta;
    %        x2_xi   x2_eta];
    %
    %   detJ = x1_xi * x2_eta - x1_eta * x2_xi
    %
    %   J^{-1} = (1/detJ) * [ x2_eta   -x1_eta;
    %                        -x2_xi     x1_xi];
    %----------------------------------------------------------------------

    detJ = x1_s1 .* x2_s2 - x1_s2 .* x2_s1;

    Jinv11 =  x2_s2 ./ detJ;
    Jinv12 = -x1_s2 ./ detJ;
    Jinv21 = -x2_s1 ./ detJ;
    Jinv22 =  x1_s1 ./ detJ;

    %----------------------------------------------------------------------
    % 3. Metric M = J^{-T} * J^{-1} in physical coordinates
    %
    %   Let Jinv = [a b; c d].
    %   Then M = Jinv^T * Jinv gives:
    %       M11 = a^2 + c^2
    %       M12 = a*b + c*d
    %       M22 = b^2 + d^2
    %----------------------------------------------------------------------

    M11 = Jinv11.^2 + Jinv21.^2;
    M12 = Jinv11 .* Jinv12 + Jinv21 .* Jinv22;
    M22 = Jinv12.^2 + Jinv22.^2;
end
