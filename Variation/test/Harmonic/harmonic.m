clear
problem = 5;
s1 = linspace(0, 1, 100);
s2 = linspace(0, 1, 100);
[x1, x2] = meshgrid(s1, s2);


if problem == 3
    [x1, x2] = InitProb3(100, 100, 0.1);
    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,1), x1(end,end); x2(end,1), x2(end,end)];
    boundary_points.l = [x1(1,1), x1(end,1); x2(1,1), x2(end,1)]';
    boundary_points.r = [x1(1,end), x1(end,end); x2(1,end), x2(end,end)]';
end
if problem == 4
    [x1, x2] = InitProb4(100, 100);
    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,1), x1(end,end); x2(end,1), x2(end,end)];
    boundary_points.l = [x1(1,1), x1(end,1); x2(1,1), x2(end,1)]';
    boundary_points.r = [x1(1,end), x1(end,end); x2(1,end), x2(end,end)]';
end
if problem == 5
    cf.new_airfoil = 1;
    cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
    cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
    cf.Nx1 = 201; cf.Nx2 = 51;
    cf.alpha = 1.005; cf.append_trail = 0;
    [x1, x2, M_samp] = InitProb5(cf);
    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,:); x2(end,:)];
    boundary_points.l = [x1(:,1), x2(:,1)];
    boundary_points.r = [x1(:,end), x2(:,end)];
    cf.Nx1 = size(x1,2);
    cf.Nx2 = size(x1,1);
    N = cf.Nx1 * cf.Nx2;
end

if problem == 8
    [x1, x2] = InitProb8(80, 40);
    x1 = flip(x1,2);
    x2 = flip(x2,2);
    cf.Nt1 = size(x1,2);
    cf.Nx2 = size(x1,1);
end
boundary_points.b = [x1(1,:); x2(1,:)];
boundary_points.t = [x1(end,:); x2(end,:)];
boundary_points.l = [x1(:,1), x2(:,1)];
boundary_points.r = [x1(:,end), x2(:,end)];

figure
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');


%%

[Nx2, Nx1] = size(x1);
N = Nx1*Nx2;
t1 = linspace(0,1,Nx1);
t2 = linspace(0,1,Nx2);
[t1, t2] = meshgrid(t1,t2);

dt1 = 1 / Nx1;
dt2 = 1 / Nx2;

for j=1:20
[dx1dt1, dx2dt2] = DCentral(x1, x2, dt1, dt2);
[dx2dt1, dx1dt2] = DCentral(x2, x1, dt1, dt2);
a = dx1dt2.^2 + dx2dt2.^2;
b = dx1dt1.*dx1dt2 + dx2dt1.*dx2dt2;
c = dx1dt1.^2 + dx2dt1.^2;

e = ones(N,1);
D2x = spdiags([e, -2*e, e], [-1, 0, 1], Nx1, Nx1);
D2y = spdiags([e, -2*e, e], [-1, 0, 1], Nx2, Nx2);

L = kron(D2x,speye(Nx2)).*a(:) / dt1^2 + kron(speye(Nx1),D2y).*c(:) / dt2^2;
L_mix = kron(spdiags([-e/2, e/2], [-1, 1], Nx1, Nx1), spdiags([-e/2, e/2], [-1, 1], Nx2, Nx2)).*b(:) / dt1 / dt2;
A = L - 2*L_mix;
A = blkdiag(A, A);


% Apply boundary conditions -----------------------------------------------
id = GetIndex(Nx1, Nx2);
t_bottom = GetBoundaryTangent(x1(1,:), x2(1,:), 1);
t_top = GetBoundaryTangent(x1(end,:), x2(end,:), 1);
t_left = GetBoundaryTangent(x1(:,1), x2(:,1), 1);
t_right = GetBoundaryTangent(x1(:,end), x2(:,end), 1);

n_bottom = [-t_bottom(2,:); t_bottom(1,:)];
n_top = [-t_top(2,:); t_top(1,:)];
n_left = [-t_left(2,:); t_left(1,:)];
n_right = [-t_right(2,:); t_right(1,:)];

I = [id.l, id.l, N+id.l, N+id.l, N+id.l, N+id.l, N+id.l, N+id.l, ...
     id.r, id.r, N+id.r, N+id.r, N+id.r, N+id.r, N+id.r, N+id.r, ...
     id.b, id.b, N+id.b, N+id.b, N+id.b, N+id.b, N+id.b, N+id.b, ...
     id.t, id.t, N+id.t, N+id.t, N+id.t, N+id.t, N+id.t, N+id.t];

J = [id.l, N+id.l, id.l, N+id.l, 1*Nx2+id.l, N+1*Nx2+id.l, 2*Nx2+id.l, N+2*Nx2+id.l, ...
     id.r, N+id.r, id.r, N+id.r, -1*Nx2+id.r, N-1*Nx2+id.r, -2*Nx2+id.r, N-2*Nx2+id.r, ...
     id.b, N+id.b, id.b, N+id.b, 1+id.b, N+1+id.b, 2+id.b, N+2+id.b, ...
     id.t, N+id.t, id.t, N+id.t, -1+id.t, N-1+id.t, -2+id.t, N-2+id.t];

V = [n_left(1,2:end-1), n_left(2,2:end-1), t_left(1,2:end-1), t_left(2,2:end-1), ...
     -4/3*t_left(1,2:end-1), -4/3*t_left(2,2:end-1), 1/3*t_left(1,2:end-1), 1/3*t_left(2,2:end-1), ...
     n_right(1,2:end-1), n_right(2,2:end-1), t_right(1,2:end-1), t_right(2,2:end-1), ...
     -4/3*t_right(1,2:end-1), -4/3*t_right(2,2:end-1), 1/3*t_right(1,2:end-1), 1/3*t_right(2,2:end-1), ...
     n_bottom(1,2:end-1), n_bottom(2,2:end-1), t_bottom(1,2:end-1), t_bottom(2,2:end-1), ...
     -4/3*t_bottom(1,2:end-1), -4/3*t_bottom(2,2:end-1), 1/3*t_bottom(1,2:end-1), 1/3*t_bottom(2,2:end-1), ...
     n_top(1,2:end-1), n_top(2,2:end-1), t_top(1,2:end-1), t_top(2,2:end-1), ...
     -4/3*t_top(1,2:end-1), -4/3*t_top(2,2:end-1), 1/3*t_top(1,2:end-1), 1/3*t_top(2,2:end-1)];

B = sparse(I,J,V,size(A,1),size(A,2));
A(id.l,:) = 0;
A(id.l+N,:) = 0;
A(id.r,:) = 0;
A(id.r+N,:) = 0;
A(id.b,:) = 0;
A(id.b+N,:) = 0;
A(id.t,:) = 0;
A(id.t+N,:) = 0;
A = A + B;

for i = id.corner; A(i, :) = 0; A(i, i) = 1; end
for i = N+id.corner; A(i, :) = 0; A(i, i) = 1; end

b1 = zeros(Nx2, Nx1);
b2 = zeros(Nx2, Nx1);
b1(:,1) = x1(:,1).*n_left(1,:)' + x2(:,1).*n_left(2,:)';
b1(:,end) = x1(:,end).*n_right(1,:)' + x2(:,end).*n_right(2,:)';
b1(1,:) = x1(1,:).*n_bottom(1,:) + x2(1,:).*n_bottom(2,:);
b1(end,:) = x1(end,:).*n_top(1,:) + x2(end,:).*n_top(2,:);
b2(:,1) = 0; b2(:,end) = 0; b2(1,:) = 0; b2(end,:) = 0;

b1(1,1) = x1(1,1); b1(end,1) = x1(end,1); b1(1,end) = x1(1,end); b1(end,end) = x1(end,end);
b2(1,1) = x2(1,1); b2(end,1) = x2(end,1); b2(1,end) = x2(1,end); b2(end,end) = x2(end,end);

% -------------------------------------------------------------------------

x_new = A \ [b1(:); b2(:)];
x1_new = reshape(x_new(1:N), Nx2, Nx1);
x2_new = reshape(x_new(N+1:end), Nx2, Nx1);
x1 = x1 + 0.5*(x1_new - x1);
x2 = x2 + 0.5*(x2_new - x2);

[x1, x2] = UpdateCorrection(x1, x2, boundary_points);

figure(2)
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k'); hold off
pause(0.1)
end
%%
% Interpolant of X(T)
[T1,T2] = ndgrid(linspace(0,1,Nx2),linspace(0,1,Nx1));
x1_in_t = griddedInterpolant(T1, T2, x1, 'linear');
x2_in_t = griddedInterpolant(T1, T2, x2, 'linear');


% Interpolant of J = dX(t)/dT
[dx1dt1_samp, dx2dt2_samp] = DCentral(x1, x2, dt1, dt2);
[dx2dt1_samp, dx1dt2_samp] = DCentral(x2, x1, dt1, dt2);

dx1dt1 = griddedInterpolant(T1, T2, dx1dt1_samp, 'linear');
dx1dt2 = griddedInterpolant(T1, T2, dx1dt2_samp, 'linear');
dx2dt1 = griddedInterpolant(T1, T2, dx2dt1_samp, 'linear');
dx2dt2 = griddedInterpolant(T1, T2, dx2dt2_samp, 'linear');

J_samp = dx1dt1_samp.*dx2dt2_samp - dx2dt1_samp.*dx1dt2_samp;

%%

Nx1 = 101; Nx2 = 101;
N = Nx1*Nx2;
S1 = linspace(0,1,Nx1);
S2 = linspace(0,1,Nx2);
t1 = linspace(0,1,Nx1);
t2 = linspace(0,1,Nx2);
[t1, t2] = meshgrid(t1, t2);


M11 = @(x1,x2) ones(size(x1));
M12 = @(x1,x2) zeros(size(x1));
M22 = @(x1,x2) ones(size(x1));

%M11 = @(t1,t2) 40000*(1+15*t1).^(-2);
%M12 = @(t1,t2) zeros(size(t1));
%M22 = @(t1,t2) ones(size(t1));


M11 = @(t1,t2) 1000 - 600*cos(pi*t1).*cos(pi*t2);
M12 = @(t1,t2) zeros(size(t1));
M22 = @(t1,t2) 1000 + 600*cos(pi*t1).*cos(pi*t2);

if problem == 5
    M11 = @(x1,x2) M_samp.F11(x1,x2);
    M12 = @(x1,x2) M_samp.F12(x1,x2);
    M22 = @(x1,x2) M_samp.F22(x1,x2);
end

for i=1:20
    [A, b, res] = AssembleLinearSystemConserve(t1, t2, dx1dt1, dx1dt2, dx2dt1, dx2dt2,...
                                                        M11, M12, M22,x1_in_t,x2_in_t);
    
    t_new = A \ b;
    t1_new = t_new(1:N); t2_new = t_new(N+1:end);
    t1_new = reshape(t1_new, Nx2, Nx1);
    t2_new = reshape(t2_new, Nx2, Nx1);

    t1 = t1 + 1*(t1_new-t1);
    t2 = t2 + 1*(t2_new-t2);
    figure(2)
    plot(t1, t2, 'k'); hold on; plot(t1', t2', 'k'); hold off
    pause(0.1)
end

figure
plot(x1_in_t(t2,t1), x2_in_t(t2,t1), 'k'); hold on; plot(x1_in_t(t2,t1)', x2_in_t(t2,t1)', 'k'); hold off
%plot(x1_in_t(t1,t2), x2_in_t(t1,t2), 'k'); hold on; plot(x1_in_t(t1,t2)', x2_in_t(t1,t2)', 'k'); hold off
%%
function [A, b, res] = AssembleLinearSystemConserve(t1, t2, dx1dt1Int, dx1dt2Int, dx2dt1Int, dx2dt2Int,...
                                                    M11Int, M12Int, M22Int,x1_in_t,x2_in_t)
[Nt2, Nt1] = size(t1);
N = Nt1*Nt2;
N_all = 2*N;

ds1 = 1 / Nt1;
ds2 = 1 / Nt2;
sigma1 = 1 / Nt1;
sigma2 = 1 / Nt2;

M11 = zeros(Nt2+2,Nt1+2); M12 = zeros(Nt2+2,Nt1+2); M22 = zeros(Nt2+2,Nt1+2);
dx1dt1 = zeros(Nt2+2,Nt1+2); dx1dt2 = zeros(Nt2+2,Nt1+2); dx2dt1 = zeros(Nt2+2,Nt1+2); dx2dt2 = zeros(Nt2+2,Nt1+2);

[dt1ds1, dt2ds2] = DCentral(t1, t2, ds1, ds2);
[dt2ds1, dt1ds2] = DCentral(t2, t1, ds1, ds2);

%M11(2:end-1,2:end-1) = M11Int(x1_in_t(t1',t2'), x2_in_t(t1',t2'));
%M12(2:end-1,2:end-1) = M12Int(x1_in_t(t1',t2'), x2_in_t(t1',t2'));
%M22(2:end-1,2:end-1) = M22Int(x1_in_t(t1',t2'), x2_in_t(t1',t2'));

%dx1dt1(2:end-1,2:end-1) = dx1dt1Int(t1',t2');
%dx1dt2(2:end-1,2:end-1) = dx1dt2Int(t1',t2');
%dx2dt1(2:end-1,2:end-1) = dx2dt1Int(t1',t2');
%dx2dt2(2:end-1,2:end-1) = dx2dt2Int(t1',t2');

%A11 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt1 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt1;
%A12 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt2 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt2;
%A22 = (dx1dt2.*M11+dx2dt2.*M12).*dx1dt2 + (dx1dt2.*M12+dx2dt2.*M22).*dx2dt2;
%A21 = (dx1dt2.*M11+dx2dt2.*M12).*dx1dt1 + (dx1dt2.*M12+dx2dt2.*M22).*dx2dt1;

%M11(2:end-1,2:end-1) = M11Int(x1_in_t(t1',t2'), x2_in_t(t1',t2'));
%M12(2:end-1,2:end-1) = M12Int(x1_in_t(t1',t2'), x2_in_t(t1',t2'));
%M22(2:end-1,2:end-1) = M22Int(x1_in_t(t1',t2'), x2_in_t(t1',t2'));

M11(2:end-1,2:end-1) = M11Int(t1,t2);
M12(2:end-1,2:end-1) = M12Int(t1,t2);
M22(2:end-1,2:end-1) = M22Int(t1,t2);
A11 = M11;
A12 = M12;
A22 = M22;

M11_f1 = (A11(:,1:end-1) + A11(:,2:end)) / 2;
M11_f2 = (A11(1:end-1,:) + A11(2:end,:)) / 2;
M12_f1 = (A12(:,1:end-1) + A12(:,2:end)) / 2;
M12_f2 = (A12(1:end-1,:) + A12(2:end,:)) / 2;
M22_f1 = (A22(:,1:end-1) + A22(:,2:end)) / 2;
M22_f2 = (A22(1:end-1,:) + A22(2:end,:)) / 2;

[dM11dt1, dM11dt2] = metric_grad(A11(2:end-1,2:end-1), t1, t2);
[dM12dt1, dM12dt2] = metric_grad(A12(2:end-1,2:end-1), t1, t2);
[dM22dt1, dM22dt2] = metric_grad(A22(2:end-1,2:end-1), t1, t2);



e = ones(N,1);

D11_coef = -(M11_f1(2:end-1,1:end-1)+M11_f1(2:end-1,2:end))*(sigma1/ds1)^2 ...
           -(M11_f2(1:end-1,2:end-1)+M11_f2(2:end,2:end-1))*(sigma2/ds2)^2;
D11 = spdiags(e, 0, N, N) .* D11_coef(:);

D12_coef = -(M12_f1(2:end-1,1:end-1)+M12_f1(2:end-1,2:end))*(sigma1/ds1)^2 ...
           -(M12_f2(1:end-1,2:end-1)+M12_f2(2:end,2:end-1))*(sigma2/ds2)^2;
D12 = spdiags(e, 0, N, N) .* D12_coef(:);

U1_11_coef = M11_f2(2:end,2:end-1) * (sigma2/ds2)^2;
U1_11 = spdiags(e, 1, N, N) .* U1_11_coef(:);

L1_11_coef = M11_f2(1:end-1,2:end-1) * (sigma2/ds2)^2;
L1_11 = spdiags(e, -1, N, N) .* L1_11_coef(:);

Ux2_11_coef = M11_f1(2:end-1,2:end) * (sigma1/ds1)^2;
Ux2_11 = spdiags(e, Nt2, N, N) .* Ux2_11_coef(:);

Lx2_11_coef = M11_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_11 = spdiags(e, -Nt2, N, N) .* Lx2_11_coef(:);

U1_12_coef = M12_f2(2:end,2:end-1) * (sigma2/ds2)^2;
U1_12 = spdiags(e, 1, N, N) .* U1_12_coef(:);

L1_12_coef = M12_f2(1:end-1,2:end-1)  * (sigma2/ds2)^2;
L1_12 = spdiags(e, -1, N, N) .* L1_12_coef(:);

Ux2_12_coef = M12_f1(2:end-1,2:end) * (sigma1/ds1)^2;
Ux2_12 = spdiags(e, Nt2, N, N) .* Ux2_12_coef(:);

Lx2_12_coef = M12_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_12 = spdiags(e, -Nt2, N, N) .* Lx2_12_coef(:);

A = sparse(N_all,N_all);
A(1:N,1:N) = D11 + U1_11 + L1_11 + Ux2_11 + Lx2_11;
A(1:N,N+1:end) = D12 + U1_12 + L1_12 + Ux2_12 + Lx2_12;



D22_coef = -(M22_f1(2:end-1,1:end-1)+M22_f1(2:end-1,2:end)) * (sigma1/ds1)^2 ...
           -(M22_f2(1:end-1,2:end-1)+M22_f2(2:end,2:end-1)) * (sigma2/ds2)^2;
D22 = spdiags(e, 0, N, N) .* D22_coef(:);

D21_coef = -(M12_f1(2:end-1,1:end-1)+M12_f1(2:end-1,2:end)) * (sigma1/ds1)^2 ...
           -(M12_f2(1:end-1,2:end-1)+M12_f2(2:end,2:end-1)) * (sigma2/ds2)^2;
D21 = spdiags(e, 0, N, N) .* D21_coef(:);

U1_22_coef = M22_f2(2:end,2:end-1) * (sigma2/ds2)^2;
U1_22 = spdiags(e, 1, N, N) .* U1_22_coef(:);

L1_22_coef = M22_f2(1:end-1,2:end-1) * (sigma2/ds2)^2;
L1_22 = spdiags(e, -1, N, N) .* L1_22_coef(:);

Ux2_22_coef = M22_f1(2:end-1,2:end) * (sigma1/ds1)^2;
Ux2_22 = spdiags(e, Nt2, N, N) .* Ux2_22_coef(:);

Lx2_22_coef = M22_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_22 = spdiags(e, -Nt2, N, N) .* Lx2_22_coef(:);

U1_21_coef = M12_f2(2:end,2:end-1) * (sigma2/ds2)^2;
U1_21 = spdiags(e, 1, N, N) .* U1_21_coef(:);

L1_21_coef = M12_f2(1:end-1,2:end-1)  * (sigma2/ds2)^2;
L1_21 = spdiags(e, -1, N, N) .* L1_21_coef(:);

Ux2_21_coef = M12_f1(2:end-1,2:end) * (sigma1/ds1)^2;
Ux2_21 = spdiags(e, Nt2, N, N) .* Ux2_21_coef(:);

Lx2_21_coef = M12_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_21 = spdiags(e, -Nt2, N, N) .* Lx2_21_coef(:);

A(N+1:end,N+1:end) = D22 + U1_22 + L1_22 + Ux2_22 + Lx2_22;
A(N+1:end,1:N) = D21 + U1_21 + L1_21 + Ux2_21 + Lx2_21;
A = 2*A;

% Apply boundary conditions
id = GetIndex(Nt1, Nt2);
   
A(id.l,:) = 0; A(id.l+N,:) = 0;
A(id.r,:) = 0; A(id.r+N,:) = 0;
A(id.b,:) = 0; 
A(id.b+N,:) = 0;
A(id.t,:) = 0;
A(id.t+N,:) = 0;

A(id.l,id.l) = eye(length(id.l));
A(id.r,id.r) = eye(length(id.r));
A(id.b+N,id.b+N) = eye(length(id.b));
A(id.t+N,id.t+N) = eye(length(id.t));
for i=id.b
    A(i,i) = 1; A(i,i+1) = -1;
end
for i=id.t
    A(i,i) = 1; A(i,i-1) = -1;
end
for i=id.l+N
    A(i,i) = 1; A(i,i+Nt2) = -1;
end
for i=id.r+N
    A(i,i) = 1; A(i,i-Nt2) = -1;
end


for i = id.corner; A(i, :) = 0; A(i, i) = 1; end
for i = N+id.corner; A(i, :) = 0; A(i, i) = 1; end

% Assemble the interior vector b
b1 = zeros(Nt2, Nt1);
b2 = zeros(Nt2, Nt1);

b1 = b1 + dM22dt1.*(dt2ds1.^2*sigma1^2 + dt2ds2.^2*sigma2^2);
b1 = b1 + 2*dM12dt1.*(dt1ds1.*dt2ds1*sigma1^2 + dt1ds2.*dt2ds2*sigma2^2);
b1 = b1 + dM11dt1.*(dt1ds1.^2*sigma1^2 + dt1ds2.^2*sigma2^2);

b2 = b2 + dM22dt2.*(dt2ds1.^2*sigma1^2 + dt2ds2.^2*sigma2^2);
b2 = b2 + 2*dM12dt2.*(dt1ds1.*dt2ds1*sigma1^2 + dt1ds2.*dt2ds2*sigma2^2);
b2 = b2 + dM11dt2.*(dt1ds1.^2*sigma1^2 + dt1ds2.^2*sigma2^2);


b1(:,1) = 0;
b1(:,end) = 1;
b1(1,:) = 0;
b1(end,:) = 0;
b2(:,1) = 0;
b2(:,end) = 0;
b2(1,:) = 0;
b2(end,:) = 1;

b1(1,1) = 0; b1(1,end) = 1; b1(end,1) = 0; b1(end,end) = 1;
b2(1,1) = 0; b2(1,end) = 0; b2(end,1) = 1; b2(end,end) = 1;


b = [b1(:); b2(:)];
res = A*[t1(:); t2(:)] - b;

end

function [dMdx1, dMdx2] = metric_grad(M, x1, x2)
%METRIC_GRADIENT  Compute ∂M/∂x₁ and ∂M/∂x₂ on a non‐uniform, skewed grid
%
%   [dMdx1, dMdx2] = metric_gradient(M, x1, x2)
%
%   Inputs:
%     M  – a 2D array of size (nrows × ncols), representing one component of the metric tensor
%     x1 – an (nrows × ncols) array of the physical x₁‐coordinates at each node
%     x2 – an (nrows × ncols) array of the physical x₂‐coordinates at each node
%
%   Outputs:
%     dMdx1 – an (nrows × ncols) array containing ∂M/∂x₁ at interior nodes (zeros on the 4‐point boundary)
%     dMdx2 – an (nrows × ncols) array containing ∂M/∂x₂ at interior nodes (zeros on the 4‐point boundary)
%
%   This function uses a 5‐point, chain‐rule stencil.  At each interior index (i,j):
%     – Δξx₁ = x₁(i+1,j) − x₁(i−1,j)
%     – Δξx₂ = x₂(i+1,j) − x₂(i−1,j)
%     – Δηx₁ = x₁(i,j+1) − x₁(i,j−1)
%     – Δηx₂ = x₂(i,j+1) − x₂(i,j−1)
%     – ΔξM  = M(i+1,j)   − M(i−1,j)
%     – ΔηM  = M(i,j+1)   − M(i,j−1)
%
%   Then one solves the 2×2 linear system:
%       [ Δξx₁   Δξx₂ ] [∂M/∂x₁]   = [ ΔξM ]
%       [ Δηx₁   Δηx₂ ] [∂M/∂x₂]     [ ΔηM ]
%
%   This is vectorized so that all interior points are handled in one shot.

[nrows, ncols] = size(M);

% Pre‐allocate output arrays (zeros on boundary)
dMdx1 = zeros(nrows, ncols);
dMdx2 = zeros(nrows, ncols);

% Define index ranges for “interior” (skip the 1‐pixel boundary)
%   i = 2:(nrows−1),  j = 2:(ncols−1)
%
% To vectorize, we form sub‐blocks of size (nrows−2)×(ncols−2):
%   x1(i+1,j) corresponds to x1(3:end,   2:end−1)
%   x1(i−1,j) corresponds to x1(1:end−2, 2:end−1)
%   x1(i,j+1) corresponds to x1(2:end−1, 3:end)
%   x1(i,j−1) corresponds to x1(2:end−1, 1:end−2)
%   etc.

% 1) Compute Δξx₁, Δξx₂, Δηx₁, Δηx₂ over all interior points at once
deta_x1  = x1(3:end,   2:end-1) - x1(1:end-2, 2:end-1);   % (i+1,j) − (i−1,j)
deta_x2  = x2(3:end,   2:end-1) - x2(1:end-2, 2:end-1);

dxi_x1 = x1(2:end-1, 3:end  ) - x1(2:end-1, 1:end-2);   % (i,j+1) − (i,j−1)
dxi_x2 = x2(2:end-1, 3:end  ) - x2(2:end-1, 1:end-2);

% 2) Compute ΔξM and ΔηM over the same interior stencil
deta_M  = M(3:end,   2:end-1) - M(1:end-2, 2:end-1);   % M(i+1,j) − M(i−1,j)
dxi_M = M(2:end-1, 3:end  ) - M(2:end-1, 1:end-2);   % M(i,j+1) − M(i,j−1)

% 3) Build the determinant of the 2×2 Jacobian for each interior point
%       detJ = Δξx₁ * Δηx₂ − Δξx₂ * Δηx₁
detJ = dxi_x1 .* deta_x2 - dxi_x2 .* deta_x1;

% 4) Compute each component of the inverse of J = [Δξx₁ Δξx₂; Δηx₁ Δηx₂]
%    inv(J) = (1/detJ) * [ Δηx₂  −Δξx₂;  −Δηx₁  Δξx₁ ]
inv11 =  deta_x2 ./ detJ;    % coefficient mapping ΔξM → ∂M/∂x₁
inv12 = -dxi_x2 ./ detJ;     % coefficient mapping ΔηM → ∂M/∂x₁
inv21 = -deta_x1 ./ detJ;    % coefficient mapping ΔξM → ∂M/∂x₂
inv22 =  dxi_x1 ./ detJ;     % coefficient mapping ΔηM → ∂M/∂x₂

% 5) Compute ∂M/∂x₁ and ∂M/∂x₂ on the interior block
dMdx1_int = inv11 .* dxi_M + inv12 .* deta_M;    % (nrows−2)×(ncols−2)
dMdx2_int = inv21 .* dxi_M + inv22 .* deta_M;

% 6) Scatter back into the full‐size arrays (leaving zeros on the 1‐pixel boundary)
dMdx1(2:end-1, 2:end-1) = dMdx1_int;
dMdx2(2:end-1, 2:end-1) = dMdx2_int;

end
