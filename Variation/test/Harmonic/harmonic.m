clear
problem = 5;
orthongonal_project = 1;

Nx1 = 151;
Nx2 = 51;

[x1,x2,Mfun] = problems.Initialization(problem,Nx1,Nx2);


boundary_points = extract_boundary_points(x1,x2);

figure
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
%[Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1, x2);

[x1, x2, info] = solve_harmonic(x1, x2, Omega=0.5, ShowPlot=true, ...
                                BoundaryPoints=boundary_points, PauseTime=0.1);
%[Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1, x2);

%%

[Nx2, Nx1] = size(x1);
% Interpolant of X(T)
[T1,T2] = ndgrid(linspace(0,1,Nx1),linspace(0,1,Nx2));
x1_in_t = griddedInterpolant(T1, T2, x1', 'cubic');
x2_in_t = griddedInterpolant(T1, T2, x2', 'cubic');


% Interpolant of J = dX(t)/dT
[dx1dt1_samp, dx2dt2_samp] = DCentral(x1, x2, 1/Nx1, 1/Nx2);
[dx2dt1_samp, dx1dt2_samp] = DCentral(x2, x1, 1/Nx1, 1/Nx2);

dx1dt1Fun = griddedInterpolant(T1, T2, dx1dt1_samp', 'cubic');
dx1dt2Fun = griddedInterpolant(T1, T2, dx1dt2_samp', 'cubic');
dx2dt1Fun = griddedInterpolant(T1, T2, dx2dt1_samp', 'cubic');
dx2dt2Fun = griddedInterpolant(T1, T2, dx2dt2_samp', 'cubic');

J_samp = dx1dt1_samp.*dx2dt2_samp - dx2dt1_samp.*dx1dt2_samp;


%%

t1 = linspace(0,1,101);
t2 = linspace(0,1,101);
[t1, t2] = ndgrid(t1, t2);

x1 = x1_in_t(t1,t2)';
x2 = x2_in_t(t1,t2)';

M = Mfun(x1, x2);

dx1dt1 = dx1dt1Fun(t1, t2)';
dx1dt2 = dx1dt2Fun(t1, t2)';
dx2dt1 = dx2dt1Fun(t1, t2)';
dx2dt2 = dx2dt2Fun(t1, t2)';

Mp.M11 = (dx1dt1.*M.M11+dx2dt1.*M.M12).*dx1dt1 + (dx1dt1.*M.M12+dx2dt1.*M.M22).*dx2dt1;
Mp.M22 = (dx1dt2.*M.M11+dx2dt2.*M.M12).*dx1dt2 + (dx1dt2.*M.M12+dx2dt2.*M.M22).*dx2dt2;
Mp.M12 = (dx1dt1.*M.M11+dx2dt1.*M.M12).*dx1dt2 + (dx1dt1.*M.M12+dx2dt1.*M.M22).*dx2dt2;


Mp11 = griddedInterpolant(t1, t2, Mp.M11', 'cubic');
Mp22 = griddedInterpolant(t1, t2, Mp.M22', 'cubic');
Mp12 = griddedInterpolant(t1, t2, Mp.M12', 'cubic');

%analysis.plot_metric(x1, x2, M);
%analysis.plot_metric(t1', t2', Mp);
%% Solve inverse instead
detJ = Mp.M11.*M.M22 - Mp.M12.^2;
Maux.M11 = Mp.M22 ./ detJ;
Maux.M22 = Mp.M11 ./ detJ;
Maux.M12 = -Mp.M12 ./ detJ;
Mp = Maux;

for i=1:10
    Mp.M11(:,1) = Mp.M11(:,2); Mp.M11(:,end) = Mp.M11(:,end-1);
    Mp.M22(1,:) = Mp.M22(2,:); Mp.M22(end,:) = Mp.M22(end-1,:);

    Mp.M11(:,2:end-1) = (Mp.M11(:,1:end-2) + Mp.M11(:,3:end)) / 2;
    Mp.M22(:,2:end-1) = (Mp.M22(:,1:end-2) + Mp.M22(:,3:end)) / 2;
    Mp.M11(2:end-1,:) = (Mp.M11(1:end-2,:) + Mp.M11(3:end,:)) / 2;
    Mp.M22(2:end-1,:) = (Mp.M22(1:end-2,:) + Mp.M22(3:end,:)) / 2;
end


Mp11 = griddedInterpolant(t1, t2, Mp.M11', 'cubic');
Mp22 = griddedInterpolant(t1, t2, Mp.M22', 'cubic');
Mp12 = griddedInterpolant(t1, t2, Mp.M12', 'cubic');
%%

Nt1 = 401; Nt2 = 201;
%Nt1 = 10; Nt2 = 10;
N = Nt1*Nt2;
S1 = linspace(0,1,Nt1);
S2 = linspace(0,1,Nt2);
[S1, S2] = meshgrid(S1, S2);
t1 = linspace(0,1,Nt1);
t2 = linspace(0,1,Nt2);
[t1, t2] = meshgrid(t1, t2);



%%
[s1, s2] = AssembleLinearSystemDual(t1, t2, Mp11, Mp12, Mp22);
%s1(:,2:end-1) = -2*s1(:,2:end-1);
%s2(2:end-1,:) = -2*s2(2:end-1,:);

s1_samp = linspace(0,1,100);
s2_samp = linspace(0,1,100);
[s1_samp,s2_samp] = meshgrid(s1_samp,s2_samp);
t1_samp = griddata(s1,s2,t1,s1_samp,s2_samp);
t2_samp = griddata(s1,s2,t2,s1_samp,s2_samp);

figure
plot(s1, s2, 'k'); hold on; plot(s1', s2', 'k'); hold off
figure
plot(t1_samp, t2_samp, 'k'); hold on; plot(t1_samp', t2_samp', 'k'); hold off
figure
plot(x1_in_t(t1_samp,t2_samp), x2_in_t(t1_samp,t2_samp), 'k'); hold on; plot(x1_in_t(t1_samp,t2_samp)', x2_in_t(t1_samp,t2_samp)', 'k'); hold off
%%


for i=1:200
    %[A, b, res] = AssembleLinearSystem(t1, t2, Mp11, Mp12, Mp22);
    [A, b, res] = AssembleLinearSystemConserve(t1, t2, dx1dt1, dx1dt2, dx2dt1, dx2dt2,...
                                                        Mp11, Mp12, Mp22,x1_in_t,x2_in_t);
    
    t_new = A \ b;
    t1_new = t_new(1:N); t2_new = t_new(N+1:end);
    t1_new = reshape(t1_new, Nt2, Nt1);
    t2_new = reshape(t2_new, Nt2, Nt1);

    t1 = t1 + 1*(t1_new-t1);
    t2 = t2 + 1*(t2_new-t2);
    figure(10)
    plot(t1, t2, 'k'); hold on; plot(t1', t2', 'k'); hold off
    figure(11)
    plot(x1_in_t(t1,t2), x2_in_t(t1,t2), 'k'); hold on; plot(x1_in_t(t1,t2)', x2_in_t(t1,t2)', 'k'); hold off
    pause(0.1)
end



%%
function [A, b, res] = AssembleLinearSystemConserve(t1, t2, dx1dt1Int, dx1dt2Int, dx2dt1Int, dx2dt2Int,...
                                                    M11Int, M12Int, M22Int,x1_in_t,x2_in_t)
[Nt2, Nt1] = size(t1);
N = Nt1*Nt2;
N_all = 2*N;

ds1 = 1 / Nt1;
ds2 = 1 / Nt2;
sigma1 = 1;
sigma2 = 1;

M11 = zeros(Nt2+2,Nt1+2); M12 = zeros(Nt2+2,Nt1+2); M22 = zeros(Nt2+2,Nt1+2);
dx1dt1 = zeros(Nt2+2,Nt1+2); dx1dt2 = zeros(Nt2+2,Nt1+2); dx2dt1 = zeros(Nt2+2,Nt1+2); dx2dt2 = zeros(Nt2+2,Nt1+2);

[dt1ds1, dt2ds2] = DCentral(t1, t2, ds1, ds2);
[dt2ds1, dt1ds2] = DCentral(t2, t1, ds1, ds2);


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

%M11_f1 = 2 ./ (1./A11(:,1:end-1) + 1./A11(:,2:end));
%M11_f2 = 2 ./ (1./A11(1:end-1,:) + 1./A11(2:end,:));
%M12_f1 = 2 ./ (1./A12(:,1:end-1) + 1./A12(:,2:end));
%M12_f2 = 2 ./ (1./A12(1:end-1,:) + 1./A12(2:end,:));
%M22_f1 = 2 ./ (1./A22(:,1:end-1) + 1./A22(:,2:end));
%M22_f2 = 2 ./ (1./A22(1:end-1,:) + 1./A22(2:end,:));


[dM11dt1, dM11dt2] = metric_grad(A11(2:end-1,2:end-1), t1, t2);
[dM12dt1, dM12dt2] = metric_grad(A12(2:end-1,2:end-1), t1, t2);
[dM22dt1, dM22dt2] = metric_grad(A22(2:end-1,2:end-1), t1, t2);

dM11dt1 = dM11dt1*0;
dM11dt2 = dM11dt2*0;
dM12dt1 = dM12dt1*0;
dM12dt2 = dM12dt2*0;
dM22dt1 = dM22dt1*0;
dM22dt2 = dM22dt2*0;



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


function [A, b, res] = AssembleLinearSystem(x1, x2, M11Int, M12Int, M22Int)

    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    M11 = M11Int(x1,x2);
    M12 = M12Int(x1,x2);
    M22 = M22Int(x1,x2);

    [dM11dx1, dM11dx2] = metric_grad(M11, x1, x2);
    [dM12dx1, dM12dx2] = metric_grad(M12, x1, x2);
    [dM22dx1, dM22dx2] = metric_grad(M22, x1, x2);

    dM11dx1 = dM11dx1*0;
    dM11dx2 = dM11dx2*0;
    dM12dx1 = dM12dx1*0;
    dM12dx2 = dM12dx2*0;
    dM22dx1 = dM22dx1*0;
    dM22dx2 = dM22dx2*0;


    [dx1ds1, dx2ds2] = DCentral(x1, x2, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, sigma1, sigma2);

    e = ones(N,1);
    L = spdiags([e, e, -4*e, e, e], [-Nx2, -1, 0, 1, Nx2], N, N);

    A = sparse(N_all, N_all);
    A(1:N,1:N) = A(1:N,1:N) - 2*L.*M11(:);
    A(N+1:end,N+1:end) = A(N+1:end,N+1:end) - 2*L.*M22(:);

    A(1:N,N+1:end) = A(1:N,N+1:end) - 2*L.*M12(:);
    A(N+1:end,1:N) = A(N+1:end,1:N) - 2*L.*M12(:);


    id = GetIndex(Nx1, Nx2);
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
        A(i,i) = 1; A(i,i+Nx2) = -1;
    end
    for i=id.r+N
        A(i,i) = 1; A(i,i-Nx2) = -1;
    end
    
    
    for i = id.corner; A(i, :) = 0; A(i, i) = 1; end
    for i = N+id.corner; A(i, :) = 0; A(i, i) = 1; end

    % Assemble the interior vector b
    b1 = zeros(Nx2, Nx1);
    b2 = zeros(Nx2, Nx1);
    

    b1 = b1 + dM11dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
    b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
    b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
    b1 = b1 + 2*dM12dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);


    b2 = b2 + dM22dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
    b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
    b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    b2 = b2 + 2*dM12dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);

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
    res = A*[x1(:); x2(:)] - b;
end


