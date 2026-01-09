clear
problem = 2;
orthongonal_project = 1;

Nx1 = 101;
Nx2 = 101;

[x1,x2,Mfun] = problems.Initialization(problem,Nx1,Nx2);


boundary_points = extract_boundary_points(x1,x2);

figure
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');

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
detJ = Mp.M11.*Mp.M22 - Mp.M12.^2;
Maux.M11 = Mp.M22 ./ detJ;
Maux.M22 = Mp.M11 ./ detJ;
Maux.M12 = -Mp.M12 ./ detJ;
Mp_inv = Maux;
%{
for i=1:20
    Mp_inv.M11(:,1) = Mp_inv.M11(:,2); Mp_inv.M11(:,end) = Mp_inv.M11(:,end-1);
    %Mp_inv.M22(1,:) = Mp_inv.M22(2,:); Mp_inv.M22(end,:) = Mp_inv.M22(end-1,:);

    Mp_inv.M11(:,2:end-1) = (Mp_inv.M11(:,1:end-2) + Mp_inv.M11(:,3:end)) / 2;
    %Mp_inv.M22(:,2:end-1) = (Mp_inv.M22(:,1:end-2) + Mp_inv.M22(:,3:end)) / 2;
    Mp_inv.M11(2:end-1,:) = (Mp_inv.M11(1:end-2,:) + Mp_inv.M11(3:end,:)) / 2;
    %Mp_inv.M22(2:end-1,:) = (Mp_inv.M22(1:end-2,:) + Mp_inv.M22(3:end,:)) / 2;
end

Mp_inv.M11 = Mp_inv.M11 + 0.01*max(Mp_inv.M11,[],'all');
Mp_inv.M22 = Mp_inv.M22 + 0.01*max(Mp_inv.M22,[],'all');
%}
Mp11_inv = griddedInterpolant(t1, t2, Mp_inv.M11', 'cubic');
Mp22_inv = griddedInterpolant(t1, t2, Mp_inv.M22', 'cubic');
Mp12_inv = griddedInterpolant(t1, t2, Mp_inv.M12', 'cubic');
%%

Nt1 = 50; Nt2 = 50;
%Nt1 = 10; Nt2 = 10;
N = Nt1*Nt2;
S1 = linspace(0,1,Nt1);
S2 = linspace(0,1,Nt2);
[S1, S2] = meshgrid(S1, S2);
t1 = linspace(0,1,Nt1);
t2 = linspace(0,1,Nt2);
[t1, t2] = meshgrid(t1, t2);



%%
[s1, s2] = AssembleLinearSystemDual(t1, t2, Mp11_inv, Mp12_inv, Mp22_inv);
%s1(:,2:end-1) = -2*s1(:,2:end-1);
%s2(2:end-1,:) = -2*s2(2:end-1,:);

s1_samp = linspace(0,1,50);
s2_samp = linspace(0,1,50);
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

Mpfun.M11 = Mp11;
Mpfun.M12 = Mp12;
Mpfun.M22 = Mp22;

Mpfun_inv.M11 = Mp11_inv;
Mpfun_inv.M12 = Mp12_inv;
Mpfun_inv.M22 = Mp22_inv;

[Lx1, Lx2] = Cost(t1_samp, t2_samp, Mpfun);
[Lx1_dual, Lx2_dual] = Cost_dual(t1_samp, t2_samp, Mpfun_inv);
figure
pcolor(x1_in_t(t1_samp,t2_samp),x2_in_t(t1_samp,t2_samp),Lx1_dual+Lx2_dual);
figure
pcolor(x1_in_t(t1_samp,t2_samp),x2_in_t(t1_samp,t2_samp),Lx1+Lx2);


%%
p = 1;
Nt1 = 50;
Nt2 = 50;
N = Nt1*Nt2;
t1_primal = linspace(0,1,Nt1);
t2_primal = linspace(0,1,Nt2);
[t1_primal, t2_primal] = meshgrid(t1_primal, t2_primal);
exact = 1;

for i=1:500
    [A, b, res] = AssembleLinearSystemConserve(t1_primal, t2_primal, dx1dt1, dx1dt2, dx2dt1, dx2dt2,...
                                                        Mp11, Mp12, Mp22,x1_in_t,x2_in_t,exact,p);
    
    t_new = A \ b;
    t1_new = t_new(1:N); t2_new = t_new(N+1:end);
    t1_new = reshape(t1_new, Nt2, Nt1);
    t2_new = reshape(t2_new, Nt2, Nt1);

    t1_primal = t1_primal + 0.05*(t1_new-t1_primal);
    t2_primal = t2_primal + 0.05*(t2_new-t2_primal);
    if mod(i,10) == 0
        figure(10)
        plot(t1_primal, t2_primal, 'k'); hold on; plot(t1_primal', t2_primal', 'k'); hold off
        figure(11)
        plot(x1_in_t(t1_primal,t2_primal), x2_in_t(t1_primal,t2_primal), 'k'); hold on; plot(x1_in_t(t1_primal,t2_primal)', x2_in_t(t1_primal,t2_primal)', 'k'); hold off
        pause(0.1)
    
        [Lx1, Lx2] = Cost(t1_primal, t2_primal, Mpfun);
        Lx = Lx1 + Lx2;
        disp(['misfit mean: ', num2str(mean(Lx(:)))]);
        disp(['misfit range: ', num2str(range(Lx(:)))]);
    end
end

figure
S = linspace(0,1,Nt1);
scatter(S, t1_primal(1,:)); hold on
plot(S, S ./ (16 - 15*S));

a = 1600*log(16)^2 / 9;
plot(S, 1/15* (-1 + exp(3*sqrt(a)*S/40)) );

%% Test with increasing p-norm
Nt1 = 50;
Nt2 = 50;
N = Nt1*Nt2;
t1_primal = linspace(0,1,Nt1);
t2_primal = linspace(0,1,Nt2);
[t1_primal, t2_primal] = meshgrid(t1_primal, t2_primal);
exact = 1;

p_list = [1, 2, 4, 8, 16];
X_S = zeros(length(p_list), Nt1);
for j=1:length(p_list)
p = p_list(j);
disp(['Now minimizing L-', num2str(p), ' norm of misfit']);
i_max = 50;
if p==1; i_max=200; end
for i=1:i_max
    [A, b, res] = AssembleLinearSystemConserve(t1_primal, t2_primal, dx1dt1, dx1dt2, dx2dt1, dx2dt2,...
                                                        Mp11, Mp12, Mp22,x1_in_t,x2_in_t,exact,p);
    
    t_new = A \ b;
    t1_new = t_new(1:N); t2_new = t_new(N+1:end);
    t1_new = reshape(t1_new, Nt2, Nt1);
    t2_new = reshape(t2_new, Nt2, Nt1);

    t1_primal = t1_primal + 0.05*(t1_new-t1_primal);
    t2_primal = t2_primal + 0.05*(t2_new-t2_primal);

    if mod(i,10) == 0
    [Lx1, Lx2] = Cost(t1_primal, t2_primal, Mpfun);
    Lx = Lx1 + Lx2;

    figure(10)
    plot(t1_primal, t2_primal, 'k'); hold on; plot(t1_primal', t2_primal', 'k'); hold off
    title(['Now minimizing L-', num2str(p), ' norm of misfit'])

    figure(11)
    pcolor(t1_primal, t2_primal, Lx1+Lx2); clim([0.5 1.5]);
    title('Local misfit L1+L2'); colorbar

    figure(12)
    contourf(t1_primal, t2_primal, Lx1+Lx2,0.5:0.1:1.5); clim([0.5 1.5]);
    title('Local misfit L1+L2'); colorbar
    pause(0.05)
    end

end
figure(12)
    contourf(t1_primal, t2_primal, Lx1+Lx2,0:0.1:2); clim([0.5 1.5]);
    title('Local misfit L1+L2'); colorbar

disp(['L-', num2str(p), ' optimization result:']);
disp(['misfit mean: ', num2str(mean(Lx(:)))]);
disp(['misfit range: ', num2str(range(Lx(:)))]);
disp(' ');
X_S(j,:) = t1_primal(1,:);
end

figure
S = linspace(0,1,Nt1);
scatter(S, X_S); hold on
plot(S, S ./ (16 - 15*S));

a = 1600*log(16)^2 / 9;
plot(S, 1/15* (-1 + exp(3*sqrt(a)*S/40)) );
legend("p1", "p2", "p4", "p8", "p16","p1 Exact", "p-inf Exact", Location='best');
%%
Nt1 = 50;
Nt2 = 50;
N = Nt1*Nt2;
x1_equi = linspace(0,1,Nt1);
x2_equi = linspace(0,1,Nt2);
[x1_equi, x2_equi] = meshgrid(x1_equi, x2_equi);

for i=1:500
    [Lx1, Lx2] = Cost(x1_equi, x2_equi, Mpfun);
    %figure(10)
    %pcolor(x1_equi, x2_equi, Lx1);
    %pause(0.1)

    [dL1dx1, dL1dx2] = gradient(Lx1);
    [dL2dx1, dL2dx2] = gradient(Lx2);
    
    dL1dx1(:,1) = 0; dL1dx1(:,end) = 0; dL1dx2(1,:) = 0; dL1dx2(end,:) = 0;
    dL2dx1(:,1) = 0; dL2dx1(:,end) = 0; dL2dx2(1,:) = 0; dL2dx2(end,:) = 0;

    x1_equi = x1_equi + 0.01*(dL1dx1);
    x2_equi = x2_equi + 0.01*(dL2dx2);
end
figure(10)
pcolor(x1_equi, x2_equi, Lx1 + Lx2);


%%
function [A, b, res] = AssembleLinearSystemConserve(t1, t2, dx1dt1Int, dx1dt2Int, dx2dt1Int, dx2dt2Int,...
                                                    M11Int, M12Int, M22Int,x1_in_t,x2_in_t,exact,p)
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

if exact
    M.M11 = M11Int(t1,t2);
    M.M12 = M12Int(t1,t2);
    M.M22 = M22Int(t1,t2);

    m1 = zeros(Nt2+2,Nt1+2); m2 = zeros(Nt2+2,Nt1+2); 
    m1(2:end-1,2:end-1) = sigma1^2*(M.M11.*dt1ds1.^2 + 2*M.M12.*dt1ds1.*dt2ds1 + M.M22.*dt2ds1.^2) / 1e3;
    m2(2:end-1,2:end-1) = sigma2^2*(M.M11.*dt1ds2.^2 + 2*M.M12.*dt1ds2.*dt2ds2 + M.M22.*dt2ds2.^2) / 1e3;
    %m1(2:end-1,2:end-1) = sigma1^2*(M.M11.*dt1ds1.^2 + 2*M.M12.*dt1ds1.*dt2ds1 + M.M22.*dt2ds1.^2) / 1e3 - 1;
    %m2(2:end-1,2:end-1) = sigma2^2*(M.M11.*dt1ds2.^2 + 2*M.M12.*dt1ds2.*dt2ds2 + M.M22.*dt2ds2.^2) / 1e3 - 1;

    m1(2:end-1,2:end-1) = m1(2:end-1,2:end-1).^(p-1);
    m2(2:end-1,2:end-1) = m2(2:end-1,2:end-1).^(p-1);

    M11_f1 = M11_f1.*(m1(:,1:end-1)/2 + m1(:,2:end)/2);
    M11_f2 = M11_f2.*(m2(1:end-1,:)/2 + m2(2:end,:)/2);
    M12_f1 = M12_f1.*(m1(:,1:end-1)/2 + m1(:,2:end)/2);
    M12_f2 = M12_f2.*(m2(1:end-1,:)/2 + m2(2:end,:)/2);
    M22_f1 = M22_f1.*(m1(:,1:end-1)/2 + m1(:,2:end)/2);
    M22_f2 = M22_f2.*(m2(1:end-1,:)/2 + m2(2:end,:)/2);
end

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
if exact
    m1 = m1(2:end-1,2:end-1);
    m2 = m2(2:end-1,2:end-1);
else
    m1 = ones(Nt2,Nt1);
    m2 = ones(Nt2,Nt1);
end

b1 = zeros(Nt2, Nt1);
b2 = zeros(Nt2, Nt1);

b1 = b1 + dM22dt1.*(m1.*dt2ds1.^2*sigma1^2 + m2.*dt2ds2.^2*sigma2^2);
b1 = b1 + 2*dM12dt1.*(m1.*dt1ds1.*dt2ds1*sigma1^2 + m2.*dt1ds2.*dt2ds2*sigma2^2);
b1 = b1 + dM11dt1.*(m1.*dt1ds1.^2*sigma1^2 + m2.*dt1ds2.^2*sigma2^2);

b2 = b2 + dM22dt2.*(m1.*dt2ds1.^2*sigma1^2 + m2.*dt2ds2.^2*sigma2^2);
b2 = b2 + 2*dM12dt2.*(m1.*dt1ds1.*dt2ds1*sigma1^2 + m2.*dt1ds2.*dt2ds2*sigma2^2);
b2 = b2 + dM11dt2.*(m1.*dt1ds1.^2*sigma1^2 + m2.*dt1ds2.^2*sigma2^2);


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


