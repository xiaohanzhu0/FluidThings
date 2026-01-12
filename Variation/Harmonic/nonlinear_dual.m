function [x1_samp, x2_samp]= nonlinear_dual()
params.problemId = 5;
params.Nx1 = 51;
params.Nx2 = 51;
params.Nt1 = 50;
params.Nt2 = 50;
params.sampleNs = 200;
params.omega1 = 0.5;
params.omega2 = 0.01;
params.max_iter = 10;
params.pauseTime = 0.1;
params.diagonal_reg = 0;
params.interpMethod = 'cubic';
params.doPlot = 1;
params.plotEvery = 1;

% requires +problems/Initialization.m
[x1,x2,M11_fun,M12_fun,M22_fun] = problems.Initialization(params.problemId,params.Nx1,params.Nx2);

% requires utils/extract_boundary_points.m
boundary_points = extract_boundary_points(x1,x2);

figure
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');

% requires solve_harmonic.m
[x1, x2, info] = solve_harmonic(x1, x2, Omega=params.omega1, ShowPlot=true, ...
                                BoundaryPoints=boundary_points, PauseTime=params.pauseTime);

%%
[params.Nx2, params.Nx1] = size(x1);
% Build interpolants for x(t) to evaluate x at non-grid (t1,t2)
[T1,T2] = ndgrid(linspace(0,1,params.Nx1),linspace(0,1,params.Nx2));
x1_in_t = griddedInterpolant(T1, T2, x1', params.interpMethod);
x2_in_t = griddedInterpolant(T1, T2, x2', params.interpMethod);

% Interpolant of J = dX(t)/dT
% requires utils/DCentral.m
[dx1dt1_samp, dx2dt2_samp] = DCentral(x1, x2, 1/(params.Nx1-1), 1/(params.Nx2-1));
[dx2dt1_samp, dx1dt2_samp] = DCentral(x2, x1, 1/(params.Nx1-1), 1/(params.Nx2-1));


%%
[T1,T2] = meshgrid(linspace(0,1,params.Nx1),linspace(0,1,params.Nx2));
t1 = linspace(0,1,50);
t2 = linspace(0,1,50);
[t1, t2] = meshgrid(t1, t2);
dx1dt1 = interp2(T1,T2,dx1dt1_samp,t1,t2, params.interpMethod);
dx1dt2 = interp2(T1,T2,dx1dt2_samp,t1,t2, params.interpMethod);
dx2dt1 = interp2(T1,T2,dx2dt1_samp,t1,t2, params.interpMethod);
dx2dt2 = interp2(T1,T2,dx2dt2_samp,t1,t2, params.interpMethod);

x1_temp = x1_in_t(t1,t2);
x2_temp = x2_in_t(t1,t2);
M11 = M11_fun(x1_temp, x2_temp);
M12 = M12_fun(x1_temp, x2_temp);
M22 = M22_fun(x1_temp, x2_temp);

% Compute pushforward metric Mp = J^T M J on the t-grid.
Mp11 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt1 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt1;
Mp22 = (dx1dt2.*M11+dx2dt2.*M12).*dx1dt2 + (dx1dt2.*M12+dx2dt2.*M22).*dx2dt2;
Mp12 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt2 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt2;

Mp11_fun = griddedInterpolant(t1', t2', Mp11', params.interpMethod);
Mp12_fun = griddedInterpolant(t1', t2', Mp12', params.interpMethod);
Mp22_fun = griddedInterpolant(t1', t2', Mp22', params.interpMethod);
%% Solve inverse instead
detJ = Mp11.*Mp22 - Mp12.^2;
Mp_inv11 = Mp22 ./ detJ;
Mp_inv22 = Mp11 ./ detJ;
Mp_inv12 = -Mp12 ./ detJ;

if params.diagonal_reg == 1
    Mp_inv11 = Mp_inv11 + 0.01*max(Mp_inv11,[],'all');
    Mp_inv22 = Mp_inv22 + 0.01*max(Mp_inv22,[],'all');
end

%%
s1 = linspace(0,1,params.Nt1);
s2 = linspace(0,1,params.Nt2);
[s1, s2] = meshgrid(s1, s2);
t1_grid = linspace(0,1,params.Nt1);
t2_grid = linspace(0,1,params.Nt2);
[t1_grid, t2_grid] = meshgrid(t1_grid, t2_grid);

%% Solve for inverse map s(t) using conservative FV discretization.
hist = struct('Lt',zeros(1,params.max_iter), ...
              'Lx',zeros(1,params.max_iter), ...
              'Theta1',zeros(1,params.max_iter), ...
              'Theta2',zeros(1,params.max_iter), ...
              'ThetaInf',zeros(1,params.max_iter));

warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
for i = 1:params.max_iter

[ds1dt1, ds2dt2] = DCentral(s1, s2, 1/(params.Nt1-1), 1/(params.Nt2-1));
[ds2dt1, ds1dt2] = DCentral(s2, s1, 1/(params.Nt1-1), 1/(params.Nt2-1));
J = abs(ds1dt1.*ds2dt2 - ds2dt1.*ds1dt2);

Mp_inv11_fun = griddedInterpolant(t1_grid', t2_grid', Mp_inv11'.*J', params.interpMethod);
Mp_inv22_fun = griddedInterpolant(t1_grid', t2_grid', Mp_inv22'.*J', params.interpMethod);
Mp_inv12_fun = griddedInterpolant(t1_grid', t2_grid', Mp_inv12'.*J', params.interpMethod);

% requires AssembleLinearSystemDual.m
[s1_new, s2_new] = AssembleLinearSystemDual(t1_grid, t2_grid, Mp_inv11_fun, Mp_inv12_fun, Mp_inv22_fun);

s1 = s1 + params.omega2*(s1_new-s1);
s2 = s2 + params.omega2*(s2_new-s2);

s1_samp = linspace(0,1,params.sampleNs);
s2_samp = linspace(0,1,params.sampleNs);
[s1_samp,s2_samp] = meshgrid(s1_samp,s2_samp);
t1_samp = griddata(s1,s2,t1_grid,s1_samp,s2_samp);
t2_samp = griddata(s1,s2,t2_grid,s1_samp,s2_samp);
x1_samp = x1_in_t(t1_samp,t2_samp);
x2_samp = x2_in_t(t1_samp,t2_samp);

if params.doPlot && mod(i, params.plotEvery)==0
    plotIterationMeshes(s1,s2,t1_samp,t2_samp,x1_samp,x2_samp);
    pause(params.pauseTime);
end

M11_samp = Mp11_fun(t1_samp, t2_samp);
M12_samp = Mp12_fun(t1_samp, t2_samp);
M22_samp = Mp22_fun(t1_samp, t2_samp);
% requires utils/Cost.m
[Lt1, Lt2] = Cost(t1_samp, t2_samp, M11_samp, M12_samp, M22_samp);
hist.Lt(i) = mean(Lt1+Lt2,'all');

M11_samp = M11_fun(x1_samp, x2_samp);
M12_samp = M12_fun(x1_samp, x2_samp);
M22_samp = M22_fun(x1_samp, x2_samp);
[Lx1, Lx2] = Cost(x1_samp, x2_samp, M11_samp, M12_samp, M22_samp);
hist.Lx(i) = mean(Lx1+Lx2,'all');

% requires +analysis/skewness.m
[Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1_samp, x2_samp);
hist.Theta1(i) = Theta_1;
hist.Theta2(i) = Theta_2;
hist.ThetaInf(i) = Theta_inf;
end
warning('on', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
end

%{
plot(Lt_list1 / Lt_list1(1)); hold on
plot(Lx_list1 / Lx_list1(1)); hold on
plot(Theta_1_list1 / Theta_1_list1(1)); hold on
plot(Theta_2_list1 / Theta_2_list1(1)); hold on
plot(Theta_inf_list1 / Theta_inf_list1(1)); hold on
legend('Lt', 'Lx', '\Theta_1', '\Theta_2', '\Theta_{\infty}')
%%
max_iter = 5000;
Lt_list2 = zeros(1,max_iter);
Lx_list2 = zeros(1,max_iter);
Theta_1_list2 = zeros(1,max_iter);
Theta_2_list2 = zeros(1,max_iter);
Theta_inf_list2 = zeros(1,max_iter);

t1_samp_new = t1_samp;
t2_samp_new = t2_samp;
for i=1:max_iter
    t1_samp_new(2:end-1,:) = t1_samp_new(1:end-2,:)/4 + t1_samp_new(2:end-1,:)/2 + t1_samp_new(3:end,:)/4;
    t1_samp_new(1,:) = t1_samp_new(1,:)/2 + t1_samp_new(2,:)/2;
    t1_samp_new(end,:) = t1_samp_new(end,:)/2 + t1_samp_new(end-1,:)/2;

    t2_samp_new(:,2:end-1) = t2_samp_new(:,1:end-2)/4 + t2_samp_new(:,2:end-1)/2 + t2_samp_new(:,3:end)/4;
    t2_samp_new(:,1) = t2_samp_new(:,1)/2 + t2_samp_new(:,2)/2;
    t2_samp_new(:,end) = t2_samp_new(:,end)/2 + t2_samp_new(:,end-1)/2;

    M11_samp = Mp11_fun(t1_samp_new, t2_samp_new);
    M12_samp = Mp11_fun(t1_samp_new, t2_samp_new);
    M22_samp = Mp11_fun(t1_samp_new, t2_samp_new);
    [Lt1, Lt2] = Cost(t1_samp_new, t2_samp_new, M11_samp, M12_samp, M22_samp);
    Lt_list2(i) = mean(Lt1+Lt2,'all');

    x1_samp = x1_in_t(t1_samp_new,t2_samp_new);
    x2_samp = x2_in_t(t1_samp_new,t2_samp_new);
    M11_samp = M11_fun(x1_samp, x2_samp);
    M12_samp = M12_fun(x1_samp, x2_samp);
    M22_samp = M22_fun(x1_samp, x2_samp);
    [Lx1, Lx2] = Cost(x1_samp, x2_samp, M11_samp, M12_samp, M22_samp);
    Lx_list2(i) = mean(Lx1+Lx2,'all');

    [Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1_samp, x2_samp);
    Theta_1_list2(i) = Theta_1;
    Theta_2_list2(i) = Theta_2;
    Theta_inf_list2(i) = Theta_inf;
end

plot(Lt_list2 / Lt_list2(1)); hold on
plot(Lx_list2 / Lx_list2(1), LineWidth=2); hold on
plot(Theta_1_list2 / Theta_1_list2(1), LineWidth=2); hold on
plot(Theta_2_list2 / Theta_2_list2(1), LineWidth=2); hold on
plot(Theta_inf_list2 / Theta_inf_list2(1), LineWidth=2); hold on
legend('Lt','Lx', '\Theta_1', '\Theta_2', '\Theta_{\infty}')
xlabel('iterations');
ylabel('percentage from initial');
%%
figure
plot([Lt_list1,Lt_list2] / Lt_list1(1)); hold on
plot([Lx_list1,Lx_list2] / Lx_list1(1)); hold on
plot([Theta_1_list1,Theta_1_list2] / Theta_1_list1(1)); hold on
plot([Theta_2_list1,Theta_2_list2] / Theta_2_list1(1)); hold on
plot([Theta_inf_list1,Theta_inf_list2] / Theta_inf_list1(1)); hold on
legend('Lt','Lx', '\Theta_1', '\Theta_2', '\Theta_{\infty}')
%% Post-processing to promote orthogonality
t1_samp_new = t1_samp;
t2_samp_new = t2_samp;

for i=1:5000
t1_samp_new(2:end-1,:) = t1_samp_new(1:end-2,:)/4 + t1_samp_new(2:end-1,:)/2 + t1_samp_new(3:end,:)/4;
t1_samp_new(1,:) = t1_samp(1,:)/2 + t1_samp(2,:)/2;
t1_samp_new(end,:) = t1_samp(end,:)/2 + t1_samp(end-1,:)/2;

t2_samp_new(:,2:end-1) = t2_samp_new(:,1:end-2)/4 + t2_samp_new(:,2:end-1)/2 + t2_samp_new(:,3:end)/4;
t2_samp_new(:,1) = t2_samp(:,1)/2 + t2_samp(:,2)/2;
t2_samp_new(:,end) = t2_samp(:,end)/2 + t2_samp(:,end-1)/2;
end

figure(11)
plot(t1_samp_new, t2_samp_new, 'r'); hold on; plot(t1_samp_new', t2_samp_new', 'r'); hold off

%%

M11_samp = Mp11_fun(t1_samp, t2_samp);
M12_samp = Mp11_fun(t1_samp, t2_samp);
M22_samp = Mp11_fun(t1_samp, t2_samp);

[Lt1, Lt2] = Cost(t1_samp, t2_samp, M11_samp, M12_samp, M22_samp);
Lt = Lt1 + Lt2;
pcolor(t1_samp,t2_samp,Lx);
disp(['misfit mean: ', num2str(mean(Lx(:)))]);
disp(['misfit range: ', num2str(range(Lx(:)))]);


%%
M11_samp = M11_fun(x1_samp, x2_samp);
M12_samp = M12_fun(x1_samp, x2_samp);
M22_samp = M22_fun(x1_samp, x2_samp);

[Lx1, Lx2] = Cost(x1_samp, x2_samp, M11_samp, M12_samp, M22_samp);
Lx = Lx1 + Lx2;
pcolor(x1_samp,x2_samp,Lx);
disp(['misfit mean: ', num2str(mean(Lx(:)))]);
disp(['misfit range: ', num2str(range(Lx(:)))]);
%%
[Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1_samp, x2_samp);



%%
% Compute optimized mesh size
K11 = Mp_inv11_fun(t1,t2); K12 = Mp_inv12_fun(t1,t2); K22 = Mp_inv22_fun(t1,t2);
p1 = (K11.*ds1dt1.^2+K22.*ds2dt1.^2+2*K12.*ds1dt1.*ds2dt1);
p2 = (K22.*ds2dt2.^2+K11.*ds1dt2.^2+2*K12.*ds1dt2.*ds2dt2);
sigma1 = sqrt(sum(p1,'all') / sum(p1.^2,'all'));
sigma2 = sqrt(sum(p2,'all') / sum(p2.^2,'all'));

downscale = sqrt(10000 / (sigma1*sigma2));
sigma1 = round(sigma1*downscale);
sigma2 = round(sigma2*downscale);

%}
function plotIterationMeshes(s1,s2,t1_samp,t2_samp,x1_samp,x2_samp)
    figure(8)
    plot(s1, s2, 'k'); hold on; plot(s1', s2', 'k'); hold off
    figure(9)
    plot(t1_samp, t2_samp, 'k'); hold on; plot(t1_samp', t2_samp', 'k'); hold off
    figure(10)
    plot(x1_samp, x2_samp, 'k'); hold on; plot(x1_samp', x2_samp', 'k'); hold off
    axis equal
end

function [x1_new, x2_new] = post_orthogonalize(x1, x2)
    x1_new = x1;
    x2_new = x2;
    x1_new(2:end-1,2:end-1) = x1(2:end-1,2:end-1)/4 + ...
             (x1(1:end-2,2:end-1)+x1(2:end-1,1:end-2)+x1(3:end,2:end-1)+x1(2:end-1,3:end))/8 + ...
             (x1(1:end-2,1:end-2)+x1(1:end-2,3:end)+x1(3:end,1:end-2)+x1(3:end,3:end))/16;
    x2_new(2:end-1,2:end-1) = x2(2:end-1,2:end-1)/4 + ...
             (x2(1:end-2,2:end-1)+x2(2:end-1,1:end-2)+x2(3:end,2:end-1)+x2(2:end-1,3:end))/8 + ...
             (x2(1:end-2,1:end-2)+x2(1:end-2,3:end)+x2(3:end,1:end-2)+x2(3:end,3:end))/16;
    x1_new(1,:) = x1(1,:)/2 + x1(2,:)/2;
    x1_new(end,:) = x1(end,:)/2 + x1(end-1,:)/2;
    x2_new(:,1) = x2(:,1)/2 + x2(:,2)/2;
    x2_new(:,end) = x2(:,end)/2 + x2(:,end-1)/2;
end