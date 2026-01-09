clear

cf.block_idx = 4;
cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
[x1, x2, M_samp, Mfun] = problems.InitProb9(cf);

%cf.new_airfoil = 1; cf.Nx1 = 151; cf.Nx2 = 51; cf.alpha = 0.05; cf.append_trail = 0;
%cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
%cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
%[x1, x2, M_samp, Mfun] = problems.InitProb5(cf);

x1 = M_samp.x_metric;
x2 = M_samp.y_metric;

figure
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
%[Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1, x2);

boundary_points = extract_boundary_points(x1,x2);
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

t1 = linspace(0,1,27);
t2 = linspace(0,1,37);
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
%%

t1 = linspace(0,1,27);
t2 = linspace(0,1,37);
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

analysis.plot_metric(x1, x2, M);
analysis.plot_metric(t1', t2', Mp);
%% Solve inverse instead
detJ = Mp.M11.*Mp.M22 - Mp.M12.^2;
Maux.M11 = Mp.M22 ./ detJ;
Maux.M22 = Mp.M11 ./ detJ;
Maux.M12 = -Mp.M12 ./ detJ;
Mp = Maux;
%%
for i=1:5
    Mp.M11(:,1) = Mp.M11(:,2); Mp.M11(:,end) = Mp.M11(:,end-1);
    Mp.M12(:,1) = Mp.M12(:,2); Mp.M12(:,end) = Mp.M12(:,end-1);
    Mp.M22(1,:) = Mp.M22(2,:); Mp.M22(end,:) = Mp.M22(end-1,:);

    Mp.M11(:,2:end-1) = (Mp.M11(:,1:end-2) + 2*Mp.M11(:,2:end-1) + Mp.M11(:,3:end)) / 4;
    Mp.M12(:,2:end-1) = (Mp.M12(:,1:end-2) + 2*Mp.M12(:,2:end-1) + Mp.M12(:,3:end)) / 4;
    Mp.M22(:,2:end-1) = (Mp.M22(:,1:end-2) + 2*Mp.M22(:,2:end-1) + Mp.M22(:,3:end)) / 4;
    Mp.M11(2:end-1,:) = (Mp.M11(1:end-2,:) + 2*Mp.M11(2:end-1,:) + Mp.M11(3:end,:)) / 4;
    Mp.M12(2:end-1,:) = (Mp.M12(1:end-2,:) + 2*Mp.M12(2:end-1,:) + Mp.M12(3:end,:)) / 4;
    Mp.M22(2:end-1,:) = (Mp.M22(1:end-2,:) + 2*Mp.M22(2:end-1,:) + Mp.M22(3:end,:)) / 4;
end

Mp.M11 = Mp.M11 + 0.01*max(Mp.M11,[],'all');
Mp.M22 = Mp.M22 + 0.01*max(Mp.M22,[],'all');

Mp11 = griddedInterpolant(t1, t2, Mp.M11', 'cubic');
Mp22 = griddedInterpolant(t1, t2, Mp.M22', 'cubic');
Mp12 = griddedInterpolant(t1, t2, Mp.M12', 'cubic');
%%

Nt1 = 750; Nt2 = 370;
N = Nt1*Nt2;
S1 = linspace(0,1,Nt1);
S2 = linspace(0,1,Nt2);
[S1, S2] = meshgrid(S1, S2);
t1 = linspace(0,1,Nt1);
t2 = linspace(0,1,Nt2);
[t1, t2] = meshgrid(t1, t2);



%%
[s1, s2] = AssembleLinearSystemDual(t1, t2, Mp11, Mp12, Mp22);

s1_samp = linspace(0,1,27);
s2_samp = linspace(0,1,37);
[s1_samp,s2_samp] = meshgrid(s1_samp,s2_samp);
t1_samp = griddata(s1,s2,t1,s1_samp,s2_samp);
t2_samp = griddata(s1,s2,t2,s1_samp,s2_samp);

figure
plot(s1, s2, 'k'); hold on; plot(s1', s2', 'k'); hold off
figure
plot(t1_samp, t2_samp, 'k'); hold on; plot(t1_samp', t2_samp', 'k'); hold off
figure
plot(x1_in_t(t1_samp,t2_samp), x2_in_t(t1_samp,t2_samp), 'k'); hold on; plot(x1_in_t(t1_samp,t2_samp)', x2_in_t(t1_samp,t2_samp)', 'k'); hold off
