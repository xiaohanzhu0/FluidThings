clear
close all
addpath('./','./utils')

cf.problem = 5;
cf.Nx1 = 217*2+1;
cf.Nx2 = 71*1;
cf.Nx1 = 41;
cf.Nx2 = 41;
cf.N = cf.Nx1*cf.Nx2;
cf.alpha = 1.005;

cf.sigma1 = 1;
cf.sigma2 = 1;

cf.nonlinear = 7;
cf.fixed_bc = 0;
cf.omega = 1;
cf.offdiag = 1;
cf.max_iter = 500;
cf.tolerance = 1e-6;
cf.smooth = 0;

cf.append_trail = 0;
cf.new_airfoil = 1;

cf.animation = 1;
cf.iter_per_frame = 1;
cf.plot_res = 1;
cf.pause_time = 0;
cf.make_gif = 0;
cf.gif_name = 'failed_example.gif';
cf.title_name = 'Alternative';
cf.save_output = 0;

if cf.problem == 1 || cf.problem == 2
    cf.C = 1;
end

if cf.nonlinear == 7
    cf.exact = 0;
    cf.repulse = 0;
end

if cf.new_airfoil == 1
    cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
    cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
end

%%
Nx1 = 41;
Nx2 = 41;
problem = 5;
alpha = 1.005;
append_trail = 0;
new_airfoil = 1;
C = 1;

N = Nx1*Nx2;
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);

if ismember(problem, [1,2])
    M_type = problem;
    [x1, x2] = meshgrid(s1, s2);
    [x1_exact, x2_exact]  = meshgrid(s1, s2);
    x1 = [0,1;0,1];
    x2 = [0,0;1,1];
elseif problem == 3
    M_type = problem;
    [x1, x2] = InitProb3(Nx1, Nx2, 0.1);
    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,:); x2(end,:)];
    boundary_points.l = [x1(:,1), x2(:,1)];
    boundary_points.r = [x1(:,end), x2(:,end)];
elseif problem == 4
    M_type = problem;
    [x1, x2] = InitProb4(Nx1, Nx2);
    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,1), x1(end,end); x2(end,1), x2(end,end)];
    boundary_points.l = [x1(1,1), x1(end,1); x2(1,1), x2(end,1)]';
    boundary_points.r = [x1(1,end), x1(end,end); x2(1,end), x2(end,end)]';
elseif problem == 5
    [x1, x2, M_type] = InitProb5(cf);
    x2(1,1) = x2(1,1) + 1e-2;
    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,:); x2(end,:)];
    boundary_points.l = [x1(:,1), x2(:,1)];
    boundary_points.r = [x1(:,end), x2(:,end)];
    Nx1 = size(x1,2);
    Nx2 = size(x1,1);
    N = Nx1 * Nx2;
elseif problem == 6
    M_type = problem;
    [x1_temp, x2_temp] = meshgrid(s1, s2);
    theta = pi/6;
    x1 = x1_temp*cos(theta) - x2_temp*sin(theta);
    x2 = x1_temp*sin(theta) + x2_temp*cos(theta);

    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,:); x2(end,:)];
    boundary_points.l = [x1(:,1), x2(:,1)];
    boundary_points.r = [x1(:,end), x2(:,end)];
elseif problem == 7
    M_type = problem;
    [x1_temp, x2_temp] = meshgrid(s1, s2);
    x1 = x1_temp + x2_temp;
    x2 = x2_temp*2;
    boundary_points.b = [x1(1,:); x2(1,:)];
    boundary_points.t = [x1(end,:); x2(end,:)];
    boundary_points.l = [x1(:,1), x2(:,1)];
    boundary_points.r = [x1(:,end), x2(:,end)];
end


boundary_points.g2 = [x1(1,:); x2(1,:)]';
boundary_points.g4 = flip([x1(end,:); x2(end,:)]',1);
boundary_points.g1 = flip([x1(:,1), x2(:,1)],1);
boundary_points.g3 = [x1(:,end), x2(:,end)];
%%
hgrad = 1;
[p, t, model] = triangulate_from_edgeArcs(boundary_points,'Hmax',0.01,'hgrad',1.05);
figure
pdemesh(model);

%%
if problem == 1
    Mfun = @(x1,x2) [40000 * (1 + 15 * x1).^(-2), 0;
                     0, 40000 * (1 + 15 * x2).^(-2)];
    bc.x0 = 2; bc.x1 = 3; bc.y0 = 1; bc.y1 = 4;
elseif problem == 2
    Mfun = @(x1,x2) [1000 + 600*sin(2*pi*x1).*sin(2*pi*x2), 0;
                     0, 1000 - 600*sin(2*pi*x1).*sin(2*pi*x2)];
    bc.x0 = 2; bc.x1 = 3; bc.y0 = 1; bc.y1 = 4;
elseif problem == 3 
    Mfun = @(x1,x2) [1,0;0,1];
    bc.x0 = 2; bc.x1 = 3; bc.y0 = 1; bc.y1 = 4;
elseif problem == 4
    Mfun = @(x1,x2) [1,0;0,1];
    bc.x0 = 2; bc.x1 = 6; bc.y0 = [1,3,4]; bc.y1 = 5;
elseif problem == 5
    %Mfun = @(x1,x2) [1,0;0,1];
    Mfun = @(x1,x2) [M_type.F11(x1,x2), M_type.F12(x1,x2);
                     M_type.F12(x1,x2), M_type.F22(x1,x2)];
    %bc.x0 = 4; bc.x1 = 2; bc.y0 = 1; bc.y1 = 3;
    bc.x0 = 3; bc.x1 = 1; bc.y0 = 2; bc.y1 = 4;
elseif problem == 6
    
end

[s1, s2, model1, model2] = solve_inverse_map_pdetbx(p, t, boundary_points, Mfun, bc);
figure
pdegplot(model1,VertexLabels="on", ...
           EdgeLabels="on", ...
           FaceLabels="on")
%%
figure
pdeplot(model1.Mesh, XYData=s1, ColorMap="jet");
figure
pdeplot(model2.Mesh, XYData=s2,  ColorMap="jet");

%%
% Given: p,t, and FEM solutions s1,s2 at nodes p
N1 = 160; N2 = 80;                     % structured grid size in s-space
[Xs, Ys, S1g, S2g, info] = inverse_map_to_grid(p, t, s1, s2, N1, N2);

% Visualize the structured grid in physical space:
figure(10); hold on;
plot(Xs, Ys, '-r'); hold on; plot(Xs', Ys', '-r'); axis equal; hold on
title('Structured grid x(s) recovered from s(x)');

[L1, L2, Linf, gtilde2, sigma1, sigma2] = Cost(Xs, Ys, M_type, C);
[L, gtilde2, sigma1, sigma2] = CostExact(Xs, Ys, M_type, C);