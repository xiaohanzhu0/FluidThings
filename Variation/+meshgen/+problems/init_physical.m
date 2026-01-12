function [x1, x2, boundary_points, cf] = init_physical(cf)
%INIT_PHYSICAL Initialize physical grid and boundary points for EL formulation.
    boundary_points = struct();
    s1 = linspace(0, 1, cf.Nx1);
    s2 = linspace(0, 1, cf.Nx2);

    if ismember(cf.problem, [1, 2])
        [x1, x2] = meshgrid(s1, s2);
    elseif cf.problem == 3
        [x1, x2] = problems.InitProb3(cf.Nx1, cf.Nx2);
        boundary_points.b = [x1(1,:); x2(1,:)];
        boundary_points.t = [x1(end,:); x2(end,:)];
        boundary_points.l = [x1(:,1), x2(:,1)];
        boundary_points.r = [x1(:,end), x2(:,end)];
    elseif cf.problem == 4
        [x1, x2] = problems.InitProb4(cf.Nx1, cf.Nx2);
        boundary_points.b = [x1(1,:); x2(1,:)];
        boundary_points.t = [x1(end,1), x1(end,end); x2(end,1), x2(end,end)];
        boundary_points.l = [x1(1,1), x1(end,1); x2(1,1), x2(end,1)]';
        boundary_points.r = [x1(1,end), x1(end,end); x2(1,end), x2(end,end)]';
    elseif cf.problem == 5
        [x1, x2] = problems.InitProb5(cf);
        boundary_points.b = [x1(1,:); x2(1,:)];
        boundary_points.t = [x1(end,:); x2(end,:)];
        boundary_points.l = [x1(:,1), x2(:,1)];
        boundary_points.r = [x1(:,end), x2(:,end)];
        cf.Nx1 = size(x1,2);
        cf.Nx2 = size(x1,1);
    elseif cf.problem == 6
        [x1_temp, x2_temp] = meshgrid(s1, s2);
        theta = pi/6;
        x1 = x1_temp*cos(theta) - x2_temp*sin(theta);
        x2 = x1_temp*sin(theta) + x2_temp*cos(theta);

        boundary_points.b = [x1(1,:); x2(1,:)];
        boundary_points.t = [x1(end,:); x2(end,:)];
        boundary_points.l = [x1(:,1), x2(:,1)];
        boundary_points.r = [x1(:,end), x2(:,end)];
    elseif cf.problem == 7
        [x1_temp, x2_temp] = meshgrid(s1, s2);
        x1 = x1_temp + x2_temp;
        x2 = x2_temp*2;
        boundary_points.b = [x1(1,:); x2(1,:)];
        boundary_points.t = [x1(end,:); x2(end,:)];
        boundary_points.l = [x1(:,1), x2(:,1)];
        boundary_points.r = [x1(:,end), x2(:,end)];
    elseif cf.problem == 8
        [x1, x2] = problems.InitProb8(cf.Nx1, cf.Nx2);
        boundary_points.b = [x1(1,:); x2(1,:)];
        boundary_points.t = [x1(end,:); x2(end,:)];
        boundary_points.l = [x1(:,1), x2(:,1)];
        boundary_points.r = [x1(:,end), x2(:,end)];
        cf.Nx1 = size(x1,2);
        cf.Nx2 = size(x1,1);
    else
        error('Unsupported problem id: %d', cf.problem);
    end
end
