function boundary_points = extract_boundary_points(x1, x2)
boundary_points.b = [x1(1,:); x2(1,:)];
boundary_points.t = [x1(end,:); x2(end,:)];
boundary_points.l = [x1(:,1), x2(:,1)];
boundary_points.r = [x1(:,end), x2(:,end)];
end