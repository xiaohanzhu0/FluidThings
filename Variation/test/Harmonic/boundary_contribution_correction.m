% This function calculates the boundary contributions to the points on
% boundaryies. So the output will be the complete matrix
% system, assuming the input rows of A and b corresponding to boundary 
% points are cleared as zero.

function [A,b] = boundary_contribution_correction(x1,x2,A,b1,b2)
[Nx2, Nx1] = size(x1);
N = Nx2*Nx1;
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

b1(:,1) = b1(:,1) + x1(:,1).*n_left(1,:)' + x2(:,1).*n_left(2,:)';
b1(:,end) = b1(:,end) + x1(:,end).*n_right(1,:)' + x2(:,end).*n_right(2,:)';
b1(1,:) = b1(1,:) + x1(1,:).*n_bottom(1,:) + x2(1,:).*n_bottom(2,:);
b1(end,:) = b1(end,:) + x1(end,:).*n_top(1,:) + x2(end,:).*n_top(2,:);
b2(:,1) = 0; b2(:,end) = 0; b2(1,:) = 0; b2(end,:) = 0;

b1(1,1) = x1(1,1); b1(end,1) = x1(end,1); b1(1,end) = x1(1,end); b1(end,end) = x1(end,end);
b2(1,1) = x2(1,1); b2(end,1) = x2(end,1); b2(1,end) = x2(1,end); b2(end,end) = x2(end,end);

b = [b1(:); b2(:)];
end

