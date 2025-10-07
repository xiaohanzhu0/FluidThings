function [A, b, res] = AssembleLinearSystemConserve(x1, x2, Mfun, param)
[Nx2, Nx1] = size(x1);
N = Nx1*Nx2;
N_all = 2*N;

ds1 = 1 / Nx1;
ds2 = 1 / Nx2;
sigma1 = 1 / Nx1; %param.sigma1;
sigma2 = 1 / Nx2; %param.sigma2;

M11 = zeros(Nx2+2,Nx1+2); M12 = zeros(Nx2+2,Nx1+2); M22 = zeros(Nx2+2,Nx1+2);

M = Mfun(x1,x2);
dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;
dM12dx1=M.dM12dx1; dM12dx2=M.dM12dx2;
[dx1ds1, dx2ds2] = DCentral(x1, x2, ds1, ds2);
[dx2ds1, dx1ds2] = DCentral(x2, x1, ds1, ds2);

m1 = zeros(Nx2+2,Nx1+2); m2 = zeros(Nx2+2,Nx1+2); 
m1(2:end-1,2:end-1) = sigma1^2*(M.M11.*dx1ds1.^2 + 2*M.M12.*dx1ds1.*dx2ds1 + M.M22.*dx2ds1.^2) - 1;
m2(2:end-1,2:end-1) = sigma2^2*(M.M11.*dx1ds2.^2 + 2*M.M12.*dx1ds2.*dx2ds2 + M.M22.*dx2ds2.^2) - 1;


M11(2:end-1,2:end-1) = M.M11; 
M12(2:end-1,2:end-1) = M.M12; 
M22(2:end-1,2:end-1) = M.M22;

M11_f1 = (M11(:,1:end-1) + M11(:,2:end)) / 2;
M11_f2 = (M11(1:end-1,:) + M11(2:end,:)) / 2;
M12_f1 = (M12(:,1:end-1) + M12(:,2:end)) / 2;
M12_f2 = (M12(1:end-1,:) + M12(2:end,:)) / 2;
M22_f1 = (M22(:,1:end-1) + M22(:,2:end)) / 2;
M22_f2 = (M22(1:end-1,:) + M22(2:end,:)) / 2;

%M11_f1 = 1 ./ ((1./M11(:,1:end-1) + 1./M11(:,2:end)) / 2);
%M11_f2 = 1 ./ ((1./M11(1:end-1,:) + 1./M11(2:end,:)) / 2);
%M12_f1 = 1 ./ ((1./M12(:,1:end-1) + 1./M12(:,2:end)) / 2);
%M12_f2 = 1 ./ ((1./M12(1:end-1,:) + 1./M12(2:end,:)) / 2);
%M22_f1 = 1 ./ ((1./M22(:,1:end-1) + 1./M22(:,2:end)) / 2);
%M22_f2 = 1 ./ ((1./M22(1:end-1,:) + 1./M22(2:end,:)) / 2);
if param.exact == 1
    q1 = max(abs(m1),[],'all') + 1;
    q2 = max(abs(m2),[],'all') + 1;
    %q1 = 1.01;
    %q2 = 1.01;
    M11_f1 = M11_f1.*(m1(:,1:end-1)/2 + m1(:,2:end)/2 + q1);
    M11_f2 = M11_f2.*(m2(1:end-1,:)/2 + m2(2:end,:)/2 + q2);
    M12_f1 = M12_f1.*(m1(:,1:end-1)/2 + m1(:,2:end)/2 + q1);
    M12_f2 = M12_f2.*(m2(1:end-1,:)/2 + m2(2:end,:)/2 + q2);
    M22_f1 = M22_f1.*(m1(:,1:end-1)/2 + m1(:,2:end)/2 + q1);
    M22_f2 = M22_f2.*(m2(1:end-1,:)/2 + m2(2:end,:)/2 + q2);
end

%{
if param.repulse == 1
    %p1 = M.M11.*dx1ds1.^2 + 2*M.M12.*dx1ds1.*dx2ds1 + M.M22.*dx2ds1.^2;
    %p2 = M.M11.*dx1ds2.^2 + 2*M.M12.*dx1ds2.*dx2ds2 + M.M22.*dx2ds2.^2;
    p1 = dx1ds1.^2 + dx2ds1.^2;
    p2 = dx1ds2.^2 + dx2ds2.^2;
    %p1_f1 = (p1(:,1:end-1) + p1(:,2:end)) / 2;
    %p1_f2 = (p1(1:end-1,:) + p1(2:end,:)) / 2;
    %p2_f1 = (p2(:,1:end-1) + p2(:,2:end)) / 2;
    %p2_f1 = (p2(1:end-1,:) + p2(2:end,:)) / 2;

    p1_f1 = 2 ./ (p1(:,1:end-1).^2 + p1(:,2:end).^2);
    p1_f2 = 2 ./ (p1(1:end-1,:).^2 + p1(2:end,:).^2);
    p2_f1 = 2 ./ (p2(:,1:end-1).^2 + p2(:,2:end).^2);
    p2_f2 = 2 ./ (p2(1:end-1,:).^2 + p2(2:end,:).^2);

    lambda = 1e-1;
    M11_f1(2:end-1,2:end-1) = M11_f1(2:end-1,2:end-1) - lambda*p1_f1;
    M11_f2(2:end-1,2:end-1) = M11_f2(2:end-1,2:end-1) - lambda*p2_f2;
    %M12_f1(2:end-1,2:end-1) = M12_f1(2:end-1,2:end-1) - lambda*p1_f1;
    %M12_f2(2:end-1,2:end-1) = M12_f2(2:end-1,2:end-1) - lambda*p2_f2;
    M22_f1(2:end-1,2:end-1) = M22_f1(2:end-1,2:end-1) - lambda*p1_f1;
    M22_f2(2:end-1,2:end-1) = M22_f2(2:end-1,2:end-1) - lambda*p2_f2;
end
%}

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
Ux2_11 = spdiags(e, Nx2, N, N) .* Ux2_11_coef(:);

Lx2_11_coef = M11_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_11 = spdiags(e, -Nx2, N, N) .* Lx2_11_coef(:);

U1_12_coef = M12_f2(2:end,2:end-1) * (sigma2/ds2)^2;
U1_12 = spdiags(e, 1, N, N) .* U1_12_coef(:);

L1_12_coef = M12_f2(1:end-1,2:end-1)  * (sigma2/ds2)^2;
L1_12 = spdiags(e, -1, N, N) .* L1_12_coef(:);

Ux2_12_coef = M12_f1(2:end-1,2:end) * (sigma1/ds1)^2;
Ux2_12 = spdiags(e, Nx2, N, N) .* Ux2_12_coef(:);

Lx2_12_coef = M12_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_12 = spdiags(e, -Nx2, N, N) .* Lx2_12_coef(:);

A = sparse(N_all,N_all);
A(1:N,1:N) = D11 + U1_11 + L1_11 + Ux2_11 + Lx2_11;
A(1:N,N+1:end) = D12 + U1_12 + L1_12 + Ux2_12 + Lx2_12;


%%

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
Ux2_22 = spdiags(e, Nx2, N, N) .* Ux2_22_coef(:);

Lx2_22_coef = M22_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_22 = spdiags(e, -Nx2, N, N) .* Lx2_22_coef(:);

U1_21_coef = M12_f2(2:end,2:end-1) * (sigma2/ds2)^2;
U1_21 = spdiags(e, 1, N, N) .* U1_21_coef(:);

L1_21_coef = M12_f2(1:end-1,2:end-1)  * (sigma2/ds2)^2;
L1_21 = spdiags(e, -1, N, N) .* L1_21_coef(:);

Ux2_21_coef = M12_f1(2:end-1,2:end) * (sigma1/ds1)^2;
Ux2_21 = spdiags(e, Nx2, N, N) .* Ux2_21_coef(:);

Lx2_21_coef = M12_f1(2:end-1,1:end-1) * (sigma1/ds1)^2;
Lx2_21 = spdiags(e, -Nx2, N, N) .* Lx2_21_coef(:);

A(N+1:end,N+1:end) = D22 + U1_22 + L1_22 + Ux2_22 + Lx2_22;
A(N+1:end,1:N) = D21 + U1_21 + L1_21 + Ux2_21 + Lx2_21;
A = 2*A;

%% Apply boundary conditions
id = GetIndex(Nx1, Nx2);
   
if param.fixed_bc == 0
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
else
A(id.l,:) = 0; A(id.l+N,:) = 0;
A(id.r,:) = 0; A(id.r+N,:) = 0;
A(id.b,:) = 0; A(id.b+N,:) = 0;
A(id.t,:) = 0; A(id.t+N,:) = 0;

A(id.l,id.l) = eye(length(id.l));
A(id.r,id.r) = eye(length(id.r));
A(id.b,id.b) = eye(length(id.b));
A(id.t,id.t) = eye(length(id.t));
A(id.l+N,id.l+N) = eye(length(id.l));
A(id.r+N,id.r+N) = eye(length(id.r));
A(id.b+N,id.b+N) = eye(length(id.b));
A(id.t+N,id.t+N) = eye(length(id.t));
end


for i = id.corner; A(i, :) = 0; A(i, i) = 1; end
for i = N+id.corner; A(i, :) = 0; A(i, i) = 1; end

%% Assemble the interior vector b
if param.exact == 1
    m1 = m1(2:end-1,2:end-1);
    m2 = m2(2:end-1,2:end-1);
else
    q1 = 0;
    q2 = 0;
    m1 = ones(Nx2,Nx1);
    m2 = ones(Nx2,Nx1);
end
b1 = zeros(Nx2, Nx1);
b2 = zeros(Nx2, Nx1);

b1 = b1 + dM22dx1.*((m1+q1).*dx2ds1.^2*sigma1^2 + (m2+q2).*dx2ds2.^2*sigma2^2);
b1 = b1 + 2*dM12dx1.*((m1+q1).*dx1ds1.*dx2ds1*sigma1^2 + (m2+q2).*dx1ds2.*dx2ds2*sigma2^2);
b1 = b1 + dM11dx1.*((m1+q1).*dx1ds1.^2*sigma1^2 + (m2+q2).*dx1ds2.^2*sigma2^2);

b2 = b2 + dM22dx2.*((m1+q1).*dx2ds1.^2*sigma1^2 + (m2+q2).*dx2ds2.^2*sigma2^2);
b2 = b2 + 2*dM12dx2.*((m1+q1).*dx1ds1.*dx2ds1*sigma1^2 + (m2+q2).*dx1ds2.*dx2ds2*sigma2^2);
b2 = b2 + dM11dx2.*((m1+q1).*dx1ds1.^2*sigma1^2 + (m2+q2).*dx1ds2.^2*sigma2^2);

if param.repulse == 1
    lambda = 1;
    p1 = M.M11.*dx1ds1.^2 + 2*M.M12.*dx1ds1.*dx2ds1 + M.M22.*dx2ds1.^2;
    p2 = M.M11.*dx1ds2.^2 + 2*M.M12.*dx1ds2.*dx2ds2 + M.M22.*dx2ds2.^2;
    %p1 = dx1ds1.^2 + dx2ds1.^2;
    %p2 = dx1ds2.^2 + dx2ds2.^2;

    aux1 = sigma1*p1.^(-2).*(M.M11.*dx1ds1 + M.M12.*dx2ds1) + sigma2*p2.^(-2).*(M.M11.*dx1ds2 + M.M12.*dx2ds2);
    aux2 = sigma1*p1.^(-2).*(M.M12.*dx1ds1 + M.M22.*dx2ds1) + sigma2*p2.^(-2).*(M.M12.*dx1ds2 + M.M22.*dx2ds2);
    %aux1 = sigma1*p1.^(-2).*(dx1ds1) + sigma2*p2.^(-2).*(dx1ds2);
    %aux2 = sigma1*p1.^(-2).*(dx2ds1) + sigma2*p2.^(-2).*(dx2ds2);

    [daux1ds1, daux1ds2] = DCentral(aux1, aux1, ds1, ds2);
    [daux2ds1, daux2ds2] = DCentral(aux2, aux2, ds1, ds2);

    b1 = b1 + lambda*(daux1ds1 + daux1ds2);
    b2 = b2 + lambda*(daux2ds1 + daux2ds2);
end

if param.fixed_bc == 0
    b1(:,1) = x1(:,1).*n_left(1,:)' + x2(:,1).*n_left(2,:)';
    b1(:,end) = x1(:,end).*n_right(1,:)' + x2(:,end).*n_right(2,:)';
    b1(1,:) = x1(1,:).*n_bottom(1,:) + x2(1,:).*n_bottom(2,:);
    b1(end,:) = x1(end,:).*n_top(1,:) + x2(end,:).*n_top(2,:);
    b2(:,1) = 0; b2(:,end) = 0; b2(1,:) = 0; b2(end,:) = 0;
else
    b1(:,1) = x1(:,1);
    b1(:,end) = x1(:,end);
    b1(1,:) = x1(1,:);
    b1(end,:) = x1(end,:);
    b2(:,1) = x2(:,1);
    b2(:,end) = x2(:,end);
    b2(1,:) = x2(1,:);
    b2(end,:) = x2(end,:);
end

b = [b1(:); b2(:)];
res = A*[x1(:); x2(:)] - b;
%%
perm = zeros(N_all,1);
perm(1:2:end) =        1:N;      % odd entries get u‑indices
perm(2:2:end) = N + (1:N);

A_p = A(perm, perm);

end