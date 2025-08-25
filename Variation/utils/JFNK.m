function [dx, Res] = JFNK(x1, x2, param)
    % Make preconditioner
    x = [x1(:); x2(:)];
    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;
    perm = zeros(N_all,1);
    perm(1:2:end) =        1:N;      % odd entries get u‑indices
    perm(2:2:end) = N + (1:N);      % even entries get v‑indices
    invperm = zeros(size(perm));
    invperm(perm) = 1:numel(perm);

    %x = x(perm);
    Res = residual(x, param, perm, invperm);

    e = ones(N,1);
    L = spdiags([e, e, -4*e, e, e], [-Nx2, -1, 0, 1, Nx2], N, N);
    
    A = sparse(N_all, N_all);
    A(1:N,1:N) = A(1:N,1:N) - 2*L;
    A(N+1:end,N+1:end) = A(N+1:end,N+1:end) - 2*L;
    A(1:N,N+1:end) = A(1:N,N+1:end) - 2*L;
    A(N+1:end,1:N) = A(N+1:end,1:N) - 2*L;
    %A = A(perm, perm);

    alpha = max(sum(abs(A),2)./diag(A))-2;
    L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',alpha));

    restart = 10;
    tolKrylov = 1e-3;
    epsFD = sqrt(eps) * max(1, norm(x)) * 1000 ;
    Afun = @(v) jtv(x, v, Res, epsFD, param, perm, invperm);
    [dx,flag,relres,iter,resvec] = gmres(Afun, -Res, restart, tolKrylov, 200, L, L');
    dx = dx(invperm);
end


function jv = jtv(x, v, Res0, epsFD, param, perm, invperm)
    Res1 = residual(x+epsFD*v, param, perm, invperm);
    jv = (Res1 - Res0) / epsFD;
end

function Res = residual(x, param, perm, invperm)
    %x = x(invperm);
    x1 = reshape(x(1:param.N), param.Nx2, param.Nx1);
    x2 = reshape(x(param.N+1:end), param.Nx2, param.Nx1);
    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    M = GetM(x1, x2, param.problem, param.C);
    M11=M.M11; M22=M.M22; 
    dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;

    M12=M.M12;
    dM12dx1=M.dM12dx1; dM12dx2=M.dM12dx2;

    [dx1ds1, dx2ds2] = DCentral(x1, x2, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, sigma1, sigma2);

    param.dxi = sigma1; param.deta = sigma2;
    ddx1dds1 = CentralD2(x1, 'horizontal', param);
    ddx1dds2 = CentralD2(x1, 'vertical', param);
    ddx2dds1 = CentralD2(x2, 'horizontal', param);
    ddx2dds2 = CentralD2(x2, 'vertical', param);
    id = GetIndex(Nx1, Nx2);

        
    t_bottom = GetBoundaryTangent(x1(1,:), x2(1,:), 1);
    t_top = GetBoundaryTangent(x1(end,:), x2(end,:), 1);
    t_left = GetBoundaryTangent(x1(:,1), x2(:,1), 1);
    t_right = GetBoundaryTangent(x1(:,end), x2(:,end), 1);

    n_bottom = [-t_bottom(2,:); t_bottom(1,:)];
    n_top = [-t_top(2,:); t_top(1,:)];
    n_left = [-t_left(2,:); t_left(1,:)];
    n_right = [-t_right(2,:); t_right(1,:)];

    e = ones(N,1);
    L = spdiags([e, e, -4*e, e, e], [-Nx2, -1, 0, 1, Nx2], N, N);
    A = sparse(N_all, N_all);
    A(1:N,1:N) = A(1:N,1:N) - 2*L.*M11(:);
    A(N+1:end,N+1:end) = A(N+1:end,N+1:end) - 2*L.*M22(:);
    A(1:N,N+1:end) = A(1:N,N+1:end) - 2*L.*M12(:);
    A(N+1:end,1:N) = A(N+1:end,1:N) - 2*L.*M12(:);

    % Boundary setup
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
    
    b1 = dM11dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
    b1 = b1 + 2*dM12dx2.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
    b1 = b1 + 2*dM11dx1.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
    b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);

    b2 = dM22dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
    b2 = b2 + 2*dM12dx1.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
    b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);

    % RHS changes due to Dirichlet on boundary
    b1(:,1) = 0;
    b1(:,end) = 0;
    b1(1,:) = 0;
    b1(end,:) = 0;

    b1(:,1) = x1(:,1).*n_left(1,:)' + x2(:,1).*n_left(2,:)';
    b1(:,end) = x1(:,end).*n_right(1,:)' + x2(:,end).*n_right(2,:)';
    b1(1,:) = x1(1,:).*n_bottom(1,:) + x2(1,:).*n_bottom(2,:);
    b1(end,:) = x1(end,:).*n_top(1,:) + x2(end,:).*n_top(2,:);
    b2(:,1) = 0; b2(:,end) = 0; b2(1,:) = 0; b2(end,:) = 0;

    % Fix corners
    for i = id.corner; A(i, :) = 0; A(i, i) = 1; end
    for i = N+id.corner; A(i, :) = 0; A(i, i) = 1; end
    b1(1,1) = x1(1,1); b1(1,end) = x1(1,end); b1(end,1) = x1(end,1); b1(end,end) = x1(end,end);
    b2(1,1) = x2(1,1); b2(1,end) = x2(1,end); b2(end,1) = x2(end,1); b2(end,end) = x2(end,end);

    if param.smooth ~= 0
        e = ones(N,1);
        % stride in i-direction is Nx2
        D1 = spdiags([-e, e], [0, Nx2], N, N) / sigma1;
        % stride in j-direction is 1
        D2 = spdiags([-e, e], [0, 1   ], N, N) / sigma2;
        A_lap = D1' * D1  +  D2' * D2;
        A_lap = blkdiag(A_lap,A_lap);
        A = A + param.smooth*A_lap;
    end
    

    if param.smooth ~= 0
        b = [b1(:); b2(:)];
        b = b - param.smooth*A_lap*[x1(:); x2(:)];
        b1 = reshape(b(1:N), Nx2, Nx1);
        b2 = reshape(b(N+1:end), Nx2, Nx1);
    end

    b = [b1(:); b2(:)];
    x = [x1(:); x2(:)];
    Res = A*x - b;
    %Res = Res(perm);
end