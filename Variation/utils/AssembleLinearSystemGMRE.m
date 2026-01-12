function [dx, Res] = AssembleLinearSystemGMRE(x1, x2, Mfun, param)

    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    M = Mfun(x1,x2);
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

    perm = zeros(N_all,1);
    perm(1:2:end) =        1:N;      % odd entries get u‑indices
    perm(2:2:end) = N + (1:N);      % even entries get v‑indices
    invperm = zeros(size(perm));
    invperm(perm) = 1:numel(perm);

    e = ones(N,1);
    L = spdiags([e, e, -4*e, e, e], [-Nx2, -1, 0, 1, Nx2], N, N);
    %%
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
    
    %M = M_old;
    %M = GetM(x1, x2, param.problem, param.C);
    dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;
    dM12dx1=M.dM12dx1; dM12dx2=M.dM12dx2;
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
    % RHS changes due to Neumann on boundary
    %{
    b2(:,1) = -x1(:,1).*t_left(1,:)' - x2(:,1).*t_left(2,:)' ...
              +1/3*( 4*(x1(:,2).*t_left(1,:)'+x2(:,2).*t_left(2,:)') ...
                      -(x1(:,3).*t_left(1,:)'+x2(:,3).*t_left(2,:)') );

    b2(:,end) = -x1(:,end).*t_right(1,:)' - x2(:,end).*t_right(2,:)' ...
              +1/3*( 4*(x1(:,end-1).*t_right(1,:)'+x2(:,end-1).*t_right(2,:)') ...
                      -(x1(:,end-2).*t_right(1,:)'+x2(:,end-2).*t_right(2,:)') );

    b2(1,:) = -x1(1,:).*t_bottom(1,:) - x2(1,:).*t_bottom(2,:) ...
              +1/3*( 4*(x1(2,:).*t_bottom(1,:)+x2(2,:).*t_bottom(2,:)) ...
                      -(x1(3,:).*t_bottom(1,:)+x2(3,:).*t_bottom(2,:)) );

    b2(end,:) = -x1(end,:).*t_top(1,:) - x2(end,:).*t_top(2,:) ...
              +1/3*( 4*(x1(end-1,:).*t_top(1,:)+x2(end-1,:).*t_top(2,:)) ...
                      -(x1(end-2,:).*t_top(1,:)+x2(end-2,:).*t_top(2,:)) );
    %}

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

    %{
    for k = 1:max_newton_iter
        % 1. Calculate the residual F(u) for the current solution u.
        %    F(u) = A*u - lambda*h^2*exp(u), where A is the discrete Laplacian.
        F = compute_residual(u, n, lambda, h);
    
        % 2. Check for convergence.
        normF = norm(F);
        fprintf('%-4d | %-14.6e |', k, normF);
    
        if normF < newton_tol
            fprintf('\nConvergence achieved!\n');
            break;
        end
    
        % 3. Solve the linear system J*s = -F using GMRES.
        %    J is the Jacobian of F. We don't form J explicitly. Instead, we
        %    provide a function handle (@(v) Jv_product(...)) that computes
        %    the Jacobian-vector product J*v.
        Jv_handle = @(v) Jv_product(v, u, n, lambda, h);
    
        [s, gmres_flag, gmres_relres, gmres_iter] = gmres(Jv_handle, -F, ...
            gmres_restart, gmres_tol, max_gmres_iter, L, L');
    
        fprintf(' %-10d | %-15.4e | %d(%d)\n', ...
            gmres_flag, gmres_relres, gmres_iter(1), gmres_iter(2));
    
        % 4. Update the solution.
        u = u + s;
    end
    %}


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

    %A = A(perm,perm);
    %b = b(perm);
    %x = x(perm);

    Res = A*x - b;
    %Res = Res(invperm);

    restart = 10;
    tolKrylov = 1e-4;
    epsFD = sqrt(eps) * max(1, norm(x));

    Afun = @(v) jtv(x, v, M, Res, epsFD, param, perm, invperm);
    opts.type  = 'ilutp'; 
    opts.udiag = 1;
    [L,U] = ilu(A, opts);
    
    [dx,flag,relres,iter,resvec] = gmres(Afun, -Res, restart, tolKrylov, 500, L, U);
    
    % ---- backtracking line search ----
    alpha = 1;
    c1  = 1e-4;   % Armijo constant
    rho = 0.5;    % backtracking factor
    min_alpha = 1e-3;
    while alpha > min_alpha
      x_trial = x + alpha*dx;
      Res_trial = residual(x_trial, M, param, perm, invperm);
      if norm(Res_trial) <= (1 - c1*alpha)*norm(Res)
        break;
      end
      alpha = rho * alpha;
    end
    if alpha <= min_alpha
      warning("Line‐search: step size got very small (α=%.2e)", alpha);
    end
    dx = alpha*dx;
end



function jv = jtv(x, v, M, Res0, epsFD, param, perm, invperm)
    vnorm = norm(v);
    if vnorm < 1e-14
        jv = zeros(size(v));  
        return;
    end
    epsFD = sqrt(eps) * max(1, norm(x)) / vnorm;
    Res1 = residual(x+epsFD*v, M, param, perm, invperm);
    jv = (Res1 - Res0) / epsFD;
end

function Res = residual(x, M, param, perm, invperm)
    %x = x(invperm);
    x1 = reshape(x(1:param.N), param.Nx2, param.Nx1);
    x2 = reshape(x(param.N+1:end), param.Nx2, param.Nx1);
    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    %M_old = M;
    %M = GetM(x1, x2, param.problem, param.C);
    M11=M.M11; M22=M.M22; 
    dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;
    M12=M.M12;
    dM12dx1=M.dM12dx1; dM12dx2=M.dM12dx2;

    [dx1ds1, dx2ds2] = DCentral(x1, x2, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, sigma1, sigma2);

    param.dxi = sigma1; param.deta = sigma2;
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
    
    %M = M_old;
    %M = GetM(x1, x2, param.problem, param.C);
    dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;
    dM12dx1=M.dM12dx1; dM12dx2=M.dM12dx2;

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

    %perm(Res);
end