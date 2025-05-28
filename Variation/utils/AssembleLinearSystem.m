function [A, b] = AssembleLinearSystem(x1, x2, M, param)
    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    M11=M.M11; M22=M.M22; 
    dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;

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
    A = blkdiag(L, L);

    Mii = [M11(:); M22(:)];
    A = -2*A.*Mii;

    if param.nonlinear == 4
        D = sparse(N_all, N_all);
        D1 = spdiags([-e/2, e/2], [-Nx2, Nx2], N, N);
        D2 = spdiags([-e/2, e/2], [-1, 1], N, N);

        D(1:N,1:N) = D(1:N,1:N) - ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) )*2.*dM11dx1(:);
        D(1:N,1:N) = D(1:N,1:N) - ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) )*2.*dM11dx2(:);
        D(1:N,N+1:end) = D(1:N,N+1:end) - ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) )*2.*dM11dx2(:);
        D(1:N,N+1:end) = D(1:N,N+1:end) + ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) )*2.*dM22dx1(:);

        D(N+1:end,N+1:end) = D(N+1:end,N+1:end) - ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) )*2.*dM22dx2(:);
        D(N+1:end,N+1:end) = D(N+1:end,N+1:end) - ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) )*2.*dM22dx1(:);
        D(N+1:end,1:N) = D(N+1:end,1:N) - ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) )*2.*dM22dx1(:);
        D(N+1:end,1:N) = D(N+1:end,1:N) + ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) )*2.*dM11dx2(:);

        A = A + D;
    end
    
    % Do upwind bias for first derivative
    if param.nonlinear == 2 || (param.nonlinear == 3 && param.forDx == 1)
        D = sparse(N_all, N_all);
        D1 = spdiags([-e/2, e/2], [-Nx2, Nx2], N, N);
        D2 = spdiags([-e/2, e/2], [-1, 1], N, N);

        D(1:N,1:N) = -sigma1*D1.*dx1ds1(:).*dM11dx1(:) - sigma2*D2.*dx1ds2(:).*dM11dx1(:);
        D(N+1:end,N+1:end) = -sigma2*D2.*dx2ds2(:).*dM22dx2(:) - sigma1*D1.*dx2ds1(:).*dM22dx2(:);

        if  param.nonlinear == 3 && param.forDx == 1
            D = D * 2;
        end

        A = A + D;
    end
    
    A(id.l,:) = 0;
    A(sub2ind(size(A), id.l, id.l))   = n_left(1,2:end-1);
    A(sub2ind(size(A), id.l, N+id.l)) = n_left(2,2:end-1);
    A(id.l+N,:) = 0;
    A(sub2ind(size(A), N+id.l, id.l)) = t_left(1,2:end-1);
    A(sub2ind(size(A), N+id.l, N+id.l)) = t_left(2,2:end-1);
    A(sub2ind(size(A), N+id.l, 1*Nx2+id.l)) = -4/3*t_left(1,2:end-1);
    A(sub2ind(size(A), N+id.l, N+1*Nx2+id.l)) = -4/3*t_left(2,2:end-1);
    A(sub2ind(size(A), N+id.l, 2*Nx2+id.l)) = 1/3*t_left(1,2:end-1);
    A(sub2ind(size(A), N+id.l, N+2*Nx2+id.l)) = 1/3*t_left(2,2:end-1);
    
    A(id.r,:) = 0;
    A(sub2ind(size(A), id.r, id.r))   = n_right(1,2:end-1);
    A(sub2ind(size(A), id.r, N+id.r)) = n_right(2,2:end-1);
    A(id.r+N,:) = 0;
    A(sub2ind(size(A), N+id.r, id.r)) = t_right(1,2:end-1);
    A(sub2ind(size(A), N+id.r, N+id.r)) = t_right(2,2:end-1);
    A(sub2ind(size(A), N+id.r, -1*Nx2+id.r)) = -4/3*t_right(1,2:end-1);
    A(sub2ind(size(A), N+id.r, N-1*Nx2+id.r)) = -4/3*t_right(2,2:end-1);
    A(sub2ind(size(A), N+id.r, -2*Nx2+id.r)) = 1/3*t_right(1,2:end-1);
    A(sub2ind(size(A), N+id.r, N-2*Nx2+id.r)) = 1/3*t_right(2,2:end-1);

    A(id.b,:) = 0;
    A(sub2ind(size(A), id.b, id.b))   = n_bottom(1,2:end-1);
    A(sub2ind(size(A), id.b, N+id.b)) = n_bottom(2,2:end-1);
    A(id.b+N,:) = 0;
    A(sub2ind(size(A), N+id.b, id.b)) = t_bottom(1,2:end-1);
    A(sub2ind(size(A), N+id.b, N+id.b)) = t_bottom(2,2:end-1);
    A(sub2ind(size(A), N+id.b, 1+id.b)) = -4/3*t_bottom(1,2:end-1);
    A(sub2ind(size(A), N+id.b, N+1+id.b)) = -4/3*t_bottom(2,2:end-1);
    A(sub2ind(size(A), N+id.b, 2+id.b)) = 1/3*t_bottom(1,2:end-1);
    A(sub2ind(size(A), N+id.b, N+2+id.b)) = 1/3*t_bottom(2,2:end-1);

    A(id.t,:) = 0;
    A(sub2ind(size(A), id.t, id.t))   = n_top(1,2:end-1);
    A(sub2ind(size(A), id.t, N+id.t)) = n_top(2,2:end-1);
    A(id.t+N,:) = 0;
    A(sub2ind(size(A), N+id.t, id.t)) = t_top(1,2:end-1);
    A(sub2ind(size(A), N+id.t, N+id.t)) = t_top(2,2:end-1);
    A(sub2ind(size(A), N+id.t, -1+id.t)) = -4/3*t_top(1,2:end-1);
    A(sub2ind(size(A), N+id.t, N-1+id.t)) = -4/3*t_top(2,2:end-1);
    A(sub2ind(size(A), N+id.t, -2+id.t)) = 1/3*t_top(1,2:end-1);
    A(sub2ind(size(A), N+id.t, N-2+id.t)) = 1/3*t_top(2,2:end-1);
    
    for i = id.corner; A(i, :) = 0; A(i, i) = 1; end
    for i = N+id.corner; A(i, :) = 0; A(i, i) = 1; end

    % Assemble the interior vector b
    b1 = 0;
    b2 = 0;
    

    if param.nonlinear == 4
        b1 = b1 + 2*M11.*(ddx1dds1*sigma1^2 + ddx1dds2*sigma2^2);
        b1 = b1 + dM11dx1.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
        b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
        b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);

        b2 = b2 + 2*M22.*(ddx2dds1*sigma1^2 + ddx2dds2*sigma2^2);
        b2 = b2 + dM22dx2.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
        b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
        b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    else
        if param.nonlinear == 1 || (param.nonlinear == 3 && param.forDx == 0)
            b1 = b1 + dM11dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
        end
        b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
        b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
    
        if param.nonlinear == 1 || (param.nonlinear == 3 && param.forDx == 0)
            b2 = b2 + dM22dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
        end
        b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
        b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    end

    % RHS changes due to Dirichlet on boundary
    if param.nonlinear == 4
        b1(:,1) = 0;
        b1(:,end) = 0;
        b1(1,:) = 0;
        b1(end,:) = 0;
         % RHS changes due to Neumann on boundary
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
    else
        b1(:,1) = x1(:,1).*n_left(1,:)' + x2(:,1).*n_left(2,:)';
        b1(:,end) = x1(:,end).*n_right(1,:)' + x2(:,end).*n_right(2,:)';
        b1(1,:) = x1(1,:).*n_bottom(1,:) + x2(1,:).*n_bottom(2,:);
        b1(end,:) = x1(end,:).*n_top(1,:) + x2(end,:).*n_top(2,:);

        b2(:,1) = 0; b2(:,end) = 0; b2(1,:) = 0; b2(end,:) = 0;
    end

    b = [b1(:); b2(:)];
end

