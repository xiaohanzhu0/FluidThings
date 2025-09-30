function [A, b, res] = AssembleLinearSystem(x1, x2, M, param)

    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    M11=M.M11; M22=M.M22; 
    dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;
    

    if param.offdiag == 1
        M12=M.M12;
        dM12dx1=M.dM12dx1; dM12dx2=M.dM12dx2;
    end

    %dM11dx1=0*M.dM11dx1; dM11dx2=0*M.dM11dx2; dM22dx1=0*M.dM22dx1; dM22dx2=0*M.dM22dx2;
    %dM12dx1=0*M.dM12dx1; dM12dx2=0*M.dM12dx2;

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
    if param.offdiag == 1
    A(1:N,N+1:end) = A(1:N,N+1:end) - 2*L.*M12(:);
    A(N+1:end,1:N) = A(N+1:end,1:N) - 2*L.*M12(:);
    end

    if param.repulse == 1
        detJ = dx1ds1.*dx2ds2 - dx1ds2.*dx2ds1;
        JinvT11 = dx2ds2 ./ detJ;  JinvT12 = -dx1ds2 ./ detJ;
        JinvT21 = -dx2ds1 ./ detJ; JinvT22 =  dx1ds1 ./ detJ;
        coef = 1;

        D1 = spdiags([-e/2, e/2], [-Nx2, Nx2], N, N);
        D2 = spdiags([-e/2, e/2], [-1, 1], N, N);
        D11 = D1*( spdiags(coef*JinvT11(:),0,N,N) ) + D2*( spdiags(coef*JinvT12(:),0,N,N) );
        D12 = D1*( spdiags(coef*JinvT21(:),0,N,N) ) + D2*( spdiags(coef*JinvT22(:),0,N,N) );
        D21 = D1*( spdiags(coef*JinvT11(:),0,N,N) ) + D2*( spdiags(coef*JinvT12(:),0,N,N) );
        D22 = D1*( spdiags(coef*JinvT21(:),0,N,N) ) + D2*( spdiags(coef*JinvT22(:),0,N,N) );  % check your indexing
        A(1:N,   1:N)     = A(1:N,1:N)    + D11;
        A(1:N,   N+1:end) = A(1:N,N+1:end) + D12;
        A(N+1:end,1:N)    = A(N+1:end,1:N) + D21;
        A(N+1:end,N+1:end)= A(N+1:end,N+1:end) + D22;
    elseif param.repulse == 2
        coef = 1e-1;
        reg =  0;
        Q1 = M11 .* (dx1ds1.^2) ...
           + 2*M12 .* (dx1ds1 .* dx2ds1) ...
           + M22 .* (dx2ds1.^2);
        
        Q2 = M11 .* (dx1ds2.^2) ...
           + 2*M12 .* (dx1ds2 .* dx2ds2) ...
           + M22 .* (dx2ds2.^2);
        %if use_inverse
          w1 = 1./(sigma1*1^2*Q1.^2 + reg);
          w2 = 1./(sigma2*2^2*Q2.^2 + reg);
        %else             % negative-log
        %  w1 = 1./Q1;
        %  w2 = 1./Q2;
        %end

        W1 = spdiags(w1(:),0,N,N);
        W2 = spdiags(w2(:),0,N,N);
        M11d = spdiags(M11(:) + reg,0,N,N);
        M12d = spdiags(M12(:) + reg,0,N,N);
        M22d = spdiags(M22(:) + reg,0,N,N);
        
        D1 = spdiags([-e/2, e/2], [-Nx2, Nx2], N, N);
        D2 = spdiags([-e/2, e/2], [-1, 1], N, N);
        D11 = 2*( D1*(W1*M11d) + D2*(W2*M11d) );   % ∂R₁/∂x₁
        D12 = 2*( D1*(W1*M12d) + D2*(W2*M12d) );   % ∂R₁/∂x₂
        D21 = 2*( D1*(W1*M12d) + D2*(W2*M12d) );   % ∂R₂/∂x₁  (M₂₁=M₁₂)
        D22 = 2*( D1*(W1*M22d) + D2*(W2*M22d) );   % ∂R₂/∂x₂

        A(1:N,   1:N)     = A(1:N,1:N)    - coef*D11;
        A(1:N,   N+1:end) = A(1:N,N+1:end) - coef*D12;
        A(N+1:end,1:N)    = A(N+1:end,1:N) - coef*D21;
        A(N+1:end,N+1:end)= A(N+1:end,N+1:end) - coef*D22;
    end

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

        % Appending off-diagonal metrix terms
        if param.offdiag == 1
        D(1:N,N+1:end) = D(1:N,N+1:end) - ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) )*4.*dM12dx2(:);
        D(N+1:end,1:N) = D(N+1:end,1:N) - ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) )*4.*dM12dx1(:);
        end

        A = A + D;
    end

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
    
    % Do upwind bias for first derivative
    if param.nonlinear == 2 || (param.nonlinear == 3 && param.forDx == 1)
        D = sparse(N_all, N_all);
        D1 = spdiags([-e/2, e/2], [-Nx2, Nx2], N, N);
        D2 = spdiags([-e/2, e/2], [-1, 1], N, N);

        D(1:N,1:N) = -sigma1*D1.*dx1ds1(:).*dM11dx1(:) - sigma2*D2.*dx1ds2(:).*dM11dx1(:);
        D(N+1:end,N+1:end) = -sigma2*D2.*dx2ds2(:).*dM22dx2(:) - sigma1*D1.*dx2ds1(:).*dM22dx2(:);

        D(1:N,1:N) = D(1:N,1:N) - ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) ).*dM11dx2(:)*2;
        D(1:N,N+1:end) = D(1:N,N+1:end) - ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) ).*dM12dx2(:)*2;
        D(1:N,N+1:end) = D(1:N,N+1:end) + ( sigma1*D1.*dx2ds1(:) + sigma2*D2.*dx2ds2(:) ).*dM22dx1(:);

        D(N+1:end,N+1:end) = D(N+1:end,N+1:end) - ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) ).*dM22dx1(:)*2;
        D(N+1:end,1:N) = D(N+1:end,1:N) - ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) ).*dM12dx1(:)*2;
        D(N+1:end,1:N) = D(N+1:end,1:N) + ( sigma1*D1.*dx1ds1(:) + sigma2*D2.*dx1ds2(:) ).*dM11dx2(:);
        

        if  param.nonlinear == 3 && param.forDx == 1
            D = D * 2;
        end

        A = A + D;
    end
    
    %{
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
    %}

    
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

    % Assemble the interior vector b
    b1 = zeros(Nx2, Nx1);
    b2 = zeros(Nx2, Nx1);
    
    if param.nonlinear == 4
        b1 = b1 + 2*M11.*(ddx1dds1*sigma1^2 + ddx1dds2*sigma2^2);
        b1 = b1 + dM11dx1.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
        b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
        b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);

        b2 = b2 + 2*M22.*(ddx2dds1*sigma1^2 + ddx2dds2*sigma2^2);
        b2 = b2 + dM22dx2.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
        b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
        b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);

        if param.offdiag == 1
        b1 = b1 + 2*M12.*(ddx2dds1*sigma1^2 + ddx2dds2*sigma2^2);
        b1 = b1 + 2*dM12dx2.*(dx2ds1.^2*sigma1^2);
        b1 = b1 + 2*dM12dx2.*(dx2ds2.^2*sigma2^2);

        b2 = b2 + 2*M12.*(ddx1dds1*sigma1^2 + ddx1dds2*sigma2^2);
        b2 = b2 + 2*dM12dx1.*(dx1ds1.^2*sigma1^2);
        b2 = b2 + 2*dM12dx1.*(dx1ds2.^2*sigma2^2);
        end
    else

        if param.nonlinear == 1
            b1 = b1 + dM11dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
            b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
            b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
            if param.offdiag == 1
                b1 = b1 + 2*dM12dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
                %b1 = b1 + 2*dM12dx2 .* (dx2ds1.^2 * sigma1^2 + dx2ds2.^2 * sigma2^2);
                %b1 = b1 + 2*dM12dx1 .* (dx1ds1.*dx2ds1 * sigma1^2 + dx1ds2.*dx2ds2 * sigma2^2);
            end
        end
        %b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
        %b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
    
        if param.nonlinear == 1
            b2 = b2 + dM22dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
            b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
            b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
            if param.offdiag == 1
                b2 = b2 + 2*dM12dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
                %b2 = b2 + 2*dM12dx1 .* (dx1ds1.^2 * sigma1^2 + dx1ds2.^2 * sigma2^2);
                %b2 = b2 + 2*dM12dx2 .* (dx2ds1.*dx1ds1 * sigma1^2 + dx2ds2.*dx1ds2 * sigma2^2);
            end
        end
        %b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
        %b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    end

    if param.smooth ~= 0
        b = [b1(:); b2(:)];
        b = b - param.smooth*A_lap*[x1(:); x2(:)];
        b1 = reshape(b(1:N), Nx2, Nx1);
        b2 = reshape(b(N+1:end), Nx2, Nx1);
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
    res = A*[x1(:); x2(:)] - b;
end


