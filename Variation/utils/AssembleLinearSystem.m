function [A, b, res] = AssembleLinearSystem(x1, x2, M, param)
    [Nx2, Nx1] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    M11=M.M11; M22=M.M22; 
    dM11dx1=M.dM11dx1; dM11dx2=M.dM11dx2; dM22dx1=M.dM22dx1; dM22dx2=M.dM22dx2;

    [dx1ds1, dx2ds2] = DCentral(x1, x2, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, sigma1, sigma2);
    id = GetIndex(Nx1, Nx2);

    A = sparse(N_all, N_all);
    
    if param.method == 2 % For approximate cost function
        % coef(k,alpha,i)
        coef111 = -2*sigma1^2*M11.*dx1ds1.*M11.*dx1ds1; coef111 = coef111(:);
        coef112 = -2*sigma1^2*M11.*dx1ds1.*M22.*dx2ds1; coef112 = coef112(:);
        coef121 = -2*sigma2^2*M11.*dx1ds2.*M11.*dx1ds2; coef121 = coef121(:);
        coef122 = -2*sigma2^2*M11.*dx1ds2.*M22.*dx2ds2; coef122 = coef122(:);

        coef211 = -2*sigma1^2*M22.*dx2ds1.*M11.*dx1ds1; coef211 = coef211(:);
        coef212 = -2*sigma1^2*M22.*dx2ds1.*M22.*dx2ds1; coef212 = coef212(:);
        coef221 = -2*sigma2^2*M22.*dx2ds2.*M11.*dx1ds2; coef221 = coef221(:);
        coef222 = -2*sigma2^2*M22.*dx2ds2.*M22.*dx2ds2; coef222 = coef222(:);

        % first component of s
        for i = [id.b, id.inner, id.t]
            % Same way to handle ds1
            A(i,i) = A(i,i) - 2*coef111(i);  % (i,i) for x1 contribution on x1
            A(i,i-Nx2) = A(i,i-Nx2) + coef111(i);
            A(i,i+Nx2) = A(i,i+Nx2) + coef111(i);

            A(i,i+N) = A(i,i+N) - 2*coef112(i); % (i,i+N) for x2 contribution on x1
            A(i,i+N-Nx2) = A(i,i+N-Nx2) + coef112(i);
            A(i,i+N+Nx2) = A(i,i+N+Nx2) + coef112(i);

            % Different way to handle ds2
            if ismember(i,id.inner) % Central difference at interior
                A(i,i) = A(i,i) - 2*coef121(i);
                A(i,i-1) = A(i,i-1) + coef121(i);
                A(i,i+1) = A(i,i+1) + coef121(i);
                A(i,i+N) = A(i,i+N) - 2*coef122(i);
                A(i,i+N-1) = A(i,i+N-1) + coef122(i);
                A(i,i+N+1) = A(i,i+N+1) + coef122(i);
            elseif ismember(i,id.b) % One sided difference at boundary
                A(i,i) = A(i,i) - 2*coef121(i);
                A(i,i+1) = A(i,i+1) + 2*coef121(i);

                A(i,i+N) = A(i,i+N) + 2*coef122(i);
                A(i,i+N+1) = A(i,i+N+1) - 5*coef122(i);
                A(i,i+N+2) = A(i,i+N+2) + 4*coef122(i);
                A(i,i+N+3) = A(i,i+N+3) - 1*coef122(i);
            elseif ismember(i,id.t)
                A(i,i) = A(i,i) - 2*coef121(i);
                A(i,i-1) = A(i,i-1) + 2*coef121(i);

                A(i,i+N) = A(i,i+N) + 2*coef122(i);
                A(i,i+N-1) = A(i,i+N-1) - 5*coef122(i);
                A(i,i+N-2) = A(i,i+N-2) + 4*coef122(i);
                A(i,i+N-3) = A(i,i+N-3) - 1*coef122(i);
            end
        end


        % second component of interior
        for i = N+[id.l, id.inner, id.r]
            A(i,i) = A(i,i) - 2*coef222(i-N); % (i,i-N) for x1 contribution on x2
            A(i,i-1) = A(i,i-1) + coef222(i-N);
            A(i,i+1) = A(i,i+1) + coef222(i-N);

            A(i,i-N) = A(i,i-N) - 2*coef221(i-N);  % (i,i) for x2 contribution on x2
            A(i,i-N-1) = A(i,i-N-1) + coef221(i-N);
            A(i,i-N+1) = A(i,i-N+1) + coef221(i-N);

            if ismember(i-N,id.inner)
                A(i,i-N) = A(i,i-N) - 2*coef211(i-N);
                A(i,i-N-Nx2) = A(i,i-N-Nx2) + coef211(i-N);
                A(i,i-N+Nx2) = A(i,i-N+Nx2) + coef211(i-N);
    
                A(i,i) = A(i,i) - 2*coef212(i-N);
                A(i,i-Nx2) = A(i,i-Nx2) + coef212(i-N);
                A(i,i+Nx2) = A(i,i+Nx2) + coef212(i-N);
            elseif ismember(i-N,id.l) % One sided difference for dx1/ds1
                A(i,i) = A(i,i) - 2*coef212(i-N);
                A(i,i+Nx2) = A(i,i+Nx2) + 2*coef212(i-N);

                A(i,i-N) = A(i,i-N) + 2*coef211(i-N);
                A(i,i-N+Nx2) = A(i,i-N+Nx2) - 5*coef211(i-N);
                A(i,i-N+2*Nx2) = A(i,i-N+2*Nx2) + 4*coef211(i-N);
                A(i,i-N+3*Nx2) = A(i,i-N+3*Nx2) - 1*coef211(i-N);
            elseif ismember(i-N,id.r)
                A(i,i) = A(i,i) - 2*coef212(i-N);
                A(i,i-Nx2) = A(i,i-Nx2) + 2*coef212(i-N);

                A(i,i-N) = A(i,i-N) + 2*coef211(i-N);
                A(i,i-N-Nx2) = A(i,i-N-Nx2) - 5*coef211(i-N);
                A(i,i-N-2*Nx2) = A(i,i-N-2*Nx2) + 4*coef211(i-N);
                A(i,i-N-3*Nx2) = A(i,i-N-3*Nx2) - 1*coef211(i-N);
            end
        end
    
        %A = A + 1*A_orth;
        
        % For Dirichlet conditions
        for i =   [id.l, id.r, id.corner]; A(i,:) = 0; A(i,i) = 1; end
        for i = N+[id.b, id.t, id.corner]; A(i,:) = 0; A(i,i) = 1; end
    
        b_aux1 = dM11dx1.*dx1ds1+dM11dx2.*dx2ds1;
        b_aux2 = dM22dx1.*dx1ds1+dM22dx2.*dx2ds1;
        b_aux3 = dM11dx1.*dx1ds2+dM11dx2.*dx2ds2;
        b_aux4 = dM22dx1.*dx1ds2+dM22dx2.*dx2ds2;
        b1 = sigma1^4 * M11.*dx1ds1.*(b_aux1.*dx1ds1.^2 + b_aux2.*dx2ds1.^2) + ...
             sigma2^4 * M11.*dx1ds2.*(b_aux3.*dx1ds2.^2 + b_aux4.*dx2ds2.^2);
        b2 = sigma1^4 * M22.*dx2ds1.*(b_aux1.*dx1ds1.^2 + b_aux2.*dx2ds1.^2) + ...
             sigma2^4 * M22.*dx2ds2.*(b_aux3.*dx1ds2.^2 + b_aux4.*dx2ds2.^2);

    elseif param.method == 1 % For alternative cost function
        
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
        
        % Do upwind bias for first derivative
        if param.nonlinear == 2 || (param.nonlinear == 3 && param.forDx == 1)
            D = sparse(N_all, N_all);
            %D1 = spdiags([-e, e], [-Nx2, 0], N, N);

            %D2 = spdiags([-e, e], [-1, 0], N, N);
            D1 = spdiags([-e/2, e/2], [-Nx2, Nx2], N, N);
            D2 = spdiags([-e/2, e/2], [-1, 1], N, N);
    
            D(1:N,1:N) = -sigma1*D1.*dx1ds1(:).*dM11dx1(:) - sigma2*D2.*dx1ds2(:).*dM11dx1(:);
            D(N+1:end,N+1:end) = -sigma2*D2.*dx2ds2(:).*dM22dx2(:) - sigma1*D1.*dx2ds1(:).*dM22dx2(:);

            if  param.nonlinear == 3 && param.forDx == 1
                D = D * 2;
            end

            A = A + D;
        end

        %A(id.inner,:) = -2*A(id.inner,:).*Mii(id.inner,:);
        %A(N+id.inner,:) = -2*A(N+id.inner,:).*Mii(N+id.inner,:);
        
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
        if param.nonlinear == 1 || (param.nonlinear == 3 && param.forDx == 0)
            b1 = b1 + dM11dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
        end
        b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
        b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
    
        b2 = 0;
        if param.nonlinear == 1 || (param.nonlinear == 3 && param.forDx == 0)
            b2 = b2 + dM22dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
        end
        b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
        b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    end

    % RHS changes due to Dirichlet on boundary
    b1(:,1) = x1(:,1).*n_left(1,:)' + x2(:,1).*n_left(2,:)';
    b1(:,end) = x1(:,end).*n_right(1,:)' + x2(:,end).*n_right(2,:)';
    b1(1,:) = x1(1,:).*n_bottom(1,:) + x2(1,:).*n_bottom(2,:);
    b1(end,:) = x1(end,:).*n_top(1,:) + x2(end,:).*n_top(2,:);

    % RHS changes due to Neumann on boundary
    b2(:,1) = 0; b2(:,end) = 0; b2(1,:) = 0; b2(end,:) = 0;

    b = [b1(:); b2(:)];

    res = norm(A*[x1(:);x2(:)] - b);
end

