function [A, b, res] = AssembleLinearSystemApprox(x1, x2, Mfun, param)
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


    dxds1 = cat(3, dx1ds1, dx2ds1);
    dxds2 = cat(3, dx1ds2, dx2ds2);
    dMi1dx1 = cat(3, dM11dx1, dM12dx1); dMi2dx1 = cat(3, dM12dx1, dM22dx1);
    dMdx1 = cat(4, dMi1dx1, dMi2dx1);
    dMi1dx2 = cat(3, dM11dx2, dM12dx2); dMi2dx2 = cat(3, dM12dx2, dM22dx2);
    dMdx2 = cat(4, dMi1dx2, dMi2dx2);
    dMdx = cat(5, dMdx1, dMdx2);

    Mi1 = cat(3, M11, M12); Mi2 = cat(3, M12, M22); 
    M = cat(4, Mi1, Mi2);

    A111 = -8*sigma1^2 * dot(squeeze(M(:,:,1,:)),dxds1,3) .* dot(squeeze(M(:,:,:,1)),dxds1,3);
    A112 = -8*sigma2^2 * dot(squeeze(M(:,:,1,:)),dxds2,3) .* dot(squeeze(M(:,:,:,1)),dxds2,3);
    A121 = -8*sigma1^2 * dot(squeeze(M(:,:,1,:)),dxds1,3) .* dot(squeeze(M(:,:,:,2)),dxds1,3);
    A122 = -8*sigma2^2 * dot(squeeze(M(:,:,1,:)),dxds2,3) .* dot(squeeze(M(:,:,:,2)),dxds2,3);
    A211 = -8*sigma1^2 * dot(squeeze(M(:,:,2,:)),dxds1,3) .* dot(squeeze(M(:,:,:,1)),dxds1,3);
    A212 = -8*sigma2^2 * dot(squeeze(M(:,:,2,:)),dxds2,3) .* dot(squeeze(M(:,:,:,1)),dxds2,3);
    A221 = -8*sigma1^2 * dot(squeeze(M(:,:,2,:)),dxds1,3) .* dot(squeeze(M(:,:,:,2)),dxds1,3);
    A222 = -8*sigma2^2 * dot(squeeze(M(:,:,2,:)),dxds2,3) .* dot(squeeze(M(:,:,:,2)),dxds2,3);

    e = ones(N,1);
    L1 = spdiags([e, -2*e, e], [-Nx2, 0, Nx2], N, N);
    L2 = spdiags([e, -2*e, e], [-1, 0, 1], N, N);

    A = sparse(N_all, N_all);
    A(1:N,1:N) = A(1:N,1:N) + L1.*A111(:) + L2.*A112(:);
    A(N+1:end,N+1:end) = A(N+1:end,N+1:end) + L1.*A221(:) + L2.*A222(:);

    A(1:N,N+1:end) = A(1:N,N+1:end) + L1.*A121(:) + L2.*A122(:);
    A(N+1:end,1:N) = A(N+1:end,1:N) + L1.*A211(:) + L2.*A212(:);

    clear A111 A112 A121 A122 A211 A212 A221 A222 L1 L2 e Mi1 Mi2 dM11dx1 dM11dx2 dM12dx1 dM12dx2 dM22dx1 dM22dx2
    clear dMdx1 dMdx2 dMi1dx1 dMi1dx2 dMi2dx1 dMi2dx2 

    id = GetIndex(Nx1, Nx2);
    [A_orth, A_lap] = AssembleOrtho(dx1ds1, dx1ds2, dx2ds1, dx2ds2, sigma1, sigma2, id);
    %A = A + 50*A_orth + 1e2*A_lap;

    clear A_orth A_lap
 
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
    clear B I J V

    for i = id.corner; A(i, :) = 0; A(i, i) = 1; end
    for i = N+id.corner; A(i, :) = 0; A(i, i) = 1; end

   aux1 = sum(sum(sum(dMdx .* reshape(dxds1,[Nx2,Nx1,2,1,1]) .*...
                              reshape(dxds1,[Nx2,Nx1,1,2,1]) .*...
                              reshape(dxds1,[Nx2,Nx1,1,1,2]), 5), 4), 3);
   aux2 = sum(sum(sum(dMdx .* reshape(dxds2,[Nx2,Nx1,2,1,1]) .*...
                              reshape(dxds2,[Nx2,Nx1,1,2,1]) .*...
                              reshape(dxds2,[Nx2,Nx1,1,1,2]), 5), 4), 3);

   b1 = 4*sigma1^4 * dot(squeeze(M(:,:,1,:)), dxds1, 3) .* aux1 + ...
          sigma2^4 * dot(squeeze(M(:,:,1,:)), dxds2, 3) .* aux2;

   b2 = 4*sigma1^4 * dot(squeeze(M(:,:,2,:)), dxds1, 3) .* aux1 + ...
          sigma2^4 * dot(squeeze(M(:,:,2,:)), dxds2, 3) .* aux2;


    b1(:,1) = x1(:,1).*n_left(1,:)' + x2(:,1).*n_left(2,:)';
    b1(:,end) = x1(:,end).*n_right(1,:)' + x2(:,end).*n_right(2,:)';
    b1(1,:) = x1(1,:).*n_bottom(1,:) + x2(1,:).*n_bottom(2,:);
    b1(end,:) = x1(end,:).*n_top(1,:) + x2(end,:).*n_top(2,:);
    b2(:,1) = 0; b2(:,end) = 0; b2(1,:) = 0; b2(end,:) = 0;


   b = [b1(:); b2(:)];


   res = norm(A*[x1(:);x2(:)] - b);
end


function [A_orth, A_lap] = AssembleOrtho(dx1ds1, dx1ds2, dx2ds1, dx2ds2, s1, s2, id)
    [Nx1, Nx2] = size(dx1ds1);
    N = Nx1*Nx2;
    N_all = N*2;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;
    cellArea = sigma1*sigma2;

    g12 = dx1ds1.*dx1ds2 + dx2ds1.*dx2ds2;
    g11 = dx1ds1.*dx1ds1 + dx2ds1.*dx2ds1;
    g22 = dx1ds2.*dx1ds2 + dx2ds2.*dx2ds2;
    gtilde = g12 ./ sqrt(g11 .* g22);
    e = ones(N,1);

    D1 = spdiags([-e, e], [0, Nx2], N, N) / sigma1;
    D2 = spdiags([-e, e], [0, 1   ], N, N) / sigma2;

    G = spdiags( cellArea * gtilde(:), 0, N, N );
    A_orth = - ( D1' * G * D2  +  D2' * G * D1 );
    A_lap = D1' * D1  +  D2' * D2;
    A_orth = blkdiag(A_orth, A_orth);
    A_lap = blkdiag(A_lap, A_lap);
end

