function A_orth = AssembleOrtho(dx1ds1, dx1ds2, dx2ds1, dx2ds2, sigma1, sigma2, id)
    [Nx1, Nx2] = size(dx1ds1);
    N = Nx1*Nx2;
    N_all = N*2;
    A_orth1 = sparse(N_all, N_all);
    A_orth2 = sparse(N_all, N_all);
    A_orth3 = sparse(N_all, N_all);

    D = dx1ds1.*dx1ds2 + dx2ds1.*dx2ds2;
    [dDds1, dDds2] = DCentral(D, D, sigma1, sigma2);
    
    for i = [id.inner, N+id.inner]
        A_orth1(i,i+Nx2+1) = 1/(4*sigma1*sigma2);
        A_orth1(i,i-Nx2+1) = -1/(4*sigma1*sigma2);
        A_orth1(i,i+Nx2-1) = -1/(4*sigma1*sigma2);
        A_orth1(i,i-Nx2-1) = 1/(4*sigma1*sigma2);
    end

    for i = N+[id.l, id.r, -N+id.inner, id.inner]
        A_orth2(i,i+1) = 1/(2*sigma2);
        A_orth2(i,i-1) = -1/(2*sigma2);
    end

    for i = [id.b, id.t, N+id.inner, id.inner]
        A_orth3(i,i+Nx2) = 1/(2*sigma1);
        A_orth3(i,i-Nx2) = -1/(2*sigma1);
    end

    A_orth = -4*A_orth1.*[D(:);D(:)] - 2*A_orth2.*[dDds1(:);dDds1(:)] - 2*A_orth3.*[dDds2(:);dDds2(:)];
end