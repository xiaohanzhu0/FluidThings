function [T1, T2, M11_fun,M12_fun,M22_fun] = InitProb9(cf)
        [data, ~, ~, ~, ~] = readTurtleFields(cf.metric_datapath);
        [x_cv, ~, ~, ~, ~, ~, ~, ~, ~, ~] = readTurtleGrid(cf.airfoil_datapath);

        M_samp.x_metric = x_cv{1,cf.block_idx}(:,:,1,1);
        M_samp.y_metric = x_cv{1,cf.block_idx}(:,:,1,2);

        M = data{1,1}{1,cf.block_idx}(:,:,1,1:2,1:2);
        dim = size(M);
        M = reshape(M, dim(2), dim(1), dim(4), dim(5));

        M_samp.metric(:,:,1) = M(:,:,1,1);
        M_samp.metric(:,:,2) = M(:,:,1,2);
        M_samp.metric(:,:,3) = M(:,:,2,2);

        M_samp.x_metric = M_samp.x_metric';
        M_samp.y_metric = M_samp.y_metric';
    
    %%

    beta = 1.;
    nbr = [-1 0;   % left
            1 0;   % right
            0 -1;  % down
            0  1];

    x1 = M_samp.x_metric;
    x2 = M_samp.y_metric;
    [Nx1, Nx2] = size(x1);
    x = cat(1,reshape(x1,1,Nx1,Nx2), reshape(x2,1,Nx1,Nx2));
    M = zeros(2,2,Nx1,Nx2);
    M(1,1,:,:) = M_samp.metric(:,:,1);
    M(1,2,:,:) = M_samp.metric(:,:,2);
    M(2,1,:,:) = M_samp.metric(:,:,2);
    M(2,2,:,:) = M_samp.metric(:,:,3);
    M_in = M;
    M_out = zeros(2,2,Nx1,Nx2);
    for iter=1:0
    for i=1:Nx1
        for j=1:Nx2
            p = x(:,i,j);
            M_p = M(:,:,i,j);
    
            for k = 1:4
                ii = i + nbr(k,1);
                jj = j + nbr(k,2);
        
                if ii < 1 || ii > Nx1 || jj < 1 || jj > Nx2
                    continue
                end
    
                q = x(:,ii,jj);
                M_q = M(:,:,ii,jj);
                pq = q-p;
    
                eta = (1 + sqrt(pq'*M_p*pq)*log(beta))^(-2);
                M_xp = eta*M_p;
    
                M_inter = inter(M_xp,M_q);
                M_q = M_inter;
                M(:,:,ii,jj) = M(:,:,ii,jj) + 0.5*(real(M_q)-M(:,:,ii,jj));
            end
        end
    end

    
    M_samp.metric(:,:,1) = M(1,1,:,:);
    M_samp.metric(:,:,2) = M(1,2,:,:);
    M_samp.metric(:,:,3) = M(2,2,:,:);
    end

    M11_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
    M12_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
    M22_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));

    T1 = M_samp.x_metric;
    T2 = M_samp.y_metric;

end

function M_inter = inter(M_xp, M_q)
    %N = M_q \ M_xp;
    N = M_xp \ M_q;
    [V,D] = eig(N);
    V(:,1) = V(:,1) / norm(V(:,1));
    V(:,2) = V(:,2) / norm(V(:,2));
    lam_xp = V'*M_xp*V;
    lam_q = V'*M_q*V;
    lam1_inter = max(lam_xp(1,1), lam_q(1,1));
    lam2_inter = max(lam_xp(2,2), lam_q(2,2));

    M_inter = inv(V')*diag([lam1_inter, lam2_inter])*inv(V);
end