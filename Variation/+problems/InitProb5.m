function [T1, T2, M11_fun, M12_fun, M22_fun] = InitProb5(cf)

    [data, ~, ~, ~, ~] = readTurtleFields(cf.metric_datapath);
    [x_cv, ~, ~, ~, ~, ~, ~, ~, ~, ~] = readTurtleGrid(cf.airfoil_datapath);
    M_samp.x_metric = [x_cv{1,9}(:,:,1,1), x_cv{1,8}(:,:,1,1);
    x_cv{1,10}(:,:,1,1),x_cv{1,1}(:,:,1,1);
    x_cv{1,11}(:,:,1,1),x_cv{1,2}(:,:,1,1);
    x_cv{1,12}(:,:,1,1),x_cv{1,3}(:,:,1,1);
    x_cv{1,13}(:,:,1,1),x_cv{1,4}(:,:,1,1);
    x_cv{1,14}(:,:,1,1),x_cv{1,5}(:,:,1,1);
    x_cv{1,15}(:,:,1,1),x_cv{1,6}(:,:,1,1);
    x_cv{1,16}(:,:,1,1),x_cv{1,7}(:,:,1,1)];

    M_samp.y_metric = [x_cv{1,9}(:,:,1,2), x_cv{1,8}(:,:,1,2);
    x_cv{1,10}(:,:,1,2),x_cv{1,1}(:,:,1,2);
    x_cv{1,11}(:,:,1,2),x_cv{1,2}(:,:,1,2);
    x_cv{1,12}(:,:,1,2),x_cv{1,3}(:,:,1,2);
    x_cv{1,13}(:,:,1,2),x_cv{1,4}(:,:,1,2);
    x_cv{1,14}(:,:,1,2),x_cv{1,5}(:,:,1,2);
    x_cv{1,15}(:,:,1,2),x_cv{1,6}(:,:,1,2);
    x_cv{1,16}(:,:,1,2),x_cv{1,7}(:,:,1,2)];

    aux1 = cat(1,data{1,1}{1,9}(:,:,1,1:2,1:2), data{1,1}{1,10}(:,:,1,1:2,1:2), ...
        data{1,1}{1,11}(:,:,1,1:2,1:2), data{1,1}{1,12}(:,:,1,1:2,1:2), ...
        data{1,1}{1,13}(:,:,1,1:2,1:2), data{1,1}{1,14}(:,:,1,1:2,1:2), ...
        data{1,1}{1,15}(:,:,1,1:2,1:2), data{1,1}{1,16}(:,:,1,1:2,1:2));
    aux2 = cat(1,data{1,1}{1,8}(:,:,1,1:2,1:2), data{1,1}{1,1}(:,:,1,1:2,1:2), ...
        data{1,1}{1,2}(:,:,1,1:2,1:2), data{1,1}{1,3}(:,:,1,1:2,1:2), ...
        data{1,1}{1,4}(:,:,1,1:2,1:2), data{1,1}{1,5}(:,:,1,1:2,1:2), ...
        data{1,1}{1,6}(:,:,1,1:2,1:2), data{1,1}{1,7}(:,:,1,1:2,1:2));
    aux = cat(2,aux1,aux2);
    dim = size(aux);
    M = reshape(aux, dim(1), dim(2), dim(4), dim(5));
    M_samp.metric(:,:,1) = M(:,:,1,1);
    M_samp.metric(:,:,2) = M(:,:,1,2);
    M_samp.metric(:,:,3) = M(:,:,2,2);
    
    
    %% Gradation to regularize the metric field
    if cf.grade
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
    for iter=1:1
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
    end
    
    M_samp.metric(:,:,1) = M(1,1,:,:);
    M_samp.metric(:,:,2) = M(1,2,:,:);
    M_samp.metric(:,:,3) = M(2,2,:,:);
    end

%%
    M_samp.F11 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
    M_samp.F12 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
    M_samp.F22 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));
    Mfun = @(x1,x2) Prob5Metric(x1, x2, M_samp);
    function M = Prob5Metric(x1,x2,M_samp)
        M.M11 = M_samp.F11(x1,x2);
        M.M12 = M_samp.F12(x1,x2);
        M.M22 = M_samp.F22(x1,x2);
    end
    M11_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
    M12_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
    M22_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));

    gd = Hyperbolic(cf.Nx1, cf.Nx2, cf.alpha, cf.append_trail);
    T1 = gd.x';
    T2 = gd.y';

end

%% Apply TFI
function [X, Y] = transfiniteInterp(b1,b2,b3,b4, bi)
%TRANSFINITEINTERP  Structured grid via transfinite (Coons) interpolation
%  [X,Y] = TRANSFINITEINTERP(b1,b2,b3,b4)
%  [X,Y] = TRANSFINITEINTERP(b1,b2,b3,b4,bi)
%
%  b1: 2×Nx1 bottom (ξ∈[0,1])
%  b2: 2×Nx2 right  (η∈[0,1])
%  b3: 2×Nx1 top    (ξ∈[0,1])
%  b4: 2×Nx2 left   (η∈[0,1])
%  bi: 2×Nx2 interior layer (optional; η∈[0,1])
%
%  If bi is provided, we split at ξ = 0.5 (mid‐column) into:
%    • Left patch:   (b1L, bi, b3L, b4)
%    • Right patch:  (b1R, b2, b3R, bi)
%  and then stitch them, dropping the duplicate seam column.

  narginchk(4,5);
  [~, Nx1] = size(b1);
  [~, Nx2] = size(b2);
  assert(all(size(b3)==[2 Nx1]), 'b3 must be 2×Nx1');
  assert(all(size(b4)==[2 Nx2]), 'b4 must be 2×Nx2');
  if nargin==5
    assert(all(size(bi)==[2 Nx2]), 'bi must be 2×Nx2');
  end

  %— BASE CASE: no interior curve → standard Coons patch
  if nargin<5
    %— parametric coords
    xi  = linspace(0,1,Nx1);
    eta = linspace(0,1,Nx2);
    [XI, ETA] = meshgrid(xi, eta);    % size Nx2×Nx1

    %— unpack
    b1x = b1(1,:);  b1y = b1(2,:);
    b3x = b3(1,:);  b3y = b3(2,:);
    b4x = b4(1,:)'; b4y = b4(2,:)';
    b2x = b2(1,:)'; b2y = b2(2,:)';

    %— corners
    Pbl = b1(:,1);   Pbr = b1(:,end);
    Ptl = b3(:,1);   Ptr = b3(:,end);

    %— blends
    P1x = (1-ETA).*b1x(ones(Nx2,1),:) + ETA.*b3x(ones(Nx2,1),:);
    P1y = (1-ETA).*b1y(ones(Nx2,1),:) + ETA.*b3y(ones(Nx2,1),:);
    P2x = (1-XI).*b4x + XI.*b2x;
    P2y = (1-XI).*b4y + XI.*b2y;

    %— bilinear corner correction
    Cx = (1-XI).*(1-ETA)*Pbl(1) + XI.*(1-ETA)*Pbr(1) ...
       + XI.*ETA*Ptr(1)     + (1-XI).*ETA*Ptl(1);
    Cy = (1-XI).*(1-ETA)*Pbl(2) + XI.*(1-ETA)*Pbr(2) ...
       + XI.*ETA*Ptr(2)     + (1-XI).*ETA*Ptl(2);

    %— assemble
    X = P1x + P2x - Cx;
    Y = P1y + P2y - Cy;

    return
  end

  %— RECURSIVE CASE: interior curve supplied → two patches
  midcol = ceil(Nx1/2);

  % split bottom/top into left/right chunks
  b1L = b1(:,     1:midcol);    b3L = b3(:,     1:midcol);
  b1R = b1(:, midcol: end);     b3R = b3(:, midcol: end);

  % left patch   uses (b1L,  bi, b3L, b4)
  [Xl, Yl] = transfiniteInterp(b1L, bi, b3L, b4);

  % right patch  uses (b1R, b2, b3R, bi)
  [Xr, Yr] = transfiniteInterp(b1R, b2, b3R, bi);

  % stitch seams (drop duplicate interior column at Xr(:,1)==bi)
  X = [ Xl(:,1:midcol),  Xr(:,2:end) ];
  Y = [ Yl(:,1:midcol),  Yr(:,2:end) ];
end


function M_inter = inter(M_xp, M_q)
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