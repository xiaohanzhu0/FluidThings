function gd = Hyperbolic(JMAX, KMAX, alpha, append_trail)

params.JMAX = JMAX;
params.KMAX = KMAX;
params.dxi  = 1 / (params.JMAX-1);
params.deta = 1 / (params.KMAX-1);
params.append_trail = append_trail;
params.trail_factor = 1.005;
params.alpha = alpha;

addpath('~/Files/data/Mesh_generation/Airfoil/foil1');
foil = readmatrix('A-airfoil.txt', 'NumHeaderLines', 1);
foil = foil(:,1:2);
N = ceil(params.JMAX/2);
j = 1:N;

%x = 0.5 - 0.5*cos(pi*(N-j)/(N-1));
x = 1 - cos(pi*(N-j)/(N-1)/2);
[~,i_head] = min(foil(:,1));
foil_up = foil(i_head:end,:);
foil_down = foil(1:i_head,:);
y_up = interp1(foil_up(:,1),foil_up(:,2),x,'spline');
y_down = interp1(foil_down(:,1),foil_down(:,2),x,'spline');

Nup = N+0;
j = 1:Nup;
xup = 1 - cos(pi*(Nup-j)/(Nup-1)/2);
[~,i_head] = min(foil(:,1));
foil_up = foil(i_head:end,:);
y_up = interp1(foil_up(:,1),foil_up(:,2),xup,'spline');
Ndown = N-0;
j = 1:Ndown;
xdown = 1 - cos(pi*(Ndown-j)/(Ndown-1)/2);
foil_down = foil(1:i_head,:);
y_down = interp1(foil_down(:,1),foil_down(:,2),xdown,'spline');


if params.append_trail
    counter = 0;
    while x(1) <= 2
        x = [x(1)+(x(1)-x(2))*params.trail_factor, x];
        y_up = [y_up(1), y_up];
        y_down = [y_down(1), y_down];
        counter = counter + 1;
    end
    params.JMAX = params.JMAX + counter*2;
end

gd.x = zeros(params.JMAX, params.KMAX);
gd.y = zeros(params.JMAX, params.KMAX);
if params.append_trail
    x = [x(1:end-1), flip(x)];
else
    x = [xup(1:end-1), flip(xdown)];
end
y = [y_up(1:end-1), flip(y_down)];
gd.x(:,1) = x;
gd.y(:,1) = y;

gd = solve_hyperbolic(gd, params);


figure
for i=1:size(gd.x,1)
    plot(gd.x(i,:),gd.y(i,:)); hold on
end
for i=1:size(gd.x,2)
    plot(gd.x(:,i),gd.y(:,i))
end
end

%% FUNCTION DEFINITIONS
function gd = solve_hyperbolic(gd, params)
    dxdxi = [gd.x(2,1)-gd.x(end,1); gd.x(3:end,1)-gd.x(1:end-2,1); gd.x(1,1)-gd.x(end-1,1)]/(2*params.dxi);
    dydxi = [gd.y(2,1)-gd.y(end,1); gd.y(3:end,1)-gd.y(1:end-2,1); gd.y(1,1)-gd.y(end-1,1)]/(2*params.dxi);
    [nx, ny] = compute_wall_normals(gd.x(:,1), gd.y(:,1));
    dn = [nx, ny] * params.deta;
    x0 = gd.x(:,1) - dn(:,1);
    y0 = gd.y(:,1) - dn(:,2);
    dxdeta = (gd.x(:,1) - x0)/params.deta;
    dydeta = (gd.y(:,1) - y0)/params.deta;

    for k = 2:params.KMAX
        a = (dxdxi.*dxdeta - dydxi.*dydeta) ./ (dxdxi.^2 + dydxi.^2);
        b = (dxdxi.*dydeta + dydxi.*dxdeta) ./ (dxdxi.^2 + dydxi.^2);
        lambda = sqrt(a.^2 + b.^2);

        params.deta = params.alpha^k*params.dxi/max(lambda);

        dxdxi = [gd.x(2,k-1)-gd.x(end,k-1); gd.x(3:end,k-1)-gd.x(1:end-2,k-1); gd.x(1,k-1)-gd.x(end-1,k-1)]/(2*params.dxi);
        dydxi = [gd.y(2,k-1)-gd.y(end,k-1); gd.y(3:end,k-1)-gd.y(1:end-2,k-1); gd.y(1,k-1)-gd.y(end-1,k-1)]/(2*params.dxi);
        if k == 2
            dxdxi = [gd.x(2,1)-gd.x(end,1); gd.x(3:end,1)-gd.x(1:end-2,1); gd.x(1,1)-gd.x(end-1,1)]/(2*params.dxi);
            dydxi = [gd.y(2,1)-gd.y(end,1); gd.y(3:end,1)-gd.y(1:end-2,1); gd.y(1,1)-gd.y(end-1,1)]/(2*params.dxi);
            [nx, ny] = compute_wall_normals(gd.x(:,1), gd.y(:,1));
            dn = [nx, ny] * params.deta;
            x0 = gd.x(:,1) - dn(:,1);
            y0 = gd.y(:,1) - dn(:,2);
            dxdeta = (gd.x(:,1) - x0)/params.deta;
            dydeta = (gd.y(:,1) - y0)/params.deta;
        else
            [nx, ny] = compute_wall_normals(gd.x(:,k-1), gd.y(:,k-1));
            dn = [nx, ny] * params.deta;
            x0 = gd.x(:,k-1) - dn(:,1);
            y0 = gd.y(:,k-1) - dn(:,2);
            dxdeta = (gd.x(:,k-1) - x0)/params.deta;
            dydeta = (gd.y(:,k-1) - y0)/params.deta;
        end
        
        Jinv = dxdxi.*dydeta - dxdeta.*dydxi;
        A(1,1,:) = dxdeta; A(1,2,:) = dydeta; A(2,1,:) = dydeta; A(2,2,:) = -dxdeta;
        B(1,1,:) = dxdxi; B(1,2,:) = dydxi; B(2,1,:) = -dydxi; B(2,2,:) = dxdxi;
        f(2,:) = 2*Jinv; f(1,:) = 0;

        M = block_tridiag(-A/(2*params.dxi), B/params.deta, A/(2*params.dxi));
        Bmat = block_tridiag(0*A, B/params.deta, 0*A);
        w_old = [gd.x(:,k-1), gd.y(:,k-1)]';

        RHS = f(:) + Bmat*w_old(:);

        % Dirichlet boundary condition
        RHS(1) = gd.x(1,1); M(1,:) = 0; M(1,1) = 1;
        RHS(end-1) = gd.x(end,1); M(end-1,:) = 0; M(end-1,end-1) = 1;
        RHS(2) = 0; M(2,:) = 0; M(2,2) = 1; M(2,4) = -1;
        RHS(end) = 0; M(end,:) = 0; M(end,end) = 1; M(end,end-2) = -1;

        w = M \ RHS;
        w = reshape(w, 2, params.JMAX)';

        %if k > 5
        %    w = [w(1,:); w(1:end-2,:)/3 + w(2:end-1,:)/3 + w(3:end,:)/3 ;w(end,:)];
        %end
        
        gd.x(:,k) = w(:,1);
        gd.y(:,k) = w(:,2);

    end
end

function [nx, ny] = compute_wall_normals(x, y)
  Ni = numel(x);
  %--- build index arrays for periodic central difference
  ip = [2:Ni,    1   ];   % i+1 with wrap
  im = [  Ni, 1:Ni-1 ];   % i-1 with wrap

  %--- compute tangent vectors T = (dx/dξ, dy/dξ)
  dx = (x(ip) - x(im)) * (1/2);   % assume uniform Δξ=1
  dy = (y(ip) - y(im)) * (1/2);

  %--- rotate +90°:  N = (Dy, -Dx)
  nxs =  dy;
  nys = -dx;

  %--- normalize
  L = sqrt(nxs.^2 + nys.^2);
  nx = nxs ./ L;
  ny = nys ./ L;
end

function M = block_tridiag(L, D, U)
    n = size(D, 3);
    N2 = 2*n;

    % local offsets within each 2×2 block
    local_i = [1; 2; 1; 2];
    local_j = [1; 1; 2; 2];

    %% Diagonal blocks D
    % replicate block indices 1..n in a 4×n array
    blocksD = repmat(1:n, 4, 1);
    % global row, col for each of the 4 entries per block
    rowD = 2*(blocksD - 1) + local_i;
    colD = 2*(blocksD - 1) + local_j;
    valD = reshape(D, 4, n);

    %% Sub-diagonal blocks L at block-rows 2..n
    blocksL = repmat(2:n, 4, 1);
    rowL = 2*(blocksL - 1) + local_i;
    colL = 2*(blocksL - 2) + local_j;
    valL = reshape(L, 4, n);
    % only keep columns 2..n
    rowL = rowL(:);
    colL = colL(:);
    valL = valL(:,2:end);
    valL = valL(:);

    %% Super-diagonal blocks U at block-rows 1..n–1
    blocksU = repmat(1:n-1, 4, 1);
    rowU = 2*(blocksU - 1) + local_i;
    colU = 2*(blocksU    ) + local_j;
    valU = reshape(U, 4, n);
    % only keep columns 1..n-1
    rowU = rowU(:);
    colU = colU(:);
    valU = valU(:,1:end-1);
    valU = valU(:);

    %% Combine all entries
    rowD = rowD(:);
    colD = colD(:);
    valD = valD(:);

    I = [rowD; rowL; rowU];
    J = [colD; colL; colU];
    V = [valD; valL; valU];

    %% Build sparse matrix
    M = sparse(I, J, V, N2, N2);
end

