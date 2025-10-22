function [s1, s2] = AssembleLinearSystemDual(x1, x2, K11Int, K12Int, K22Int)
% AssembleLinearSystemDual (Conservative finite-volume, node-centered)
%   Solves div( K grad s_k ) = 0 on [0,1]^2 with mixed BCs via face flux balance.
% BCs:
%   s1: s=0 on left, s=1 on right; Neumann(0) on bottom & top
%   s2: s=0 on bottom, s=1 on top;  Neumann(0) on left & right
%
% INPUT
%   x1, x2  : Ny-by-Nx node coordinate arrays (x1 increases along columns, x2 along rows)
%   K*Int   : @(x,y) interpolants for K11, K12, K22 (face samples)
% OUTPUT
%   s1, s2  : Ny-by-Nx node solutions
%
% Notes
% - Uniform grid assumed.
% - Strictly conservative: sum of outgoing face fluxes per control volume is zero.
% - Face fluxes:
%     East face (normal +x):  n·K∇s = K11 * s_x + K12 * s_y at the face center
%     West face (normal −x):  n·K∇s = −K11 * s_x − K12 * s_y at the face center
%     North face(normal +y):  n·K∇s = K12 * s_x + K22 * s_y at the face center
%     South face(normal −y):  n·K∇s = −K12 * s_x − K22 * s_y at the face center
%   Flux = (n·K∇s) * face_length; Sum(E+W+N+S)=0.
%
    [Ny, Nx] = size(x1);
    if any(size(x2) ~= [Ny, Nx]), error('x1 and x2 must be same size'); end
    if Nx < 2 || Ny < 2, error('Grid must be at least 2x2'); end

    % Uniform steps and face lengths
    hx = x1(1,2) - x1(1,1);
    hy = x2(2,1) - x2(1,1);
    Le = hy; Lw = hy; Ln = hx; Ls = hx;  % face lengths

    % Face centers (for coefficient sampling)
    xe = 0.5*(x1(:,1:end-1) + x1(:,2:end));   ye = 0.5*(x2(:,1:end-1) + x2(:,2:end));  % Ny x (Nx-1)
    xn = 0.5*(x1(1:end-1,:) + x1(2:end,:));   yn = 0.5*(x2(1:end-1,:) + x2(2:end,:));  % (Ny-1) x Nx

    % Sample K at faces
    K11e = K11Int(xe, ye);  K12e = K12Int(xe, ye);  K22e = K22Int(xe, ye);
    K11n = K11Int(xn, yn);  K12n = K12Int(xn, yn);  K22n = K22Int(xn, yn);

    % Solve the two scalar problems with their BC sets
    s1 = solve_scalar('s1', x1, x2, hx, hy, ...
                      K11e, K12e, K22e, K11n, K12n, K22n, ...
                      Le, Lw, Ln, Ls);
    s2 = solve_scalar('s2', x1, x2, hx, hy, ...
                      K11e, K12e, K22e, K11n, K12n, K22n, ...
                      Le, Lw, Ln, Ls);
end

% -------------------------------------------------------------------------
function S = solve_scalar(which, x1, x2, hx, hy, ...
                          K11e, K12e, K22e, K11n, K12n, K22n, ...
                          Le, Lw, Ln, Ls)
% Assemble face-by-face conservative FV operator for one scalar field.
    [Ny, Nx] = size(x1);
    id = zeros(Ny, Nx);   % unknown mapping

    % Unknown set: exclude Dirichlet sides only (Neumann sides remain unknown)
    isUnknown = true(Ny, Nx);
    switch which
        case 's1' % Dirichlet on left/right
            isUnknown(:,1)   = false;
            isUnknown(:,end) = false;
        case 's2' % Dirichlet on bottom/top
            isUnknown(1,:)   = false;
            isUnknown(end,:) = false;
        otherwise
            error('which must be ''s1'' or ''s2''');
    end
    id(isUnknown) = 1:nnz(isUnknown);
    N = nnz(isUnknown);

    % Triplet storage (9-pt stencil worst case)
    I = zeros(9*N,1); J = zeros(9*N,1); V = zeros(9*N,1); nz = 0;
    rhs = zeros(N,1);

    % helpers
    function addA(r,c,val)
        nz = nz + 1; I(nz)=r; J(nz)=c; V(nz)=val;
    end
    function add_to_row(row, i2, j2, coeff)
        if coeff==0, return; end
        if isUnknown(i2,j2)
            addA(row, id(i2,j2), coeff);
        else
            g = dirichlet_value(which, i2, j2, Nx, Ny);  % must be Dirichlet on that node
            rhs(row) = rhs(row) - coeff * g;
        end
    end

    % Main loop: one CV per unknown node
    for i = 1:Ny
        for j = 1:Nx
            if ~isUnknown(i,j), continue; end
            row = id(i,j);

            % Accumulate coefficients for this row (start at zero)
            diagC = 0.0;

            % -------------------- EAST face (between j and j+1), n=(+1,0)
            if j < Nx
                k11 = K11e(i,j); k12 = K12e(i,j); % k22 not needed in x-face normal
                % s_x at east face
                cxE = Le * k11 / hx;
                % s_y at east face (centered, uses rows i±1 at columns j and j+1)
                syE = Le * k12 / (4*hy);

                % Coeffs: +cxE*S(i,j+1) - cxE*S(i,j)
                if isDirichletNode(which, i, j+1, Nx, Ny)
                    rhs(row) = rhs(row) - cxE * dirichlet_value(which, i, j+1, Nx, Ny);
                else
                    add_to_row(row, i, j+1, +cxE);
                end
                diagC = diagC - cxE;

                % s_y term: +syE*[ S(i+1,j+1)-S(i-1,j+1)+S(i+1,j)-S(i-1,j) ]
                im = max(1, i-1); ip = min(Ny, i+1);
                % (i+1, j+1)
                if j+1<=Nx, add_to_row(row, ip, j+1, +syE); end
                % (i-1, j+1)
                if j+1<=Nx, add_to_row(row, im, j+1, -syE); end
                % (i+1, j)
                add_to_row(row, ip, j, +syE);
                % (i-1, j)
                add_to_row(row, im, j, -syE);
            else
                % Right boundary face:
                if strcmp(which,'s1') % Dirichlet s1=1 at right -> flux uses boundary value
                    % n·K∇s at right boundary face with s=1
                    k11 = K11e(i,j-1); k12 = K12e(i,j-1); %#ok<NASGU>
                    % To keep strict Dirichlet by rows/cols elimination, we handle it at node level:
                    % rightmost node is not unknown, so this block is not reached for unknown rows.
                else % which='s2', Neumann(0) at right => zero flux: skip
                end
            end

            % -------------------- WEST face (between j-1 and j), n=(−1,0)
            if j > 1
                k11 = K11e(i,j-1); k12 = K12e(i,j-1);
                % For west face: n·K∇s = −K11*s_x − K12*s_y
                cxW = Lw * (-k11) / hx;
                syW = Lw * (-k12) / (4*hy);

                % s_x: cxW*[ S(i,j) - S(i,j-1) ]
                diagC = diagC + cxW;
                if isDirichletNode(which, i, j-1, Nx, Ny)
                    rhs(row) = rhs(row) - (-cxW) * dirichlet_value(which,i,j-1,Nx,Ny);
                else
                    add_to_row(row, i, j-1, -cxW);
                end

                % s_y at west: syW*[ S(i+1,j)-S(i-1,j)+S(i+1,j-1)-S(i-1,j-1) ]
                im = max(1, i-1); ip = min(Ny, i+1);
                add_to_row(row, ip, j,    +syW);
                add_to_row(row, im, j,    -syW);
                add_to_row(row, ip, j-1,  +syW);
                add_to_row(row, im, j-1,  -syW);
            else
                % Left boundary
                if strcmp(which,'s1') % Dirichlet s1=0 at left
                    % leftmost node not unknown -> this block not reached for unknown rows
                else % which='s2', Neumann(0) at left => zero flux: skip
                end
            end

            % -------------------- NORTH face (between i and i+1), n=(0,+1)
            if i < Ny
                k12 = K12n(i,j); k22 = K22n(i,j);
                % n·K∇s = K12*s_x + K22*s_y
                sxN = Ln * k12 / (4*hx);
                cyN = Ln * k22 / hy;

                % s_y: cyN*[ S(i+1,j) - S(i,j) ]
                if isDirichletNode(which, i+1, j, Nx, Ny)
                    rhs(row) = rhs(row) - cyN * dirichlet_value(which, i+1, j, Nx, Ny);
                else
                    add_to_row(row, i+1, j, +cyN);
                end
                diagC = diagC - cyN;

                % s_x at north: sxN*[ S(i+1,j) - S(i+1,j-1) + S(i,j) - S(i,j-1) ]
                jm = max(1, j-1); % clamp for boundary in tangential direction
                add_to_row(row, i+1, j, +sxN);
                add_to_row(row, i+1, jm, -sxN);
                add_to_row(row, i,   j, +sxN);
                add_to_row(row, i,   jm, -sxN);
            else
                % Top boundary
                if strcmp(which,'s2') % Dirichlet s2=1 at top
                    % top nodes not unknown -> not reached
                else % which='s1', Neumann(0) at top => zero flux: skip
                end
            end

            % -------------------- SOUTH face (between i-1 and i), n=(0,−1)
            if i > 1
                k12 = K12n(i-1,j); k22 = K22n(i-1,j);
                % n·K∇s = −K12*s_x − K22*s_y
                sxS = Ls * (-k12) / (4*hx);
                cyS = Ls * (-k22) / hy;

                % s_y: cyS*[ S(i,j) - S(i-1,j) ]
                diagC = diagC + cyS;
                if isDirichletNode(which, i-1, j, Nx, Ny)
                    rhs(row) = rhs(row) - (-cyS) * dirichlet_value(which, i-1, j, Nx, Ny);
                else
                    add_to_row(row, i-1, j, -cyS);
                end

                % s_x at south: sxS*[ S(i,j) - S(i,j-1) + S(i-1,j) - S(i-1,j-1) ]
                jm = max(1, j-1);
                add_to_row(row, i,   j,   +sxS);
                add_to_row(row, i,   jm,  -sxS);
                add_to_row(row, i-1, j,   +sxS);
                add_to_row(row, i-1, jm,  -sxS);
            else
                % Bottom boundary
                if strcmp(which,'s2') % Dirichlet s2=0 at bottom
                    % bottom nodes not unknown -> not reached
                else % which='s1', Neumann(0) at bottom => zero flux: skip
                end
            end

            % Put diagonal
            addA(row, id(i,j), diagC);
        end
    end

    % Build sparse and solve
    I = I(1:nz); J = J(1:nz); V = V(1:nz);
    A = sparse(I, J, V, N, N);
    u = A \ rhs;

    % Scatter to full grid and set Dirichlet nodes
    S = zeros(Ny, Nx);
    S(isUnknown) = u;
    switch which
        case 's1'
            S(:,1)   = 0.0;  % left
            S(:,end) = 1.0;  % right
        case 's2'
            S(1,:)   = 0.0;  % bottom
            S(end,:) = 1.0;  % top
    end
end

% -------------------------------------------------------------------------
function tf = isDirichletNode(which, i, j, Nx, Ny)
% True if the node (i,j) lies on a Dirichlet boundary for the given scalar.
    switch which
        case 's1', tf = (j==1) || (j==Nx);
        case 's2', tf = (i==1) || (i==Ny);
        otherwise,  tf = false;
    end
end

function g = dirichlet_value(which, i, j, Nx, Ny)
% Dirichlet values as specified.
    switch which
        case 's1'
            if j==1,   g = 0.0; elseif j==Nx, g = 1.0; else, g = NaN; end
        case 's2'
            if i==1,   g = 0.0; elseif i==Ny, g = 1.0; else, g = NaN; end
        otherwise
            g = NaN;
    end
end



%{
function [s1, s2] = AssembleLinearSystemDual(x1, x2, K11Int, K12Int, K22Int)
% AssembleLinearSystemDual
% Conservative FD (finite-volume) solver for the dual inverse-map system:
%   div( K * grad s_k ) = 0,  k = 1,2,  on [0,1]^2
% Boundary conditions:
%   s1: Dirichlet 0 on left (j=1), Dirichlet 1 on right (j=Nx),
%       Neumann (zero normal derivative) on top/bottom (i=1 and i=Ny).
%   s2: Dirichlet 0 on bottom (i=1), Dirichlet 1 on top (i=Ny),
%       Neumann (zero normal derivative) on left/right (j=1 and j=Nx).
%
% INPUT:
%   x1, x2 : Ny-by-Nx uniform grids (x1 increases along columns, x2 along rows)
%            e.g., [x1,x2] = meshgrid(linspace(0,1,Nx), linspace(0,1,Ny));
%   K11Int, K12Int, K22Int : function handles @(x,y) for components of K (symmetric)
% OUTPUT:
%   s1, s2 : Ny-by-Nx solutions on the nodes.

    [Ny, Nx] = size(x1);
    if any(size(x2) ~= [Ny, Nx]), error('x1 and x2 must be same size'); end

    % Uniform spacings (assumed)
    if Nx < 2 || Ny < 2, error('Grid must be at least 2x2'); end
    hx = x1(1,2) - x1(1,1);
    hy = x2(2,1) - x2(1,1);

    % Face-center coordinates for flux evaluation
    % East/West faces: centers between columns j and j+1 (size Ny x (Nx-1))
    xe = 0.5*(x1(:,1:end-1) + x1(:,2:end));
    ye = 0.5*(x2(:,1:end-1) + x2(:,2:end));
    % North/South faces: centers between rows i and i+1 (size (Ny-1) x Nx)
    xn = 0.5*(x1(1:end-1,:) + x1(2:end,:));
    yn = 0.5*(x2(1:end-1,:) + x2(2:end,:));

    % Sample K at face centers
    K11e = K11Int(xe, ye);  K12e = K12Int(xe, ye);  K22e = K22Int(xe, ye);   % Ny x (Nx-1)
    K11n = K11Int(xn, yn);  K12n = K12Int(xn, yn);  K22n = K22Int(xn, yn);   % (Ny-1) x Nx

    % Assemble and solve for s1 and s2 with their respective BCs
    s1 = solve_scalar('s1', Nx, Ny, hx, hy, x1, x2, K11e, K12e, K22e, K11n, K12n, K22n);
    s2 = solve_scalar('s2', Nx, Ny, hx, hy, x1, x2, K11e, K12e, K22e, K11n, K12n, K22n);
end

% -------------------------------------------------------------------------
function S = solve_scalar(which, Nx, Ny, hx, hy, x1, x2, K11e, K12e, K22e, K11n, K12n, K22n)
% Build conservative linear system div(K grad S) = 0 with mixed BCs.
%
% which = 's1' or 's2' to select boundary conditions.

    % Unknown mask depends on where Dirichlet BCs are applied
    isUnknown = true(Ny, Nx);
    switch which
        case 's1'
            % s1: Dirichlet on left/right columns
            isUnknown(:,1)   = false;
            isUnknown(:,end) = false;
        case 's2'
            % s2: Dirichlet on bottom/top rows
            isUnknown(1,:)   = false;
            isUnknown(end,:) = false;
        otherwise
            error('which must be s1 or s2');
    end

    % Mapping from (i,j) -> eqn index
    id = zeros(Ny, Nx);
    id(isUnknown) = 1:nnz(isUnknown);
    N = nnz(isUnknown);

    I = zeros(9*N,1); J = zeros(9*N,1); V = zeros(9*N,1);  % over-alloc (9-pt stencil)
    rhs = zeros(N,1);
    nz  = 0;   % fill pointer

    % Helper to add matrix entry
    function add_entry(ri, cj, val)
        nz = nz + 1;
        I(nz) = ri; J(nz) = cj; V(nz) = val;
    end
    % Helper to add to RHS
    function add_rhs(ri, val), rhs(ri) = rhs(ri) + val; end

    % Dirichlet boundary value function
    function val = dirichlet_value(i,j)
        switch which
            case 's1'
                if j==1, val = 0.0; elseif j==Nx, val = 1.0; else, val = NaN; end
            case 's2'
                if i==1, val = 0.0; elseif i==Ny, val = 1.0; else, val = NaN; end
        end
    end

    % Loop over unknown control volumes (nodes)
    for i = 1:Ny
        for j = 1:Nx
            if ~isUnknown(i,j), continue; end
            row = id(i,j);

            aP = 0.0; % central coefficient

            % -----------------------------
            % EAST face (between j and j+1)
            if j < Nx
                % Face K at (i, j+1/2)
                k11 = K11e(i,j); k12 = K12e(i,j); k22 = K22e(i,j);

                % Normal derivative approx at face center
                % s_x at east face:
                %   (S(i,j+1) - S(i,j)) / hx
                % Tangential derivative s_y at east face (average of the two columns):
                %   ((S(ip,j) - S(im,j)) + (S(ip,j+1) - S(im,j+1))) / (4*hy)
                im = (i>1)* (i-1) + (i==1)*i;      % clamp for Neumann at bottom
                ip = (i<Ny)* (i+1) + (i==Ny)*i;    % clamp for Neumann at top

                % Flux through EAST face: F_e = k11*s_x + k12*s_y
                % Contribution pattern (nine-point):
                cE =  k11/hx;  cC = -k11/hx;    % from s_x term

                % s_y contributions (1/(4hy)) * k12 * [ +S(ip,j) - S(im,j) + S(ip,j+1) - S(im,j+1) ]
                sy_fac = k12/(4*hy);

                % Accumulate to matrix / rhs for EAST face (added with + sign)
                % Handle Dirichlet at (i,j+1)
                if j+1<=Nx && isUnknown(i,j+1)
                    add_entry(row, id(i,j+1),  cE);  aP = aP - cE;
                else
                    val = dirichlet_value(i,j+1);
                    if ~isnan(val), add_rhs(row,  cE*val); else, % Neumann side (for s2)
                        % No Dirichlet there, but if Neumann was specified on this side,
                        % do nothing extra (unknown remains if applicable). For s1, sides are Dirichlet already.
                    end
                end
                % Center node
                add_entry(row, id(i,j), cC);      aP = aP - cC;

                % s_y at EAST face uses four nodes: (im,j), (ip,j), (im,j+1), (ip,j+1)
                % (im,j)
                if isUnknown(im,j)
                    add_entry(row, id(im,j), -sy_fac);
                else
                    val = dirichlet_value(im,j);
                    if ~isnan(val), add_rhs(row, -sy_fac*val); end
                end
                % (ip,j)
                if isUnknown(ip,j)
                    add_entry(row, id(ip,j),  sy_fac);
                else
                    val = dirichlet_value(ip,j);
                    if ~isnan(val), add_rhs(row,  sy_fac*val); end
                end
                % (im,j+1)
                if j+1<=Nx && isUnknown(im,j+1)
                    add_entry(row, id(im,j+1), -sy_fac);
                else
                    if j+1<=Nx
                        val = dirichlet_value(im,j+1);
                        if ~isnan(val), add_rhs(row, -sy_fac*val); end
                    end
                end
                % (ip,j+1)
                if j+1<=Nx && isUnknown(ip,j+1)
                    add_entry(row, id(ip,j+1),  sy_fac);
                else
                    if j+1<=Nx
                        val = dirichlet_value(ip,j+1);
                        if ~isnan(val), add_rhs(row,  sy_fac*val); end
                    end
                end
            end

            % -----------------------------
            % WEST face (between j-1 and j)
            if j > 1
                % Face K at (i, j-1/2)   -> indices (i, j-1) in *_e arrays
                k11 = K11e(i,j-1); k12 = K12e(i,j-1); k22 = K22e(i,j-1);

                % s_x at west face: (S(i,j) - S(i,j-1))/hx
                cC =  k11/hx;  cW = -k11/hx;

                im = (i>1)* (i-1) + (i==1)*i;
                ip = (i<Ny)* (i+1) + (i==Ny)*i;
                sy_fac = k12/(4*hy); % same form, now for west column pair (j-1, j)

                % WEST face enters with minus sign in divergence: -(F_w)
                % So we ADD (-cW) to neighbor, (-cC) to center, etc.

                % (i,j-1)
                if isUnknown(i,j-1)
                    add_entry(row, id(i,j-1), -cW);  aP = aP + cW;
                else
                    val = dirichlet_value(i,j-1);
                    if ~isnan(val), add_rhs(row, -cW*val); end
                end
                % center
                add_entry(row, id(i,j), -cC);       aP = aP + cC;

                % s_y at WEST face: +S(ip,j-1)-S(im,j-1)+S(ip,j)-S(im,j)
                % multiplied by sy_fac, then whole WEST flux is subtracted.
                % => coefficients get a negative sign.
                % (im,j-1)
                if isUnknown(im,j-1)
                    add_entry(row, id(im,j-1),  sy_fac);   % minus * ( -1 ) = +sy_fac
                else
                    val = dirichlet_value(im,j-1);
                    if ~isnan(val), add_rhs(row,  sy_fac*val); end
                end
                % (ip,j-1)
                if isUnknown(ip,j-1)
                    add_entry(row, id(ip,j-1), -sy_fac);   % minus * (+1) = -sy_fac
                else
                    val = dirichlet_value(ip,j-1);
                    if ~isnan(val), add_rhs(row, -sy_fac*val); end
                end
                % (im,j)
                if isUnknown(im,j)
                    add_entry(row, id(im,j),  sy_fac);
                else
                    val = dirichlet_value(im,j);
                    if ~isnan(val), add_rhs(row,  sy_fac*val); end
                end
                % (ip,j)
                if isUnknown(ip,j)
                    add_entry(row, id(ip,j), -sy_fac);
                else
                    val = dirichlet_value(ip,j);
                    if ~isnan(val), add_rhs(row, -sy_fac*val); end
                end
            end

            % -----------------------------
            % NORTH face (between i and i+1)
            if i < Ny
                % Face K at (i+1/2, j) -> indices (i, j) in *_n arrays
                k11 = K11n(i,j); k12 = K12n(i,j); k22 = K22n(i,j);

                % s_y at north face: (S(i+1,j) - S(i,j))/hy
                cN =  k22/hy;  cC = -k22/hy;

                jm = (j>1)* (j-1) + (j==1)*j;      % clamp for Neumann at left
                jp = (j<Nx)* (j+1) + (j==Nx)*j;    % clamp for Neumann at right

                % tangential derivative s_x at north face (avg of the two rows):
                % ((S(i,jp)-S(i,jm)) + (S(i+1,jp)-S(i+1,jm))) / (4*hx)
                sx_fac = k12/(4*hx);

                % NORTH face enters with plus sign
                % (i+1,j)
                if isUnknown(i+1,j)
                    add_entry(row, id(i+1,j),  cN);   aP = aP - cN;
                else
                    val = dirichlet_value(i+1,j);
                    if ~isnan(val), add_rhs(row,  cN*val); end
                end
                % center
                add_entry(row, id(i,j), cC);         aP = aP - cC;

                % s_x contributors: (i,jm), (i,jp), (i+1,jm), (i+1,jp)
                % (i,jm)
                if isUnknown(i,jm)
                    add_entry(row, id(i,jm), -sx_fac);
                else
                    val = dirichlet_value(i,jm);
                    if ~isnan(val), add_rhs(row, -sx_fac*val); end
                end
                % (i,jp)
                if isUnknown(i,jp)
                    add_entry(row, id(i,jp),  sx_fac);
                else
                    val = dirichlet_value(i,jp);
                    if ~isnan(val), add_rhs(row,  sx_fac*val); end
                end
                % (i+1,jm)
                if isUnknown(i+1,jm)
                    add_entry(row, id(i+1,jm), -sx_fac);
                else
                    val = dirichlet_value(i+1,jm);
                    if ~isnan(val), add_rhs(row, -sx_fac*val); end
                end
                % (i+1,jp)
                if isUnknown(i+1,jp)
                    add_entry(row, id(i+1,jp),  sx_fac);
                else
                    val = dirichlet_value(i+1,jp);
                    if ~isnan(val), add_rhs(row,  sx_fac*val); end
                end
            end

            % -----------------------------
            % SOUTH face (between i-1 and i)
            if i > 1
                % Face K at (i-1/2, j) -> indices (i-1, j) in *_n arrays
                k11 = K11n(i-1,j); k12 = K12n(i-1,j); k22 = K22n(i-1,j);

                % s_y at south face: (S(i,j) - S(i-1,j))/hy
                cC =  k22/hy;  cS = -k22/hy;

                jm = (j>1)* (j-1) + (j==1)*j;
                jp = (j<Nx)* (j+1) + (j==Nx)*j;
                sx_fac = k12/(4*hx);

                % SOUTH face enters with minus sign in divergence
                % Neighbor (i-1,j)
                if isUnknown(i-1,j)
                    add_entry(row, id(i-1,j), -cS);  aP = aP + cS;
                else
                    val = dirichlet_value(i-1,j);
                    if ~isnan(val), add_rhs(row, -cS*val); end
                end
                % center
                add_entry(row, id(i,j), -cC);        aP = aP + cC;

                % s_x contributors in SOUTH (minus sign applied):
                % (i-1,jm), (i-1,jp), (i,jm), (i,jp)
                % (i-1,jm)
                if isUnknown(i-1,jm)
                    add_entry(row, id(i-1,jm),  sx_fac);
                else
                    val = dirichlet_value(i-1,jm);
                    if ~isnan(val), add_rhs(row,  sx_fac*val); end
                end
                % (i-1,jp)
                if isUnknown(i-1,jp)
                    add_entry(row, id(i-1,jp), -sx_fac);
                else
                    val = dirichlet_value(i-1,jp);
                    if ~isnan(val), add_rhs(row, -sx_fac*val); end
                end
                % (i,jm)
                if isUnknown(i,jm)
                    add_entry(row, id(i,jm),  sx_fac);
                else
                    val = dirichlet_value(i,jm);
                    if ~isnan(val), add_rhs(row,  sx_fac*val); end
                end
                % (i,jp)
                if isUnknown(i,jp)
                    add_entry(row, id(i,jp), -sx_fac);
                else
                    val = dirichlet_value(i,jp);
                    if ~isnan(val), add_rhs(row, -sx_fac*val); end
                end
            end

            % Finally, put the central coefficient aP on the diagonal row balance
            % (Note: aP should sum the negative of all off-diagonal contributions
            %  added via the normal-derivative terms above.)
            add_entry(row, id(i,j), -aP);
        end
    end

    % Shrink triplets and build sparse matrix
    I = I(1:nz); J = J(1:nz); V = V(1:nz);
    A = sparse(I, J, V, N, N);

    % Solve
    u = A \ rhs;

    % Scatter back to full grid
    S = zeros(Ny, Nx);
    S(~isUnknown) = 0;  % will overwrite Dirichlet nodes explicitly next
    S(isUnknown) = u;

    % Set Dirichlet boundaries explicitly
    switch which
        case 's1'
            S(:,1)   = 0.0;
            S(:,end) = 1.0;
        case 's2'
            S(1,:)   = 0.0;
            S(end,:) = 1.0;
    end
end
%}
