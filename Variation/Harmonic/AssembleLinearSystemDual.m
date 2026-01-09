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

    % All nodes are unknowns now
    N  = Ny * Nx;
    id = reshape(1:N, Ny, Nx);  % mapping (i,j) -> global index

    % Masks for Dirichlet and Neumann nodes
    isDirichlet = false(Ny, Nx);
    isNeumann   = false(Ny, Nx);

    switch which
        case 's1'
            % Dirichlet on left/right
            isDirichlet(:,1)   = true;
            isDirichlet(:,end) = true;

            % Neumann(0) on bottom/top (exclude corners so they stay pure Dirichlet)
            if Nx > 2
                isNeumann(1,   2:Nx-1) = true;  % bottom
                isNeumann(Ny,  2:Nx-1) = true;  % top
            end

        case 's2'
            % Dirichlet on bottom/top
            isDirichlet(1,:)   = true;
            isDirichlet(end,:) = true;

            % Neumann(0) on left/right (exclude corners)
            if Ny > 2
                isNeumann(2:Ny-1, 1  ) = true; % left
                isNeumann(2:Ny-1, Nx ) = true; % right
            end

        otherwise
            error('which must be ''s1'' or ''s2''');
    end

    % Triplet storage (9-pt stencil worst case)
    I = zeros(9*N,1); J = zeros(9*N,1); V = zeros(9*N,1); nz = 0;
    rhs = zeros(N,1);

    % helper to add matrix entry
    function addA(r,c,val)
        nz = nz + 1;
        I(nz) = r;
        J(nz) = c;
        V(nz) = val;
    end

    % Main loop: one CV / equation per node
    for i = 1:Ny
        for j = 1:Nx
            row = id(i,j);

            % ================= DIRICHLET ROW =================
            if isDirichlet(i,j)
                % s(i,j) = g(i,j)
                addA(row, row, 1.0);
                rhs(row) = dirichlet_value(which, i, j, Nx, Ny);
                continue;
            end
            % =================================================

            % ================= NEUMANN ROW ===================
            if isNeumann(i,j)
                % A - B = 0, where B is neighbor in normal direction
                switch which
                    case 's1'
                        % Neumann on bottom/top, normal ±y
                        if i == 1        % bottom -> neighbor above
                            in = i + 1; jn = j;
                        elseif i == Ny   % top -> neighbor below
                            in = i - 1; jn = j;
                        else
                            error('isNeumann(i,j) true at non-Neumann row for s1');
                        end
                    case 's2'
                        % Neumann on left/right, normal ±x
                        if j == 1        % left -> neighbor to the right
                            in = i; jn = j + 1;
                        elseif j == Nx   % right -> neighbor to the left
                            in = i; jn = j - 1;
                        else
                            error('isNeumann(i,j) true at non-Neumann column for s2');
                        end
                    otherwise
                        error('Unknown scalar field');
                end

                addA(row, id(i,j),  1.0);   % A
                addA(row, id(in,jn), -1.0); % -B
                % RHS = 0 by default
                continue;
            end
            % =================================================

            % Interior / non-BC node: assemble conservative flux balance
            diagC = 0.0;

            % -------------------- EAST face (between j and j+1), n=(+1,0)
            if j < Nx
                k11 = K11e(i,j); k12 = K12e(i,j);
                cxE = Le * k11 / hx;
                syE = Le * k12 / (4*hy);

                % s_x: +cxE*S(i,j+1) - cxE*S(i,j)
                addA(row, id(i, j+1), +cxE);
                diagC = diagC - cxE;

                % s_y term: +syE*[ S(i+1,j+1)-S(i-1,j+1)+S(i+1,j)-S(i-1,j) ]
                im = max(1,  i-1);
                ip = min(Ny, i+1);
                addA(row, id(ip, j+1), +syE);
                addA(row, id(im, j+1), -syE);
                addA(row, id(ip, j   ), +syE);
                addA(row, id(im, j   ), -syE);
            end

            % -------------------- WEST face (between j-1 and j), n=(−1,0)
            if j > 1
                k11 = K11e(i,j-1); k12 = K12e(i,j-1);
                % n·K∇s = −K11*s_x − K12*s_y
                cxW = Lw * (-k11) / hx;
                syW = Lw * (-k12) / (4*hy);

                % s_x: cxW*[ S(i,j) - S(i,j-1) ]
                diagC = diagC + cxW;
                addA(row, id(i, j-1), -cxW);

                % s_y: syW*[ S(i+1,j)-S(i-1,j)+S(i+1,j-1)-S(i-1,j-1) ]
                im = max(1,  i-1);
                ip = min(Ny, i+1);
                addA(row, id(ip, j  ),  +syW);
                addA(row, id(im, j  ),  -syW);
                addA(row, id(ip, j-1),  +syW);
                addA(row, id(im, j-1),  -syW);
            end

            % -------------------- NORTH face (between i and i+1), n=(0,+1)
            if i < Ny
                k12 = K12n(i,j); k22 = K22n(i,j);
                % n·K∇s = K12*s_x + K22*s_y
                sxN = Ln * k12 / (4*hx);
                cyN = Ln * k22 / hy;

                % s_y: cyN*[ S(i+1,j) - S(i,j) ]
                addA(row, id(i+1, j), +cyN);
                diagC = diagC - cyN;

                % s_x: sxN*[ S(i+1,j) - S(i+1,j-1) + S(i,j) - S(i,j-1) ]
                jm = max(1, j-1);
                addA(row, id(i+1, j ), +sxN);
                addA(row, id(i+1, jm), -sxN);
                addA(row, id(i,   j ), +sxN);
                addA(row, id(i,   jm), -sxN);
            end

            % -------------------- SOUTH face (between i-1 and i), n=(0,−1)
            if i > 1
                k12 = K12n(i-1,j); k22 = K22n(i-1,j);
                % n·K∇s = −K12*s_x − K22*s_y
                sxS = Ls * (-k12) / (4*hx);
                cyS = Ls * (-k22) / hy;

                % s_y: cyS*[ S(i,j) - S(i-1,j) ]
                diagC = diagC + cyS;
                addA(row, id(i-1, j), -cyS);

                % s_x: sxS*[ S(i,j) - S(i,j-1) + S(i-1,j) - S(i-1,j-1) ]
                jm = max(1, j-1);
                addA(row, id(i,   j ),  +sxS);
                addA(row, id(i,   jm),  -sxS);
                addA(row, id(i-1, j ),  +sxS);
                addA(row, id(i-1, jm),  -sxS);
            end

            % Put diagonal contribution
            addA(row, id(i,j), diagC);
        end
    end

    % Build sparse and solve
    I = I(1:nz); J = J(1:nz); V = V(1:nz);
    A = sparse(I, J, V, N, N);

    u = A \ rhs;

    % Reshape back to grid
    S = reshape(u, Ny, Nx);

    % (Optional) enforce exact Dirichlet values again for cleanliness
    switch which
        case 's1'
            S(:,1)   = 0.0;
            S(:,end) = 1.0;
        case 's2'
            S(1,:)   = 0.0;
            S(end,:) = 1.0;
    end
end


%{
function S = solve_scalar(which, x1, x2, hx, hy, ...
                          K11e, K12e, K22e, K11n, K12n, K22n, ...
                          Le, Lw, Ln, Ls)
% Assemble face-by-face conservative FV operator for one scalar field.
    [Ny, Nx] = size(x1);
    id = zeros(Ny, Nx);   % unknown mapping

    % Unknown set: exclude Dirichlet sides only (Neumann sides remain unknown)
    isUnknown  = true(Ny, Nx);
    isNeumann  = false(Ny, Nx);  % mark Neumann boundary nodes (unknowns only)

    switch which
        case 's1' % Dirichlet on left/right; Neumann(0) on bottom/top
            isUnknown(:,1)   = false;
            isUnknown(:,end) = false;

            % Neumann boundaries: bottom and top, interior columns only
            if Ny > 1 && Nx > 2
                isNeumann(1,   2:Nx-1) = true;  % bottom edge (excluding corners)
                isNeumann(Ny,  2:Nx-1) = true;  % top edge    (excluding corners)
            end

        case 's2' % Dirichlet on bottom/top; Neumann(0) on left/right
            isUnknown(1,:)   = false;
            isUnknown(end,:) = false;

            % Neumann boundaries: left and right, interior rows only
            if Nx > 1 && Ny > 2
                isNeumann(2:Ny-1, 1   ) = true; % left edge (excluding corners)
                isNeumann(2:Ny-1, Nx  ) = true; % right edge (excluding corners)
            end

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

            % ================== EXPLICIT NEUMANN ROW ==================
            if isNeumann(i,j)
                % Rewrite row as A - B = 0, where B is neighbor in normal direction
                switch which
                    case 's1'
                        % Neumann on bottom/top edges, normal in ±y
                        if i == 1        % bottom: neighbor just above
                            in = i + 1; jn = j;
                        elseif i == Ny   % top: neighbor just below
                            in = i - 1; jn = j;
                        else
                            error('isNeumann(i,j) true at non-Neumann row for s1');
                        end
                    case 's2'
                        % Neumann on left/right edges, normal in ±x
                        if j == 1        % left: neighbor just to the right
                            in = i; jn = j + 1;
                        elseif j == Nx   % right: neighbor just to the left
                            in = i; jn = j - 1;
                        else
                            error('isNeumann(i,j) true at non-Neumann column for s2');
                        end
                    otherwise
                        error('Unknown scalar field');
                end

                % A - B = 0
                addA(row, id(i,j), 1.0);
                add_to_row(row, in, jn, -1.0);

                % Skip flux-balance assembly for this Neumann node
                continue;
            end
            % ==========================================================

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
                    % Node is Dirichlet (not an unknown), so we don't get here.
                else % which='s2', Neumann(0) at right: now handled by explicit row if unknown; here j==Nx for interior i, but those i are Neumann and caught above.
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
                if strcmp(which,'s1')
                    % Dirichlet s1=0 at left: node not unknown
                else
                    % which='s2', Neumann(0) at left: unknown boundary node now handled
                    % by explicit A-B=0 row above when j==1.
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
                if strcmp(which,'s2')
                    % Dirichlet s2=1 at top: node not unknown
                else
                    % which='s1', Neumann(0) at top: unknown boundary node handled
                    % by explicit A-B=0 row above when i==Ny.
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
                if strcmp(which,'s2')
                    % Dirichlet s2=0 at bottom: node not unknown
                else
                    % which='s1', Neumann(0) at bottom: unknown boundary node handled
                    % by explicit A-B=0 row above when i==1.
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

%}


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
