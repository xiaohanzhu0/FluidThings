%% separable_metric_meshgen_fixedpoint.m
% Separable (perfectly orthogonal) mesh generation on [0,1]^2
% using alternating projected Newton/LM on the discrete EL objective:
%
%   J1[x1] = ∫_0^1 ∫_0^1 (Delta1^2 M11(x1(s1), x2(s2)) (x1')^2 - 1)^2 ds2 ds1
%   J2[x2] = ∫_0^1 ∫_0^1 (Delta2^2 M22(x1(s1), x2(s2)) (x2')^2 - 1)^2 ds1 ds2
%
% The coupling enters through the inner integral over the other coordinate.
% Each subproblem is solved with a damped Newton step and monotonicity projection.
%
% REQUIREMENTS:
%   Define function handles M11(x1,x2), M22(x1,x2) in workspace or below.
%
% Notes:
%   - Under the separable orthogonal ansatz, M12 does not enter.
%   - We include simple damping to help preserve monotonicity (no folding).
%   - Delta1/Delta2 updated each outer iteration using ∬ a / ∬ a^2 style.

clear; clc;

%% ---------------- User-provided metric handles ----------------
% Example (uncomment to test):
M11 = @(x1,x2) 40000 * (1 + 15*x1).^(-2);
M22 = @(x1,x2) 40000 * (1 + 15*x2).^(-2);

if ~exist('M11','var') || ~exist('M22','var')
    error('Please define function handles M11 and M22 before running.');
end

%% ---------------- Parameters ----------------
N1 = 51; N2 = 51;
h1 = 1/(N1-1);
h2 = 1/(N2-1);

% Quadrature weights (composite trapezoid) on [0,1]
w1 = h1*ones(N1,1); w1(1)=h1/2; w1(end)=h1/2;
w2 = h2*ones(N2,1); w2(1)=h2/2; w2(end)=h2/2;

% Outer alternation
maxOuter = 20;
tolOuter  = 1e-8;

% Newton/LM settings for each 1D subproblem
maxNewton1 = 5;
maxNewton2 = 5;
tolNewton = 1e-8;
gradEps = 1e-6;
lambda0 = 1e-4;
lambdaFactor = 10.0;
maxLM = 10;

% Damping to preserve monotonicity
omega = 0.5;                 % 0<omega<=1, smaller = safer
minSlopeFactor = 1e-8;       % enforce diff(x) > minSlopeFactor*h

%% ---------------- Initialize maps ----------------
s1 = linspace(0,1,N1)';
s2 = linspace(0,1,N2)';

x1 = s1;
x2 = s2;

Delta1 = 1; Delta2 = 1;

Jprev = inf;

%% ---------------- Outer loop (robust EL Newton with coupling) ----------------
for outer = 1:maxOuter
    x1_old_outer = x1;
    x2_old_outer = x2;

    [x1, Delta1] = optimize_x1_newton(x1, x2, M11, w1, w2, h1, ...
        gradEps, maxNewton1, tolNewton, lambda0, lambdaFactor, maxLM, ...
        omega, minSlopeFactor*h1);

    Delta1 = updateDelta1(x1,x2,M11,w1,w2,h1);
    Delta2 = updateDelta2(x1,x2,M22,w1,w2,h2);
    [x2, Delta2] = optimize_x2_newton(x2, x1, M22, w2, w1, h2, ...
        gradEps, maxNewton2, tolNewton, lambda0, lambdaFactor, maxLM, ...
        omega, minSlopeFactor*h2);

    Delta2 = updateDelta2(x1,x2,M22,w1,w2,h2);

    % Monitor objective
    J = fullObjective(x1,x2,M11,M22,Delta1,Delta2,h1,h2);
    relDec = abs(Jprev - J) / max(1,abs(Jprev));

    fprintf('Outer %2d | J=%.6e | relDec=%.3e | ||dx1||=%.3e | ||dx2||=%.3e\n', ...
        outer, J, relDec, norm(x1-x1_old_outer,2), norm(x2-x2_old_outer,2));

    if relDec < tolOuter
        break;
    end
    Jprev = J;
end

%% ---------------- Plot mesh ----------------
[X1, X2] = meshgrid(x1, x2);

figure; hold on; axis equal tight;
plot(X1, X2, 'k'); plot(X1', X2', 'k');
title('Separable orthogonal mesh (EL Newton/LM, coupled)');
xlabel('x_1'); ylabel('x_2');

%% ===================== Local functions =====================

function x = enforceMonotone(x, minDx)
    % Enforce strict monotonicity by "pushing up" any non-increasing steps.
    dx = diff(x);
    if all(dx > minDx), return; end
    for i = 2:numel(x)
        if x(i) <= x(i-1) + minDx
            x(i) = x(i-1) + minDx;
        end
    end
    % Renormalize end point to 1, start to 0 (preserve monotonicity)
    x = (x - x(1)) / max(x(end)-x(1), eps);
    x(1) = 0; x(end) = 1;
end

function [x1, Delta1] = optimize_x1_newton(x1, x2, M11, w1, w2, h1, ...
    gradEps, maxNewton, tolNewton, lambda0, lambdaFactor, maxLM, omega, minSlope)
    for it = 1:maxNewton
        [Jcur, Delta1] = objective_x1(x1, x2, M11, w1, w2, h1);
        grad = grad_fd_x1(x1, x2, M11, w1, w2, h1, gradEps, minSlope);
        g = grad(2:end-1);
        gnorm = norm(g, 2);
        if gnorm < tolNewton * max(1, norm(x1, 2))
            break;
        end

        H = hess_fd_x1(x1, x2, M11, w1, w2, h1, gradEps, minSlope);
        H = 0.5 * (H + H');
        lambda = lambda0;
        accepted = false;
        for lm = 1:maxLM
            d = -(H + lambda * eye(size(H))) \ g;
            x_try = x1;
            x_try(2:end-1) = x_try(2:end-1) + d;
            x_try = enforceMonotone(x_try, minSlope);
            [J_try, Delta_try] = objective_x1(x_try, x2, M11, w1, w2, h1);
            if J_try < Jcur
                accepted = true;
                break;
            end
            lambda = lambda * lambdaFactor;
        end

        if ~accepted
            break;
        end

        x1 = (1-omega) * x1 + omega * x_try;
        x1 = enforceMonotone(x1, minSlope);
        Delta1 = Delta_try;

        if abs(Jcur - J_try) / max(1, Jcur) < tolNewton
            break;
        end
    end
end

function [x2, Delta2] = optimize_x2_newton(x2, x1, M22, w2, w1, h2, ...
    gradEps, maxNewton, tolNewton, lambda0, lambdaFactor, maxLM, omega, minSlope)
    for it = 1:maxNewton
        [Jcur, Delta2] = objective_x2(x2, x1, M22, w2, w1, h2);
        grad = grad_fd_x2(x2, x1, M22, w2, w1, h2, gradEps, minSlope);
        g = grad(2:end-1);
        gnorm = norm(g, 2);
        if gnorm < tolNewton * max(1, norm(x2, 2))
            break;
        end

        H = hess_fd_x2(x2, x1, M22, w2, w1, h2, gradEps, minSlope);
        H = 0.5 * (H + H');
        lambda = lambda0;
        accepted = false;
        for lm = 1:maxLM
            d = -(H + lambda * eye(size(H))) \ g;
            x_try = x2;
            x_try(2:end-1) = x_try(2:end-1) + d;
            x_try = enforceMonotone(x_try, minSlope);
            [J_try, Delta_try] = objective_x2(x_try, x1, M22, w2, w1, h2);
            if J_try < Jcur
                accepted = true;
                break;
            end
            lambda = lambda * lambdaFactor;
        end

        if ~accepted
            break;
        end

        x2 = (1-omega) * x2 + omega * x_try;
        x2 = enforceMonotone(x2, minSlope);
        Delta2 = Delta_try;

        if abs(Jcur - J_try) / max(1, Jcur) < tolNewton
            break;
        end
    end
end

function [J, Delta1] = objective_x1(x1, x2, M11, w1, w2, h1)
    Delta1 = updateDelta1(x1, x2, M11, w1, w2, h1);
    p = d1_oneSidedCentral(x1, h1);
    [X1, X2] = meshgrid(x1, x2);
    M11_vals = M11(X1, X2);
    p2 = repmat((p.^2)', numel(x2), 1);
    r = Delta1^2 * M11_vals .* p2 - 1;
    W = w2 * w1';
    J = sum(sum(W .* (r.^2)));
end

function [J, Delta2] = objective_x2(x2, x1, M22, w2, w1, h2)
    Delta2 = updateDelta2(x1, x2, M22, w1, w2, h2);
    q = d1_oneSidedCentral(x2, h2);
    [X1, X2] = meshgrid(x1, x2);
    M22_vals = M22(X1, X2);
    q2 = repmat((q.^2), 1, numel(x1));
    r = Delta2^2 * M22_vals .* q2 - 1;
    W = w2 * w1';
    J = sum(sum(W .* (r.^2)));
end

function grad = grad_fd_x1(x1, x2, M11, w1, w2, h1, gradEps, minSlope)
    N1 = numel(x1);
    grad = zeros(N1,1);
    for i = 2:N1-1
        x1p = perturb_monotone(x1, i, gradEps, minSlope);
        x1m = perturb_monotone(x1, i, -gradEps, minSlope);
        if x1p(i) == x1m(i)
            grad(i) = 0;
            continue;
        end
        Jp = objective_x1(x1p, x2, M11, w1, w2, h1);
        Jm = objective_x1(x1m, x2, M11, w1, w2, h1);
        grad(i) = (Jp - Jm) / (x1p(i) - x1m(i));
    end
end

function grad = grad_fd_x2(x2, x1, M22, w2, w1, h2, gradEps, minSlope)
    N2 = numel(x2);
    grad = zeros(N2,1);
    for j = 2:N2-1
        x2p = perturb_monotone(x2, j, gradEps, minSlope);
        x2m = perturb_monotone(x2, j, -gradEps, minSlope);
        if x2p(j) == x2m(j)
            grad(j) = 0;
            continue;
        end
        Jp = objective_x2(x2p, x1, M22, w2, w1, h2);
        Jm = objective_x2(x2m, x1, M22, w2, w1, h2);
        grad(j) = (Jp - Jm) / (x2p(j) - x2m(j));
    end
end

function H = hess_fd_x1(x1, x2, M11, w1, w2, h1, gradEps, minSlope)
    N1 = numel(x1);
    I = 2:N1-1;
    n = numel(I);
    H = zeros(n, n);
    for k = 1:n
        idx = I(k);
        x1p = perturb_monotone(x1, idx, gradEps, minSlope);
        x1m = perturb_monotone(x1, idx, -gradEps, minSlope);
        if x1p(idx) == x1m(idx)
            continue;
        end
        gp = grad_fd_x1(x1p, x2, M11, w1, w2, h1, gradEps, minSlope);
        gm = grad_fd_x1(x1m, x2, M11, w1, w2, h1, gradEps, minSlope);
        H(:, k) = (gp(I) - gm(I)) / (x1p(idx) - x1m(idx));
    end
end

function H = hess_fd_x2(x2, x1, M22, w2, w1, h2, gradEps, minSlope)
    N2 = numel(x2);
    I = 2:N2-1;
    n = numel(I);
    H = zeros(n, n);
    for k = 1:n
        idx = I(k);
        x2p = perturb_monotone(x2, idx, gradEps, minSlope);
        x2m = perturb_monotone(x2, idx, -gradEps, minSlope);
        if x2p(idx) == x2m(idx)
            continue;
        end
        gp = grad_fd_x2(x2p, x1, M22, w2, w1, h2, gradEps, minSlope);
        gm = grad_fd_x2(x2m, x1, M22, w2, w1, h2, gradEps, minSlope);
        H(:, k) = (gp(I) - gm(I)) / (x2p(idx) - x2m(idx));
    end
end

function x_new = perturb_monotone(x, idx, delta, minSlope)
    x_new = x;
    x_new(idx) = x_new(idx) + delta;
    if x_new(idx) <= x_new(idx-1) + minSlope
        x_new(idx) = x_new(idx-1) + minSlope;
    end
    if x_new(idx) >= x_new(idx+1) - minSlope
        x_new(idx) = x_new(idx+1) - minSlope;
    end
end

function p = d1_oneSidedCentral(x,h)
    N = numel(x);
    p = zeros(N,1);
    p(1)     = (x(2)-x(1))/h;
    p(N)     = (x(N)-x(N-1))/h;
    p(2:N-1) = (x(3:N)-x(1:N-2))/(2*h);
end

function J = fullObjective(x1,x2,M11,M22,Delta1,Delta2,h1,h2)
    N1 = numel(x1); N2 = numel(x2);
    p = d1_oneSidedCentral(x1,h1);
    q = d1_oneSidedCentral(x2,h2);

    w1 = h1*ones(N1,1); w1(1)=h1/2; w1(end)=h1/2;
    w2 = h2*ones(N2,1); w2(1)=h2/2; w2(end)=h2/2;

    J = 0;
    for i = 1:N1
        for j = 1:N2
            r1 = Delta1^2 * M11(x1(i),x2(j)) * p(i)^2 - 1;
            r2 = Delta2^2 * M22(x1(i),x2(j)) * q(j)^2 - 1;
            J = J + w1(i)*w2(j) * (r1^2 + r2^2);
        end
    end
end

function Delta1 = updateDelta1(x1,x2,M11,w1,w2,h1)
    N1 = numel(x1); N2 = numel(x2);
    p = d1_oneSidedCentral(x1,h1);

    num = 0; den = 0;
    for i = 1:N1
        for j = 1:N2
            a = M11(x1(i),x2(j)) * p(i)^2;
            wij = w1(i)*w2(j);
            num = num + wij*a;
            den = den + wij*a^2;
        end
    end
    Delta1 = sqrt(max(num / max(den,eps), eps));
end

function Delta2 = updateDelta2(x1,x2,M22,w1,w2,h2)
    N1 = numel(x1); N2 = numel(x2);
    q = d1_oneSidedCentral(x2,h2);

    num = 0; den = 0;
    for i = 1:N1
        for j = 1:N2
            b = M22(x1(i),x2(j)) * q(j)^2;
            wij = w1(i)*w2(j);
            num = num + wij*b;
            den = den + wij*b^2;
        end
    end
    Delta2 = sqrt(max(num / max(den,eps), eps));
end

function [A1, b1] = buildAb_x1(x1,x2,M11,Delta1,h1,w2,fdEpsX)
    % Build A1(i), b1(i) for i=2..N1-1 (interior unknowns)
    N1 = numel(x1); N2 = numel(x2);

    p = d1_oneSidedCentral(x1,h1);

    A1 = zeros(N1,1);
    b1 = zeros(N1,1);

    for i = 2:N1-1
        xi = x1(i);
        pi = p(i);

        Ai = 0;
        bi = 0;

        for j = 1:N2
            x2j = x2(j);
            m11 = M11(xi, x2j);

            % dM11/dx1 by physical-space central FD
            m11p = M11(xi + fdEpsX, x2j);
            m11m = M11(xi - fdEpsX, x2j);
            dm11dx1 = (m11p - m11m) / (2*fdEpsX);

            r = Delta1^2 * m11 * (pi^2) - 1;

            Ai = Ai + w2(j) * (m11 * r);
            bi = bi + w2(j) * (r * (pi^2) * dm11dx1);
        end

        A1(i) = 2*Ai;
        b1(i) = bi;
    end
end

function [A2, b2] = buildAb_x2(x1,x2,M22,Delta2,h2,w1,fdEpsX)
    % Build A2(j), b2(j) for j=2..N2-1 (interior unknowns)
    N1 = numel(x1); N2 = numel(x2);

    q = d1_oneSidedCentral(x2,h2);

    A2 = zeros(N2,1);
    b2 = zeros(N2,1);

    for j = 2:N2-1
        xj = x2(j);
        qj = q(j);

        Aj = 0;
        bj = 0;

        for i = 1:N1
            x1i = x1(i);
            m22 = M22(x1i, xj);

            % dM22/dx2 by physical-space central FD
            m22p = M22(x1i, xj + fdEpsX);
            m22m = M22(x1i, xj - fdEpsX);
            dm22dx2 = (m22p - m22m) / (2*fdEpsX);

            r = Delta2^2 * m22 * (qj^2) - 1;

            Aj = Aj + w1(i) * (m22 * r);
            bj = bj + w1(i) * (r * (qj^2) * dm22dx2);
        end

        A2(j) = 2*Aj;
        b2(j) = bj;
    end
end

function x = solve_varcoef_poisson_1d(A, b, h, xL, xR)
    % Solve for x on 1D grid with Dirichlet x(1)=xL, x(N)=xR:
    %   A(i) * (x(i+1)-2x(i)+x(i-1))/h^2 = b(i),  i=2..N-1
    %
    % A, b are length N; only interior used.
    N = numel(A);
    I = 2:N-1;
    n = numel(I);

    Ai = A(I);
    bi = b(I);

    % Safeguard for small/negative coefficients:
    % if A flips sign, the elliptic-like operator changes type; we clamp magnitude.
    Ai = sign(Ai) .* max(abs(Ai), 1e-12);

    main = (-2*Ai)/(h^2);
    off  = ( Ai)/(h^2);

    % Tridiagonal matrix
    T = spdiags([off, main, off], -1:1, n, n);

    % RHS includes boundary contributions
    rhs = bi;
    rhs(1)   = rhs(1)   - off(1)   * xL; % i=2 depends on x(1)
    rhs(end) = rhs(end) - off(end) * xR; % i=N-1 depends on x(N)

    x_int = T \ rhs;

    x = zeros(N,1);
    x(1) = xL; x(end) = xR;
    x(I) = x_int;
end
