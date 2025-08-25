function [dMdx1, dMdx2] = metric_grad(M, x1, x2)
%METRIC_GRADIENT  Compute ∂M/∂x₁ and ∂M/∂x₂ on a non‐uniform, skewed grid
%
%   [dMdx1, dMdx2] = metric_gradient(M, x1, x2)
%
%   Inputs:
%     M  – a 2D array of size (nrows × ncols), representing one component of the metric tensor
%     x1 – an (nrows × ncols) array of the physical x₁‐coordinates at each node
%     x2 – an (nrows × ncols) array of the physical x₂‐coordinates at each node
%
%   Outputs:
%     dMdx1 – an (nrows × ncols) array containing ∂M/∂x₁ at interior nodes (zeros on the 4‐point boundary)
%     dMdx2 – an (nrows × ncols) array containing ∂M/∂x₂ at interior nodes (zeros on the 4‐point boundary)
%
%   This function uses a 5‐point, chain‐rule stencil.  At each interior index (i,j):
%     – Δξx₁ = x₁(i+1,j) − x₁(i−1,j)
%     – Δξx₂ = x₂(i+1,j) − x₂(i−1,j)
%     – Δηx₁ = x₁(i,j+1) − x₁(i,j−1)
%     – Δηx₂ = x₂(i,j+1) − x₂(i,j−1)
%     – ΔξM  = M(i+1,j)   − M(i−1,j)
%     – ΔηM  = M(i,j+1)   − M(i,j−1)
%
%   Then one solves the 2×2 linear system:
%       [ Δξx₁   Δξx₂ ] [∂M/∂x₁]   = [ ΔξM ]
%       [ Δηx₁   Δηx₂ ] [∂M/∂x₂]     [ ΔηM ]
%
%   This is vectorized so that all interior points are handled in one shot.

[nrows, ncols] = size(M);

% Pre‐allocate output arrays (zeros on boundary)
dMdx1 = zeros(nrows, ncols);
dMdx2 = zeros(nrows, ncols);

% Define index ranges for “interior” (skip the 1‐pixel boundary)
%   i = 2:(nrows−1),  j = 2:(ncols−1)
%
% To vectorize, we form sub‐blocks of size (nrows−2)×(ncols−2):
%   x1(i+1,j) corresponds to x1(3:end,   2:end−1)
%   x1(i−1,j) corresponds to x1(1:end−2, 2:end−1)
%   x1(i,j+1) corresponds to x1(2:end−1, 3:end)
%   x1(i,j−1) corresponds to x1(2:end−1, 1:end−2)
%   etc.

% 1) Compute Δξx₁, Δξx₂, Δηx₁, Δηx₂ over all interior points at once
deta_x1  = x1(3:end,   2:end-1) - x1(1:end-2, 2:end-1);   % (i+1,j) − (i−1,j)
deta_x2  = x2(3:end,   2:end-1) - x2(1:end-2, 2:end-1);

dxi_x1 = x1(2:end-1, 3:end  ) - x1(2:end-1, 1:end-2);   % (i,j+1) − (i,j−1)
dxi_x2 = x2(2:end-1, 3:end  ) - x2(2:end-1, 1:end-2);

% 2) Compute ΔξM and ΔηM over the same interior stencil
deta_M  = M(3:end,   2:end-1) - M(1:end-2, 2:end-1);   % M(i+1,j) − M(i−1,j)
dxi_M = M(2:end-1, 3:end  ) - M(2:end-1, 1:end-2);   % M(i,j+1) − M(i,j−1)

% 3) Build the determinant of the 2×2 Jacobian for each interior point
%       detJ = Δξx₁ * Δηx₂ − Δξx₂ * Δηx₁
detJ = dxi_x1 .* deta_x2 - dxi_x2 .* deta_x1;

% 4) Compute each component of the inverse of J = [Δξx₁ Δξx₂; Δηx₁ Δηx₂]
%    inv(J) = (1/detJ) * [ Δηx₂  −Δξx₂;  −Δηx₁  Δξx₁ ]
inv11 =  deta_x2 ./ detJ;    % coefficient mapping ΔξM → ∂M/∂x₁
inv12 = -dxi_x2 ./ detJ;     % coefficient mapping ΔηM → ∂M/∂x₁
inv21 = -deta_x1 ./ detJ;    % coefficient mapping ΔξM → ∂M/∂x₂
inv22 =  dxi_x1 ./ detJ;     % coefficient mapping ΔηM → ∂M/∂x₂

% 5) Compute ∂M/∂x₁ and ∂M/∂x₂ on the interior block
dMdx1_int = inv11 .* dxi_M + inv12 .* deta_M;    % (nrows−2)×(ncols−2)
dMdx2_int = inv21 .* dxi_M + inv22 .* deta_M;

% 6) Scatter back into the full‐size arrays (leaving zeros on the 1‐pixel boundary)
dMdx1(2:end-1, 2:end-1) = dMdx1_int;
dMdx2(2:end-1, 2:end-1) = dMdx2_int;

end
