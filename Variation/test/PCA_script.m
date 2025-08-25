clear
Nx1 = 40;
Nx2 = 40;
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);
[x1, x2] = meshgrid(s1, s2);
m = GetM(x1, x2, 1, 0.1);

[x1_temp, x2_temp] = meshgrid(s1, s2);
x1 = x1_temp*cos(pi/4) - x2_temp*sin(pi/4);
x2 = x1_temp*sin(pi/4) + x2_temp*cos(pi/4);
x1 = x1 + 0.001*randn(size(x1));
x2 = x2 + 0.001*randn(size(x2));
m = GetM(x1, x2, 6, 0.1);

M11 = m.M11;
M12 = m.M12;
M22 = m.M22;

M1 = cat(3, M11, M12);
M2 = cat(3, M12, M22);
M = cat(4, M1, M2);
[eigvals, eigvecs] = Accurate_eig(M);

anisotropy = max(eigvals, [], 3) ./ min(eigvals, [], 3);

figure()
f = pcolor(x1, x2, anisotropy);
set(f, 'EdgeColor', 'none'); colorbar;
title('anisotropy, log(1+abs(Mii)) scaled');

figure()
space = 7;
quiver(x1(1:space:end,1:space:end), x2(1:space:end,1:space:end),...
       eigvals(1:space:end,1:space:end,1).*eigvecs(1:space:end,1:space:end,1,1), eigvals(1:space:end,1:space:end,1).*eigvecs(1:space:end,1:space:end,2,1), 'AutoScale', 'on', Color='k'); hold on

quiver(x1(1:space:end,1:space:end), x2(1:space:end,1:space:end),...
       eigvals(1:space:end,1:space:end,2).*eigvecs(1:space:end,1:space:end,1,2), eigvals(1:space:end,1:space:end,2).*eigvecs(1:space:end,1:space:end,2,2), 'AutoScale', 'on', Color='k');

%quiver(x1(1:space:end,1:space:end), x2(1:space:end,1:space:end),...
%       eigvecs(1:space:end,1:space:end,1,1), eigvecs(1:space:end,1:space:end,2,1), 'AutoScale', 'on'); hold on
%quiver(x1(1:space:end,1:space:end), x2(1:space:end,1:space:end),...
%       eigvecs(1:space:end,1:space:end,1,2), eigvecs(1:space:end,1:space:end,2,2), 'AutoScale', 'on');

xlim([-1.1,1.1]); ylim([-0.1,1.1]); 
grid on
axis equal

%%
function [eigvals, eigvecs] = Accurate_eig(M)
    eigvals = zeros(size(M,1), size(M,2), 2);
    eigvecs = zeros(size(M,1), size(M,2), 2, 2);
    for i=1:size(M,1)
        for j=1:size(M,2)
            [V,D] = eig(reshape(M(i,j,:,:),2,2));
            eigvals(i,j,:) = diag(D);
            eigvecs(i,j,:,:) = V;
        end
    end
end



function [eigvals, eigvecs] = eig2x2_metric(M)
% eig2x2_metric computes eigenvalues and eigenvectors for an array of 2x2 matrices.
%
%   [eigvals, eigvecs] = eig2x2_metric(M)
%
%   Input:
%       M - a N1 x N2 x 2 x 2 array, where each 2x2 slice is a symmetric metric tensor.
%
%   Output:
%       eigvals - a N1 x N2 x 2 array of eigenvalues. For each grid point (i,j),
%                 eigvals(i,j,1) is the larger eigenvalue and eigvals(i,j,2) is the smaller one.
%       eigvecs - a N1 x N2 x 2 x 2 array of eigenvectors. For each grid point (i,j):
%                 eigvecs(i,j,:,1) is the unit eigenvector corresponding to eigvals(i,j,1),
%                 eigvecs(i,j,:,2) is the unit eigenvector corresponding to eigvals(i,j,2).
%
%   Note: The function assumes that M(:,:,1,2) equals M(:,:,2,1) (i.e., each matrix is symmetric).

% Extract matrix elements
a = M(:,:,1,1);
b = M(:,:,1,2);  % Note: Since the matrix is symmetric, M(:,:,2,1) == b.
c = M(:,:,2,2);

% Compute eigenvalues for symmetric 2x2 matrix:
%   lambda = (a+c ± sqrt((a-c)^2 + 4*b^2)) / 2
traceM = a + c;
D = sqrt((a - c).^2 + 4*b.^2);
lambda1 = (traceM + D) / 2;
lambda2 = (traceM - D) / 2;

diag_mask = b < eps;
lambda1(diag_mask) = a(diag_mask);
lambda2(diag_mask) = c(diag_mask);

% Arrange eigenvalues into a N1 x N2 x 2 array
eigvals = cat(3, lambda1, lambda2);

% Compute eigenvectors.
% For eigenvalue lambda, an eigenvector is given by [b; lambda - a].
v1_1 = b;
v1_2 = lambda1 - a;
v2_1 = b;
v2_2 = lambda2 - a;

v1_1(diag_mask) = 1;
v1_2(diag_mask) = 0;
v2_1(diag_mask) = 0;
v2_2(diag_mask) = 1;


% Compute norms for normalization
norm1 = sqrt(v1_1.^2 + v1_2.^2);
norm2 = sqrt(v2_1.^2 + v2_2.^2);

% Handle any cases where the norm is zero (degenerate case)
zero_mask1 = norm1 < eps;
zero_mask2 = norm2 < eps;
v1_1(zero_mask1) = 1;  % default eigenvector [1; 0]
v1_2(zero_mask1) = 0;
v2_1(zero_mask2) = 1;
v2_2(zero_mask2) = 0;
norm1(zero_mask1) = 1;
norm2(zero_mask2) = 1;

% Normalize eigenvectors
v1_1 = v1_1 ./ norm1;
v1_2 = v1_2 ./ norm1;
v2_1 = v2_1 ./ norm2;
v2_2 = v2_2 ./ norm2;

% Construct eigenvector array:
% For each grid point (i,j), eigenvectors(:,:, :, 1) is the eigenvector corresponding to lambda1,
% and eigenvectors(:,:, :, 2) is the eigenvector corresponding to lambda2.
eigvecs = zeros(size(M));
eigvecs(:,:,1,1) = v1_1;
eigvecs(:,:,2,1) = v1_2;
eigvecs(:,:,1,2) = v2_1;
eigvecs(:,:,2,2) = v2_2;
end
