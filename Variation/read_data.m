clear
addpath('../data');

m = load('metric_airfoil_bodyfitted.mat');
m_component = ["M11"; "M12"; "M22"];

%%
space1 = 10;
space2 = 2;
Nx1 = size(m.x_metric, 1);
Nx2 = size(m.x_metric, 2);

for i = 1:space1:Nx1
    scatter(m.x_metric(i,1:space2:end), m.y_metric(i,1:space2:end),'.','k'); hold on
end

%% Metric magnitude distribution
figure()
for i=1:3
    subplot(3,1,i);
    histogram(m.metric(:,:,i));
    set(gca,'YScale','log'); grid on;
    title(m_component(i));
end

%% 
for i=1:3
    figure()
    f = pcolor(m.x_metric, m.y_metric, log(1+abs(m.metric(:,:,i))));
    set(f, 'EdgeColor', 'none'); colorbar;
    title(strcat(m_component(i),', log(1+abs(Mii)) scaled'));
end

%%
M11 = m.metric(:,:,1);
M12 = m.metric(:,:,2);
M22 = m.metric(:,:,3);

M1 = cat(3, M11, M12);
M2 = cat(3, M12, M22);
M = cat(4, M1, M2);
clear M1 M2 M11 M12 M22

%%
[eigvals, eigvecs] = eig2x2_metric(M);
anisotropy = max(eigvals, [], 3) ./ min(eigvals, [], 3);

f = pcolor(m.x_metric, m.y_metric, log(1+abs(anisotropy)));
set(f, 'EdgeColor', 'none'); colorbar;
title('anisotropy, log(1+abs(Mii)) scaled');
%%
space = 10;
quiver(m.x_metric(1:space:end,:), m.y_metric(1:space:end,:),...
       eigvecs(1:space:end,:,1,1), eigvecs(1:space:end,:,2,1), 'AutoScale', 'on'); hold on
quiver(m.x_metric(1:space:end,:), m.y_metric(1:space:end,:),...
       eigvecs(1:space:end,:,1,2), eigvecs(1:space:end,:,2,2), 'AutoScale', 'on');

%% Look at gradient
[M11grad1, M11grad2] = DCentralUneven(m.metric(:,:,1), m.metric(:,:,1), m.y_metric, m.x_metric);

figure()
f = pcolor(m.x_metric, m.y_metric, log(1+abs(M11grad1)));
set(f, 'EdgeColor', 'none'); colorbar;
title(strcat(m_component(1),', log(1+abs(Mii)) scaled'));

%%
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

% Arrange eigenvalues into a N1 x N2 x 2 array
eigvals = cat(3, lambda1, lambda2);

% Compute eigenvectors.
% For eigenvalue lambda, an eigenvector is given by [b; lambda - a].
v1_1 = b;
v1_2 = lambda1 - a;
v2_1 = b;
v2_2 = lambda2 - a;

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
