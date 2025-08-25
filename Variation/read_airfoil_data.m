
[data, vartype, varname, ncvs, gridfile] = readTurtleFields('metricField.fields');
[x_cv, x_no, x_fa, ncvs, ...
          bndNames, bnd_block_begin_end, ...
          int_block_begin_end_A, int_block_begin_end_B, int_offset, int_angle] = readTurtleGrid('airfoil_18M_coarseIJK.grid');

aux = data{1,1}{1,6}(:,:,:,1:2,1:2);
dim = size(aux);
M = reshape(aux(:,:,1,:,:), dim(1), dim(2), dim(4), dim(5));

%%
legend()
for i = 1:16
    figure(1)
    x_metric = x_cv{1,i}(:,:,1,1);
    y_metric = x_cv{1,i}(:,:,1,2);
    plot(x_metric(:,1), y_metric(:,1),'k',LineWidth=1.5); hold on
    plot(x_metric(:,end), y_metric(:,end),'k',LineWidth=1.5); hold on
    plot(x_metric(1,:), y_metric(1,:),'k',LineWidth=1.5); hold on
    plot(x_metric(end,:), y_metric(end,:),'k',LineWidth=1.5); hold on
    for j=1:size(x_metric, 1)
        plot(x_metric(j,:), y_metric(j,:),'k'); hold on
    end
    for j=1:size(x_metric, 2)
        plot(x_metric(:,j), y_metric(:,j),'k'); hold on
    end

    text(mean(x_metric,'all'),mean(y_metric,'all'),num2str(i),Color='r'); hold on

    figure(2);
    aux = data{1,1}{1,i}(:,:,:,1:2,1:2);
    dim = size(aux);
    M = reshape(aux(:,:,1,:,:), dim(1), dim(2), dim(4), dim(5));
    pcolor(x_metric, y_metric, M(:,:,1,1)); hold on
end

%%
x_all = [x_cv{1,9}(:,:,1,1), x_cv{1,8}(:,:,1,1);
            x_cv{1,10}(:,:,1,1),x_cv{1,1}(:,:,1,1);
            x_cv{1,11}(:,:,1,1),x_cv{1,2}(:,:,1,1);
            x_cv{1,12}(:,:,1,1),x_cv{1,3}(:,:,1,1);
            x_cv{1,13}(:,:,1,1),x_cv{1,4}(:,:,1,1);
            x_cv{1,14}(:,:,1,1),x_cv{1,5}(:,:,1,1);
            x_cv{1,15}(:,:,1,1),x_cv{1,6}(:,:,1,1);
            x_cv{1,16}(:,:,1,1),x_cv{1,7}(:,:,1,1)];

y_all = [x_cv{1,9}(:,:,1,2), x_cv{1,8}(:,:,1,2);
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

pcolor(x_all, y_all, M(:,:,1,1));
%%
xi = 1:190;
eta = 1:10;
eigvals = zeros(length(xi),length(eta),2);
eigvecs = zeros(length(xi),length(eta),2,2);

for i=1:length(xi)
    for j=1:length(eta)
        [v,d] = eig(reshape(M(i,j,:,:), 2, 2));
        d = diag(d);
        eigvals(i,j,:) = d;
        eigvecs(i,j,:,:) = v;
    end
end

drawEigenEllipses(eigvecs, 5e-2./sqrt(eigvals), x_metric(xi,eta), y_metric(xi,eta)); hold on
plot(x_metric(xi,1), y_metric(xi,1),'r')


%%
scatter3(data{1, 5}{1, 5}(:,:,1), data{1, 6}{1, 5}(:,:,1), data{1, 7}{1, 5}(:,:,1),'.'); hold on
%%
clear
addpath('~/Files/data/Mesh_generation/Airfoil');

m = load('metric_airfoil_bodyfitted.mat');
m_component = ["M11"; "M12"; "M22"];


%%
ind = 100;
ds = sqrt((m.x_metric(3:end,ind) - m.x_metric(1:end-2,ind)).^2 + (m.y_metric(3:end,ind) - m.y_metric(1:end-2,ind)).^2);
displacement = cumsum(ds);
subplot(2,1,1);
plot(displacement, m.metric(2:end-1,ind,1));
title("physical space"); xlabel('x'); ylabel('M');
grid on
subplot(2,1,2);
plot(linspace(0,1,size(m.metric(2:end-1,ind,1),1)),m.metric(2:end-1,ind,1));
title("computational space"); xlabel('\xi'); ylabel('M');
grid on
%%
h = gaussian1d(100, 51);
filtered = conv(m.metric(:,30,1),h,'same');
plot(linspace(0,1,size(m.metric(:,ind,1),1)),m.metric(:,30,1)); hold on
plot(linspace(0,1,size(m.metric(:,ind,1),1)),filtered); hold on
%%
h = gaussian1d(100, 100);
filtered = conv2(h,h,m.metric(:,:,1),'same');
plot(linspace(0,1,size(m.metric(:,ind,1),1)),filtered(:,1)*200,".");hold on

filtered = conv2(h,h,filtered,'same');
plot(linspace(0,1,size(m.metric(:,ind,1),1)),filtered(:,1)*200,".");hold on

filtered = conv2(h,h,filtered,'same');
plot(linspace(0,1,size(m.metric(:,ind,1),1)),filtered(:,1)*200,".");hold on

filtered = conv2(h,h,filtered,'same');
plot(linspace(0,1,size(m.metric(:,ind,1),1)),filtered(:,1)*200,".");hold on

%%
for i=1:138
    plot(linspace(0,1,size(m.metric(:,ind,1),1)),filtered(:,i)*51,".");
    pause(0.1)
end
%%
plot(displacement, m.metric(2:end-1,ind,1)); hold on
temp = nufft(m.metric(2:end-1,ind,1), displacement);
temp(abs(temp)<1000) = 0;
plot(displacement, abs(ifft(temp)));

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
    f = surf(m.x_metric, m.y_metric, m.metric(:,:,i));
    set(f, 'EdgeColor', 'none'); v;
    title(strcat(m_component(i),', log(1+abs(Mii)) scaled'));
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

%% Expensive but accurate way of computing eigensystem
xi = 149:2:2362;
eta = 1:5;
eigvals = zeros(length(xi),length(eta),2);
eigvecs = zeros(length(xi),length(eta),2,2);
for i=1:length(xi)
    for j=1:length(eta)
        [v,d] = eig(reshape(M(i,j,:,:), 2, 2));
        d = diag(d);
        eigvals(i,j,:) = d;
        eigvecs(i,j,:,:) = v;
    end
end

drawEigenEllipses(eigvecs, 5e-2./sqrt(eigvals), m.x_metric(xi,eta), m.y_metric(xi,eta)); hold on
plot(m.x_metric(xi,1), m.y_metric(xi,1),'r')
scatter(m.x_metric(111,35), m.y_metric(111,35))
%%
[eigvals, eigvecs] = eig2x2_metric(M);
anisotropy = max(eigvals, [], 3) ./ min(eigvals, [], 3);

f = pcolor(m.x_metric, m.y_metric, abs(anisotropy));
%set(f, 'EdgeColor', 'none');
colorbar;
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
%2361
%drawEigenEllipses(eigvecs(149:150,1,:,:), 1./eigvals(149:150,1,:), m.x_metric(149:150,1), m.y_metric(149:150,1)); hold on
%drawEigenEllipses(eigvecs(2361:2362,1,:,:), 1./eigvals(2361:2362,1,:), m.x_metric(2361:2362,1), m.y_metric(2361:2362,1))
space = 5;

drawEigenEllipses(eigvecs(1:space:end,1,:,:), 1e-5.*sqrt(eigvals(1:space:end,1,:)), m.x_metric(1:space:end,1), m.y_metric(1:space:end,1)); hold on
%%
drawEigenEllipses(d, v, m.x_metric(1988,1), m.y_metric(1988,1)); hold on
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

function h = gaussian1d(sigma, ksize)
%GAUSSIAN1D  Return a normalized 1-D Gaussian kernel.
%   h = GAUSSIAN1D(sigma, ksize) returns a row vector of length ksize
%   containing samples of a Gaussian with standard deviation sigma.
%
%   Example:
%     h = gaussian1d(2, 11);
%     plot(h, '-o'); axis tight; title('1D Gaussian Kernel');

    % ensure kernel size is a positive integer
    validateattributes(ksize, {'numeric'}, {'scalar','integer','positive'});
    validateattributes(sigma, {'numeric'}, {'scalar','positive'});

    % define sample positions centered at zero
    half = (ksize-1)/2;
    if mod(ksize,2)==0
        % even length—center between the two middle samples
        x = -half:(half-1);
    else
        % odd length—symmetric around zero
        x = -floor(half):floor(half);
    end

    % compute the Gaussian
    h = exp( - (x.^2) / (2*sigma^2) );

    % normalize so sum(h) == 1
    h = h / sum(h);
end

function drawEigenEllipses(eigVecs, eigVals, X, Y, nPts, lineSpec)
% drawEigenEllipses  Draws ellipses from local 2×2 eigensystems.
%
%   drawEigenEllipses(eigVecs, eigVals, X, Y)
%   plots ellipses at each (X(i,j),Y(i,j)) with principal axes
%   given by the columns of eigVecs(i,j,:,:) and lengths
%   eigVals(i,j,1), eigVals(i,j,2).
%
%   drawEigenEllipses(..., nPts) specifies the number of points per ellipse
%   (default 50).
%
%   drawEigenEllipses(..., lineSpec) passes a MATLAB plot spec (e.g. 'r-')
%   to control color/linestyle.

if nargin < 5 || isempty(nPts),   nPts   = 50;  end
if nargin < 6 || isempty(lineSpec), lineSpec = 'b-'; end

% Precompute the unit circle
theta = linspace(0,2*pi,nPts);
unitCircle = [cos(theta); sin(theta)];  % 2×nPts

hold on;
[M,N,~,~] = size(eigVecs);
for i = 1:M
    for j = 1:N
        % extract local 2×2 eigenvector matrix
        U = squeeze(eigVecs(i,j,:,:));      % 2×2
        % build scaling matrix from eigenvalues
        L = diag(squeeze(eigVals(i,j,:)));  % 2×2
        % transform unit circle into ellipse
        ellipsePts = U * L * unitCircle;    % 2×nPts
        % translate to center (X(i,j),Y(i,j))
        xE = ellipsePts(1,:) + X(i,j);
        yE = ellipsePts(2,:) + Y(i,j);
        % plot
        plot(xE, yE, lineSpec);
    end
end
axis equal
hold off;
end