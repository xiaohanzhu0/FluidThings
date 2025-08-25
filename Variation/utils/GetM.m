% If M_type is integer, then M is supposed to be symbolically defined
% If M_type is struct {x1_samp, x2_samp, M_samp, x1, x2}, then M would be
% numerically defined

function M = GetM(x1, x2, M_type, C)
    if isnumeric(M_type)
        if M_type == 1
            M.M11 = 40000 * (1 + C*15 * x1).^(-2);
            M.M22 = 40000 * (1 + C*15 * x2).^(-2);
            M.M12 = zeros(size(x1));

            M.dM11dx1 = -C*1200000 * (1 + C*15 * x1).^(-3);
            M.dM11dx2 = zeros(size(x1));
            M.dM22dx1 = zeros(size(x1));
            M.dM22dx2 = -C*1200000 * (1 + C*15 * x2).^(-3);
            M.dM12dx1 = zeros(size(x1));
            M.dM12dx2 = zeros(size(x1));
        elseif M_type == 2
            freq = 1;
            M.M11 = 1000 + C*600*sin(freq*2*pi*x1).*sin(freq*2*pi*x2);
            M.M12 = zeros(size(x1));
            M.M22 = 1000 - C*600*sin(freq*2*pi*x1).*sin(freq*2*pi*x2);
            M.dM11dx1 = freq*C*1200*pi*cos(freq*2*pi*x1).*sin(freq*2*pi*x2);
            M.dM11dx2 = freq*C*1200*pi*sin(freq*2*pi*x1).*cos(freq*2*pi*x2);
            M.dM12dx1 = zeros(size(x1));
            M.dM12dx2 = zeros(size(x1));
            M.dM22dx1 = -freq*C*1200*pi*cos(freq*2*pi*x1).*sin(freq*2*pi*x2);
            M.dM22dx2 = -freq*C*1200*pi*sin(freq*2*pi*x1).*cos(freq*2*pi*x2);
        elseif M_type == 3
            [Nx2, Nx1] = size(x1);
            M.M11 = 2000*ones(Nx2, Nx1);
            M.M22 = 2000*ones(Nx2, Nx1);
            M.M12 = zeros(Nx2, Nx1);
            [Nx2, Nx1] = size(x1);
            M.dM11dx1 = zeros(Nx2, Nx1);
            M.dM11dx2 = zeros(Nx2, Nx1);
            M.dM22dx1 = zeros(Nx2, Nx1);
            M.dM22dx2 = zeros(Nx2, Nx1);
            M.dM12dx1 = zeros(Nx2, Nx1);
            M.dM12dx2 = zeros(Nx2, Nx1);
        elseif M_type == 4 || M_type == 8
            [Nx2, Nx1] = size(x1);
            M.M11 = 400*ones(Nx2, Nx1);
            M.M22 = 400*ones(Nx2, Nx1);
            M.M12 = zeros(Nx2, Nx1);
            [Nx2, Nx1] = size(x1);
            M.dM11dx1 = zeros(Nx2, Nx1);
            M.dM11dx2 = zeros(Nx2, Nx1);
            M.dM22dx1 = zeros(Nx2, Nx1);
            M.dM22dx2 = zeros(Nx2, Nx1);
            M.dM12dx1 = zeros(Nx2, Nx1);
            M.dM12dx2 = zeros(Nx2, Nx1);
        elseif M_type == 5
            M.M11 = 1+x1;
            M.M22 = 1+x2;
            M.M12 = 1.3*x2;

            M.dM11dx1 = ones(size(x1));
            M.dM11dx2 = zeros(size(x1));
            M.dM22dx1 = zeros(size(x1));
            M.dM22dx2 = ones(size(x1));
            M.dM12dx1 = zeros(size(x1));
            M.dM12dx2 = 1.3*ones(size(x1));
        elseif M_type == 6
            theta = pi/6;

            s1 = x1*cos(theta) + x2*sin(theta);
            s2 = -x1*sin(theta) + x2*cos(theta);

            M11_temp = 40000 * (1 + C*15 * s1).^(-2);
            M22_temp = 40000 * (1 + C*15 * s2).^(-2);
            M12_temp = zeros(size(s1));

            M.M11 = cos(theta)*(M11_temp*cos(theta)-M12_temp*sin(theta)) - sin(theta)*(M12_temp*cos(theta)-M22_temp*sin(theta));
            M.M22 = sin(theta)*(M11_temp*sin(theta)+M12_temp*cos(theta)) + cos(theta)*(M12_temp*sin(theta)+M22_temp*cos(theta));
            M.M12 = cos(theta)*(M11_temp*sin(theta)+M12_temp*cos(theta)) - sin(theta)*(M12_temp*sin(theta)+M22_temp*cos(theta));

            [M.dM11dx1, M.dM11dx2] = metric_grad(M.M11, x1, x2);
            [M.dM22dx1, M.dM22dx2] = metric_grad(M.M22, x1, x2);
            [M.dM12dx1, M.dM12dx2] = metric_grad(M.M12, x1, x2);
        elseif M_type == 7
            s1 = x1 - x2/2;
            s2 = x2/2;
            
            M11_temp = 40000 * (1 + C*15 * s1).^(-2);
            M22_temp = 40000 * (1 + C*15 * s2).^(-2);
            M12_temp = zeros(size(s1));

            M.M11 = M11_temp;
            M.M12 = M11_temp + 2*M12_temp;
            M.M22 = M11_temp + 4*(M12_temp+M22_temp);

            [M.dM11dx1, M.dM11dx2] = metric_grad(M.M11, x1, x2);
            [M.dM22dx1, M.dM22dx2] = metric_grad(M.M22, x1, x2);
            [M.dM12dx1, M.dM12dx2] = metric_grad(M.M12, x1, x2);
        end
    else
        x1_samp = M_type.x_metric;
        x2_samp = M_type.y_metric;
        M_samp = M_type.metric;
        %tic
        %M.M11 = griddata(x1_samp',x2_samp',M_samp(:,:,1)',x1,x2,"linear");
        %M.M12 = griddata(x1_samp',x2_samp',M_samp(:,:,2)',x1,x2,"linear");
        %M.M22 = griddata(x1_samp',x2_samp',M_samp(:,:,3)',x1,x2,"linear");
        %toc

        [Nx2, Nx1] = size(x1);
        M.M11 = M_type.F11(x1(:),x2(:)); M.M11 = reshape(M.M11, Nx2, Nx1);
        M.M12 = M_type.F12(x1(:),x2(:)); M.M12 = reshape(M.M12, Nx2, Nx1);
        M.M22 = M_type.F22(x1(:),x2(:)); M.M22 = reshape(M.M22, Nx2, Nx1);


        [M.dM11dx1, M.dM11dx2] = metric_grad(M.M11, x1, x2);
        [M.dM22dx1, M.dM22dx2] = metric_grad(M.M22, x1, x2);
        [M.dM12dx1, M.dM12dx2] = metric_grad(M.M12, x1, x2);

        %space = 10;
        %Metric1 = cat(3, M.M11, M.M12);
        %Metric2 = cat(3, M.M12, M.M22);
        %Metric = cat(4, Metric1, Metric2);
     
        %[eigvals, eigvecs] = eig2x2_metric(Metric);
        %pcolor(x1,x2,1./sqrt(eigvals(:,:,1).*eigvals(:,:,2))); hold on
        %figure;
        
        %drawEigenEllipses(eigvecs(1:space:end,:,:,:), 1e-1./sqrt(eigvals(1:space:end,:,:)), x1(1:space:end,:), x2(1:space:end,:)); hold on
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