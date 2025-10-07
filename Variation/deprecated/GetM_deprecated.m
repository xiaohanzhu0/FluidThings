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
