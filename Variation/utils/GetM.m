% If M_type is integer, then M is supposed to be symbolically defined
% If M_type is struct {x1_samp, x2_samp, M_samp, x1, x2}, then M would be
% numerically defined

function M = GetM(x1, x2, M_type)
    if isnumeric(M_type)
        if M_type == 1
            M.M11 = 40000 * (1 + 15 * x1).^(-2);
            M.M22 = 40000 * (1 + 15 * x2).^(-2);
            M.dM11dx1 = -1200000 * (1 + 15 * x1).^(-3);
            M.dM11dx2 = zeros(size(x1));
            M.dM22dx1 = zeros(size(x1));
            M.dM22dx2 = -1200000 * (1 + 15 * x2).^(-3);
        elseif M_type == 2
            M.M11 = 1000 + 600*sin(2*pi*x1).*sin(2*pi*x2);
            M.M22 = 1000 - 600*sin(2*pi*x1).*sin(2*pi*x2);
            M.dM11dx1 = 1200*pi*cos(2*pi*x1).*sin(2*pi*x2);
            M.dM11dx2 = 1200*pi*sin(2*pi*x1).*cos(2*pi*x2);
            M.dM22dx1 = -1200*pi*cos(2*pi*x1).*sin(2*pi*x2);
            M.dM22dx2 = -1200*pi*sin(2*pi*x1).*cos(2*pi*x2);
        elseif M_type == 3
            [Nx2, Nx1] = size(x1);
            M.M11 = 2000*ones(Nx2, Nx1);
            M.M22 = 2000*ones(Nx2, Nx1);
            [Nx2, Nx1] = size(x1);
            M.dM11dx1 = zeros(Nx2, Nx1);
            M.dM11dx2 = zeros(Nx2, Nx1);
            M.dM22dx1 = zeros(Nx2, Nx1);
            M.dM22dx2 = zeros(Nx2, Nx1);
        elseif M_type == 4
            [Nx2, Nx1] = size(x1);
            M.M11 = 400*ones(Nx2, Nx1);
            M.M22 = 400*ones(Nx2, Nx1);
            [Nx2, Nx1] = size(x1);
            M.dM11dx1 = zeros(Nx2, Nx1);
            M.dM11dx2 = zeros(Nx2, Nx1);
            M.dM22dx1 = zeros(Nx2, Nx1);
            M.dM22dx2 = zeros(Nx2, Nx1);
        end
    else
        M.M11 = M_type.FM11(x1(:), x2(:)); 
        M.M11 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.M11_samp, x1, x2, 'linear');
        M.M12 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.M12_samp, x1, x2, 'linear');
        M.M22 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.M22_samp, x1, x2, 'linear');

        M.dM11dx1 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.dM11dx1_samp, x1, x2, 'linear');
        M.dM12dx1 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.dM12dx1_samp, x1, x2, 'linear');
        M.dM22dx1 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.dM22dx1_samp, x1, x2, 'linear');

        M.dM11dx2 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.dM12dx2_samp, x1, x2, 'linear');
        M.dM12dx2 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.dM12dx2_samp, x1, x2, 'linear');
        M.dM22dx2 = interp2(M_type.x1_samp, M_type.x2_samp, M_type.dM22dx2_samp, x1, x2, 'linear');
    end
end