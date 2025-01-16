function [dM11dx1, dM11dx2, dM22dx1, dM22dx2] = dMdx(x1, x2, problem)
    if problem == 1
        dM11dx1 = -1200000 * (1 + 15 * x1).^(-3);
        dM11dx2 = zeros(size(x1));
        dM22dx1 = zeros(size(x1));
        dM22dx2 = -1200000 * (1 + 15 * x2).^(-3);
    elseif problem == 2
        dM11dx1 = 1200*pi*cos(2*pi*x1).*sin(2*pi*x2);
        dM11dx2 = 1200*pi*sin(2*pi*x1).*cos(2*pi*x2);
        dM22dx1 = -1200*pi*cos(2*pi*x1).*sin(2*pi*x2);
        dM22dx2 = -1200*pi*sin(2*pi*x1).*cos(2*pi*x2);
    end
end