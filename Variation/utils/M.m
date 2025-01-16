function [M11, M22] = M(x1, x2, problem)
    if problem == 1
        M11 = 40000 * (1 + 15 * x1).^(-2);
        M22 = 40000 * (1 + 15 * x2).^(-2);
    elseif problem == 2
        M11 = 1000 + 600*sin(2*pi*x1).*sin(2*pi*x2);
        M22 = 1000 - 600*sin(2*pi*x1).*sin(2*pi*x2);
    end
end