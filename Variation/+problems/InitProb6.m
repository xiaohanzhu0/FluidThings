function [x1,x2] = InitProb6(Nx1,Nx2)
    x1 = linspace(0,1,Nx1);
    x2 = linspace(0,1,Nx2);
    [x1_temp, x2_temp] = meshgrid(x1, x2);
    theta = pi/6;
    x1 = x1_temp*cos(theta) - x2_temp*sin(theta);
    x2 = x1_temp*sin(theta) + x2_temp*cos(theta);
end

