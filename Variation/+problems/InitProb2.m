function [x1,x2] = InitProb2(Nx1,Nx2)
    x1 = linspace(0,1,Nx1);
    x2 = linspace(0,1,Nx2);
    [x1, x2] = meshgrid(x1, x2);
end