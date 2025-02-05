function [T1, T2] = InitGuess(Nx1, Nx2, A)
    e01 = linspace(0,1,Nx1);
    e10 = linspace(1,0,Nx1);
    
    bottom = [e01; -A*sin(2*pi*e01)];
    right = [1+A*sin(2*pi*e01); e01];
    top = [e01; 1-A*sin(2*pi*e01)];
    left = [A*sin(2*pi*e01); e01];
    
    rb = right(:,1);
    rt = right(:,end);
    lt = left(:,end);
    lb = left(:,1);
    
    s1 = linspace(0,1,Nx1);
    s2 = linspace(0,1,Nx2);
    [S1, S2] = meshgrid(s1, s2);
    
    T1 = (1-S2).*bottom(1,:) + S2.*top(1,:) + (1-S1).*right(1,:)' + S1.*left(1,:)' ...
    - (1-S1).*(1-S2)*rb(1) - S1.*S2*lt(1) - S1.*(1-S2)*lb(1) - (1-S1).*S2*rt(1);
    
    T2 = (1-S2).*bottom(2,:) + S2.*top(2,:) + (1-S1).*right(2,:)' + S1.*left(2,:)' ...
    - (1-S1).*(1-S2)*rb(2) - S1.*S2*lt(2) - S1.*(1-S2)*lb(2) - (1-S1).*S2*rt(2);
end