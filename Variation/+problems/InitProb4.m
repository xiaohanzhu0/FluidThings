function [T1, T2] = InitProb4(Nx1, Nx2)
    e1 = ones(1,Nx1);
    e2 = ones(1,Nx2);
    
    N_line = floor(Nx1/3);
    N_circ = Nx1 - 2*N_line;
    bottom1 = [linspace(-1.5,1.5,Nx1); 0*e1];
    bottom2 = 1/3*[1.5*cos(linspace(pi,0,N_circ)); 1.5*sin(linspace(pi,0,N_circ))];
    bottom = [bottom1(:,1:N_line), bottom2, bottom1(:,N_line+N_circ+1:end)];
    right = [1.5*e2; linspace(0,2,Nx2)];
    top = [linspace(-1.5,1.5,Nx1); 2*e1];
    left = [-1.5*e2; linspace(0,2,Nx2)];
    
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