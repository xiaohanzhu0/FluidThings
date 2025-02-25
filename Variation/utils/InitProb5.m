function [T1, T2, T1_, T2_, M] = InitProb5(Nx1, Nx2)
    addpath('~/Files/data/Mesh_generation/Airfoil');
    
    % Processing sampled metric tensor point
    m = load('metric_airfoil_bodyfitted.mat');
    [M.dM11dx1_samp, M.dM11dx2_samp] = DCentralUneven(m.metric(:,:,1), m.metric(:,:,1), ...
                                        m.y_metric, m.x_metric);
    
    [M.dM12dx1_samp, M.dM12dx2_samp] = DCentralUneven(m.metric(:,:,2), m.metric(:,:,2), ...
                                        m.y_metric, m.x_metric);
    
    [M.dM22dx1_samp, M.dM22dx2_samp] = DCentralUneven(m.metric(:,:,3), m.metric(:,:,3), ...
                                        m.y_metric, m.x_metric);
    
    M.M11_samp = m.metric(:,:,1);
    M.M12_samp = m.metric(:,:,2);
    M.M22_samp = m.metric(:,:,3);
    M.x1_samp = m.x_metric;
    M.x2_samp = m.y_metric;

    M.FM11 = scatteredInterpolant(M.x1_samp(:),M.x2_samp(:),M.M11_samp(:));
    
    %% Processing sampled boundary points
    Ns2_raw = Nx2;
    
    b1_raw = [M.x1_samp(:,1), M.x2_samp(:,1)];
    b3_raw = [M.x1_samp(:,end), M.x2_samp(:,end)];
    b2_raw = [b1_raw(end,:); b3_raw(end,:)];
    b4_raw = [b1_raw(1,:); b3_raw(1,:)];
    b2_raw = [linspace(b2_raw(1,1),b2_raw(2,1),Ns2_raw)', linspace(b2_raw(1,2),b2_raw(2,2),Ns2_raw)'];
    b4_raw = [linspace(b4_raw(1,1),b4_raw(2,1),Ns2_raw)', linspace(b4_raw(1,2),b4_raw(2,2),Ns2_raw)'];
    
    Ns1_raw = size(b1_raw, 1);
    
    
    F = (1-abs(linspace(-1,1,Ns1_raw)')).^(1.5).*vecnorm([gradient(b1_raw(:,1)),gradient(b1_raw(:,2))],2,2);
    phi_solution = BoundaryMesh(F, Nx1);
    x1 = interp1(linspace(0,1,Ns1_raw), b1_raw(:,1), phi_solution, 'pchip');
    x2 = interp1(linspace(0,1,Ns1_raw), b1_raw(:,2), phi_solution, 'pchip');
    b1 = [x1', x2'];
    
    F = (1-abs(linspace(-1,1,Ns1_raw)')).^(1/4).*vecnorm([gradient(b3_raw(:,1)),gradient(b3_raw(:,2))],2,2);
    phi_solution = BoundaryMesh(F, Nx1);
    x1 = interp1(linspace(0,1,Ns1_raw), b3_raw(:,1), phi_solution, 'pchip');
    x2 = interp1(linspace(0,1,Ns1_raw), b3_raw(:,2), phi_solution, 'pchip');
    b3 = [x1', x2'];
    
    F = exp(2*linspace(0,1,Ns2_raw)').*vecnorm([gradient(b2_raw(:,1)),gradient(b2_raw(:,2))],2,2);
    phi_solution = BoundaryMesh(F, Nx2);
    x1 = interp1(linspace(0,1,Ns2_raw), b2_raw(:,1), phi_solution, 'pchip');
    x2 = interp1(linspace(0,1,Ns2_raw), b2_raw(:,2), phi_solution, 'pchip');
    b2 = [x1', x2'];
    
    F = exp(2*linspace(0,1,Ns2_raw)').*vecnorm([gradient(b4_raw(:,1)),gradient(b4_raw(:,2))],2,2);
    phi_solution = BoundaryMesh(F, Nx2);
    x1 = interp1(linspace(0,1,Ns2_raw), b4_raw(:,1), phi_solution, 'pchip');
    x2 = interp1(linspace(0,1,Ns2_raw), b4_raw(:,2), phi_solution, 'pchip');
    b4 = [x1', x2'];

    [T1, T2] = TFI(b1', flip(b2',2), b3', flip(b4',2));
    T1_ = T1;
    T2_ = T1;

    %scatter(b1(:,1),b1(:,2)); hold on
    %scatter(b2(:,1),b2(:,2)); hold on
    %scatter(b3(:,1),b3(:,2)); hold on
    %scatter(b4(:,1),b4(:,2)); hold on
    %plot(b3_raw(:,1),b3_raw(:,2)); hold on
    
    %[T1, T2] = TFI(b1', flip(b2',2), b3', flip(b4',2));
    %for i=1:Nx2
    %    plot(T1(i,:),T2(i,:)); hold on
    %end
    %for i=1:Nx1
    %    plot(T1(:,i),T2(:,i))
    %end
end

%% Apply TFI
function [T1, T2] = TFI(bottom, right, top, left)
    Nx1 = size(bottom, 2);
    Nx2 = size(left, 2);
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