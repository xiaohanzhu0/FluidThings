clear
addpath('./','./utils')
animation = 1;
pause_time = 0.0;
make_gif = 0;
problem = 1;
method = 1;
gif_name = 'example13.gif';
title_name = 'Curved boundary, Curved metric, Alternative Cost function';
save_output = 0;
march_type = 'regular';

C = 0.1;
Nx1 = 200;
Nx2 = 200;
max_iter = 200;
tolerance = 1e-6;
epsilon = 0.1; % Under-relaxation factor

% Grid size
N = Nx1*Nx2;

% Equispaced computational coordinates
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);

% Initial physical coordinates (equispaced)
if ismember(problem, [1,2])
    M_type = problem;
    %F = 40000 * (1 + 10*15 * s1).^(-1);
    %phi_solution = BoundaryMesh(F, Nx1);
    %x1 = interp1(linspace(0,1,Nx1), s1, phi_solution, 'pchip');
    %x2 = interp1(linspace(0,1,Nx2), s2, phi_solution, 'pchip');
    %[x1, x2] = meshgrid(x1, x2);

    %x1 = (-2./(s1-2) - 1).^4;
    %x2 = (-2./(s2-2) - 1).^4;
    %[x1, x2] = meshgrid(x1, x2);

    [x1, x2] = meshgrid(s1, s2);
    [x1_exact, x2_exact]  = meshgrid(s1, s2);
elseif problem == 3
    M_type = problem;
    [x1, x2] = InitProb3(Nx1, Nx2, C);
elseif problem == 4
    M_type = problem;
    [x1, x2] = InitProb4(Nx1, Nx2);
    [x1_exact, x2_exact] = InitProb4(50, 50);
elseif problem == 5
    [x1, x2, x1_exact, x2_exact, M_type] = InitProb5(Nx1, Nx2);
end

% Plot initial grid
%plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
%title(title_name); xlim([-30,30]); ylim([-30,30]);
%axis equal; hold off
%if make_gif; exportgraphics(gcf, gif_name); end

res_list = [];
err_list = [];
%%

for iter = 1:100
    M = GetM(x1, x2, M_type);
    [A, b, res] = AssembleLinearSystem(x1, x2, M, method);
    
    if strcmp(march_type, 'p-c')
        res = b - A*[x1(:);x2(:)];
        dx = epsilon*res;

        %err = A \ res;
        %dx = err;
        dx1 = reshape(dx(1:N), Nx2, Nx1);
        dx2 = reshape(dx(N+1:end), Nx2, Nx1);
        dx1(1,1)=0; dx1(1,end)=0; dx1(end,1)=0; dx1(end,end)=0;
        dx2(1,1)=0; dx2(1,end)=0; dx2(end,1)=0; dx2(end,end)=0;
        %prev_x1 = x1; prev_x2 = x2; 
        %k1_dx1 = epsilon * dx1; k1_dx2 = epsilon * dx2;
        x1 = x1 + dx1;
        x2 = x2 + dx2;
    end
    
    if strcmp(march_type, 'regular')
        tic
        x_star = A \ b;
        toc

        [L, U, P, Q, D] = lu(A);
        %tic
        %x_star_ = Q * (U \ (D \ (L \ (P * b))));
        %toc

        x1_star = x_star(1:N);
        x2_star = x_star(N+1:end);
        dx1 = reshape(x1_star, Nx2, Nx1) - x1;
        dx2 = reshape(x2_star, Nx2, Nx1) - x2;
        res_list = [res_list, norm(res)];

        dx1(1,1)=0; dx1(1,end)=0; dx1(end,1)=0; dx1(end,end)=0;
        dx2(1,1)=0; dx2(1,end)=0; dx2(end,1)=0; dx2(end,end)=0;
    
        x1 = x1 + epsilon * dx1;
        x2 = x2 + epsilon * dx2;
        
    end


    if problem == 1 || problem == 4 || problem == 5
        [x1, x2] = UpdateCorrection(x1, x2, x1_exact, x2_exact);
    end
    
    if animation == 1
        plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
        title(title_name); xlim([0,1]); ylim([0,1]);
        %plot(s1, x1(20,:));  xlim([0,1]); ylim([-0.1,1]); grid on
        axis equal; hold off;
        pause(pause_time);
    end
    if make_gif; exportgraphics(gcf, gif_name, Append=true); end
end

%%
figure()
semilogy(res_list); hold on
legend('Residual', 'Error');
grid on
