clear

addpath('./','./utils')
animation = 1;
pause_time = 0.1;
make_gif = 0;
problem = 2;
method = 2;
gif_name = 'example21.gif';
title_name = 'Curved metric, Alternative Cost function';

max_iter = 100;
tolerance = 1e-6;
epsilon = 0.2; % Under-relaxation factor

% Grid size
Nx1 = 40;
Nx2 = 40;
N = Nx1*Nx2;

% Equispaced computational coordinates
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);

% Initial physical coordinates (equispaced)
[x1, x2] = meshgrid(s1, s2);

% Plot initial grid
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
title(title_name); xlim([0,1]); ylim([0,1]);
axis equal; hold off
if make_gif; exportgraphics(gcf, gif_name); end

[M11, M22] = M(x1, x2, problem);

res_list = [];

for iter = 1:max_iter
    [A, b, res] = AssembleLinearSystem(x1, x2, problem, method);
    x_star = A \ b;

    x1_star = x_star(1:N);
    x2_star = x_star(N+1:end);

    dx1 = reshape(x1_star, Nx2, Nx1) - x1;
    dx2 = reshape(x2_star, Nx2, Nx1) - x2;

    res_list = [res_list, res];
    x1 = x1 + epsilon * dx1;
    x2 = x2 + epsilon * dx2;
    
    if animation == 1
        plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
        title(title_name); xlim([0,1]); ylim([0,1]);
        axis equal; hold off;
        pause(pause_time);
    end
    if make_gif; exportgraphics(gcf, gif_name, Append=true); end
end

%%
figure()
plot(1:iter, res_list);
legend('Cost', 'Residual');
grid on

%%
function [A, b, res] = AssembleLinearSystem(x1, x2, problem, method)
    [Nx1, Nx2] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;
    [M11, M22] = M(x1, x2, problem);
    [dM11dx1, dM11dx2, dM22dx1, dM22dx2] = dMdx(x1, x2, problem);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, sigma1, sigma2);
    id = GetIndex(Nx1, Nx2);

    A_orth = AssembleOrtho(dx1ds1, dx1ds2, dx2ds1, dx2ds2, sigma1, sigma2, id);
    A = sparse(N_all, N_all);
    
    if method == 2 % For approximate cost function
        % coef(k,alpha,i)
        coef111 = -2*sigma1^2*M11.*dx1ds1.*M11.*dx1ds1; coef111 = coef111(:);
        coef112 = -2*sigma1^2*M11.*dx1ds1.*M22.*dx2ds1; coef112 = coef112(:);
        coef121 = -2*sigma2^2*M11.*dx1ds2.*M11.*dx1ds2; coef121 = coef121(:);
        coef122 = -2*sigma2^2*M11.*dx1ds2.*M22.*dx2ds2; coef122 = coef122(:);

        coef211 = -2*sigma1^2*M22.*dx2ds1.*M11.*dx1ds1; coef211 = coef211(:);
        coef212 = -2*sigma1^2*M22.*dx2ds1.*M22.*dx2ds1; coef212 = coef212(:);
        coef221 = -2*sigma2^2*M22.*dx2ds2.*M11.*dx1ds2; coef221 = coef221(:);
        coef222 = -2*sigma2^2*M22.*dx2ds2.*M22.*dx2ds2; coef222 = coef222(:);

        % first component of s
        for i = [id.b, id.inner, id.t]
            % Same way to handle ds1
            A(i,i) = A(i,i) - 2*coef111(i);  % (i,i) for x1 contribution on x1
            A(i,i-Nx2) = A(i,i-Nx2) + coef111(i);
            A(i,i+Nx2) = A(i,i+Nx2) + coef111(i);

            A(i,i+N) = A(i,i+N) - 2*coef112(i); % (i,i+N) for x2 contribution on x1
            A(i,i+N-Nx2) = A(i,i+N-Nx2) + coef112(i);
            A(i,i+N+Nx2) = A(i,i+N+Nx2) + coef112(i);

            % Different way to handle ds2
            if ismember(i,id.inner) % Central difference at interior
                A(i,i) = A(i,i) - 2*coef121(i);
                A(i,i-1) = A(i,i-1) + coef121(i);
                A(i,i+1) = A(i,i+1) + coef121(i);
                A(i,i+N) = A(i,i+N) - 2*coef122(i);
                A(i,i+N-1) = A(i,i+N-1) + coef122(i);
                A(i,i+N+1) = A(i,i+N+1) + coef122(i);
            elseif ismember(i,id.b) % One sided difference at boundary
                A(i,i) = A(i,i) - 2*coef121(i);
                A(i,i+1) = A(i,i+1) + 2*coef121(i);

                A(i,i+N) = A(i,i+N) + 2*coef122(i);
                A(i,i+N+1) = A(i,i+N+1) - 5*coef122(i);
                A(i,i+N+2) = A(i,i+N+2) + 4*coef122(i);
                A(i,i+N+3) = A(i,i+N+3) - 1*coef122(i);
            elseif ismember(i,id.t)
                A(i,i) = A(i,i) - 2*coef121(i);
                A(i,i-1) = A(i,i-1) + 2*coef121(i);

                A(i,i+N) = A(i,i+N) + 2*coef122(i);
                A(i,i+N-1) = A(i,i+N-1) - 5*coef122(i);
                A(i,i+N-2) = A(i,i+N-2) + 4*coef122(i);
                A(i,i+N-3) = A(i,i+N-3) - 1*coef122(i);
            end
        end


        % second component of interior
        for i = N+[id.l, id.inner, id.r]
            A(i,i) = A(i,i) - 2*coef222(i-N); % (i,i-N) for x1 contribution on x2
            A(i,i-1) = A(i,i-1) + coef222(i-N);
            A(i,i+1) = A(i,i+1) + coef222(i-N);

            A(i,i-N) = A(i,i-N) - 2*coef221(i-N);  % (i,i) for x2 contribution on x2
            A(i,i-N-1) = A(i,i-N-1) + coef221(i-N);
            A(i,i-N+1) = A(i,i-N+1) + coef221(i-N);

            if ismember(i-N,id.inner)
                A(i,i-N) = A(i,i-N) - 2*coef211(i-N);
                A(i,i-N-Nx2) = A(i,i-N-Nx2) + coef211(i-N);
                A(i,i-N+Nx2) = A(i,i-N+Nx2) + coef211(i-N);
    
                A(i,i) = A(i,i) - 2*coef212(i-N);
                A(i,i-Nx2) = A(i,i-Nx2) + coef212(i-N);
                A(i,i+Nx2) = A(i,i+Nx2) + coef212(i-N);
            elseif ismember(i-N,id.l) % One sided difference for dx1/ds1
                A(i,i) = A(i,i) - 2*coef212(i-N);
                A(i,i+Nx2) = A(i,i+Nx2) + 2*coef212(i-N);

                A(i,i-N) = A(i,i-N) + 2*coef211(i-N);
                A(i,i-N+Nx2) = A(i,i-N+Nx2) - 5*coef211(i-N);
                A(i,i-N+2*Nx2) = A(i,i-N+2*Nx2) + 4*coef211(i-N);
                A(i,i-N+3*Nx2) = A(i,i-N+3*Nx2) - 1*coef211(i-N);
            elseif ismember(i-N,id.r)
                A(i,i) = A(i,i) - 2*coef212(i-N);
                A(i,i-Nx2) = A(i,i-Nx2) + 2*coef212(i-N);

                A(i,i-N) = A(i,i-N) + 2*coef211(i-N);
                A(i,i-N-Nx2) = A(i,i-N-Nx2) - 5*coef211(i-N);
                A(i,i-N-2*Nx2) = A(i,i-N-2*Nx2) + 4*coef211(i-N);
                A(i,i-N-3*Nx2) = A(i,i-N-3*Nx2) - 1*coef211(i-N);
            end
        end
    
        %A = A + 1*A_orth;
        
        % For Dirichlet conditions
        for i =   [id.l, id.r, id.corner]; A(i,:) = 0; A(i,i) = 1; end
        for i = N+[id.b, id.t, id.corner]; A(i,:) = 0; A(i,i) = 1; end
    
        b_aux1 = dM11dx1.*dx1ds1+dM11dx2.*dx2ds1;
        b_aux2 = dM22dx1.*dx1ds1+dM22dx2.*dx2ds1;
        b_aux3 = dM11dx1.*dx1ds2+dM11dx2.*dx2ds2;
        b_aux4 = dM22dx1.*dx1ds2+dM22dx2.*dx2ds2;
        b1 = sigma1^4 * M11.*dx1ds1.*(b_aux1.*dx1ds1.^2 + b_aux2.*dx2ds1.^2) + ...
             sigma2^4 * M11.*dx1ds2.*(b_aux3.*dx1ds2.^2 + b_aux4.*dx2ds2.^2);
        b2 = sigma1^4 * M22.*dx2ds1.*(b_aux1.*dx1ds1.^2 + b_aux2.*dx2ds1.^2) + ...
             sigma2^4 * M22.*dx2ds2.*(b_aux3.*dx1ds2.^2 + b_aux4.*dx2ds2.^2);

    elseif method == 1 % For alternative cost function
        for i = [id.inner, N+id.inner]  % N+id.inner for the second component
            A(i,i) = -4;
            A(i,i-1) = 1;
            A(i,i+1) = 1;
            A(i,i-Nx2) = 1;
            A(i,i+Nx2) = 1;
        end
    
        for i = N+[id.l, id.r]  % s2 component of left and right boundary: replace Laciancian by dy
            A(i,i) = -4;
            A(i,i-1) = 1;
            A(i,i+1) = 1;
            if ismember(i-N,id.l); A(i,i+Nx2) = 2; end
            if ismember(i-N,id.r); A(i,i-Nx2) = 2; end
        end
    
        for i = [id.b, id.t]  % s1 component of top and bottom boundary: replace Laciancian by dx
            A(i,i) = -4;
            A(i,i-Nx2) = 1;
            A(i,i+Nx2) = 1;
            if ismember(i,id.b); A(i,i+1) = 2; end
            if ismember(i,id.t); A(i,i-1) = 2; end
        end
        Mii = [M11(:); M22(:)];
        A = -2*A.*Mii;

        %A = A + 0.*A_orth;
    
        % For Dirichlet conditions
        for i =   [id.l, id.r, id.corner]; A(i, i) = 1; end
        for i = N+[id.b, id.t, id.corner]; A(i, i) = 1; end
    
        % Assemble the vector b
        b1 = dM11dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
        b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
        b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
    
        b2 = dM22dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
        b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
        b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    end
    % Apply x1=0 and x1=1 to the left and right boundary
    b1(:,1) = 0;
    b1(:,end) = 1;
    % Apply x2=0 and x2=1 to the top and bottom boundary
    b2(1,:) = 0;
    b2(end,:) = 1;
    b = [b1(:); b2(:)];

    res = norm(A*[x1(:);x2(:)] - b);
end

function A_orth = AssembleOrtho(dx1ds1, dx1ds2, dx2ds1, dx2ds2, sigma1, sigma2, id)
    [Nx1, Nx2] = size(dx1ds1);
    N = Nx1*Nx2;
    N_all = N*2;
    A_orth1 = sparse(N_all, N_all);
    A_orth2 = sparse(N_all, N_all);
    A_orth3 = sparse(N_all, N_all);
    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;

    D = dx1ds1.*dx1ds2 + dx2ds1.*dx2ds2;
    [dDds1, dDds2] = DCentral(D, D, sigma1, sigma2);
    
    for i = [id.inner, N+id.inner]
        A_orth1(i,i+Nx2+1) = 1/(4*sigma1*sigma2);
        A_orth1(i,i-Nx2+1) = -1/(4*sigma1*sigma2);
        A_orth1(i,i+Nx2-1) = -1/(4*sigma1*sigma2);
        A_orth1(i,i-Nx2-1) = 1/(4*sigma1*sigma2);
    end

    for i = N+[id.l, id.r, -N+id.inner, id.inner]
        A_orth2(i,i+1) = 1/(2*sigma2);
        A_orth2(i,i-1) = -1/(2*sigma2);
    end

    for i = [id.b, id.t, N+id.inner, id.inner]
        A_orth3(i,i+Nx2) = 1/(2*sigma1);
        A_orth3(i,i-Nx2) = -1/(2*sigma1);
    end

    A_orth = -4*A_orth1.*[D(:);D(:)] - 2*A_orth2.*[dDds1(:);dDds1(:)] - 2*A_orth3.*[dDds2(:);dDds2(:)];
end
