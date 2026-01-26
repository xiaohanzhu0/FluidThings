function [metric, update_history] = apply_metric_gradation(x_metric, y_metric, metric, params)
%APPLY_METRIC_GRADATION Apply local gradation to a structured metric field.
% Returns update_history to show whether each outer iteration changed M.
    if nargin < 4 || isempty(params)
        params = struct();
    end
    if ~isfield(params, 'beta') || isempty(params.beta)
        params.beta = 1.002;
    end
    if ~isfield(params, 'grade_iterations') || isempty(params.grade_iterations)
        params.grade_iterations = 20;
    end
    if ~isfield(params, 'omega') || isempty(params.omega)
        params.omega = 1;
    end

    tol = 1e-12;
    nbr = [-1 0; 1 0; 0 -1; 0 1];
    [Nx1, Nx2] = size(x_metric);
    x = cat(1, reshape(x_metric, 1, Nx1, Nx2), reshape(y_metric, 1, Nx1, Nx2));
    M = zeros(2, 2, Nx1, Nx2);
    M(1,1,:,:) = metric(:,:,1);
    M(1,2,:,:) = metric(:,:,2);
    M(2,1,:,:) = metric(:,:,2);
    M(2,2,:,:) = metric(:,:,3);

    i_list = randperm(Nx1);
    j_list = randperm(Nx2);
    max_iter = params.grade_iterations;
    update_history = false(1, max_iter);
    for iter = 1:max_iter
        did_update = false;
        %if mod(iter,2) == 1; i_list = 1:Nx1; j_list = 1:Nx2;
        %else i_list = Nx1:-1:1; j_list = Nx2:-1:1; end

        for i = i_list
            for j = j_list
                % Skip adjancent points to avoid redundant edge iterations
                if mod(i+j,2) == 1; continue; end
                p = x(:,i,j);

                for k = 1:size(nbr,1)
                    ii = i + nbr(k,1);
                    jj = j + nbr(k,2);

                    if ii < 1 || ii > Nx1 || jj < 1 || jj > Nx2
                        continue
                    end

                    M_p = M(:,:,i,j);
                    q = x(:,ii,jj);
                    M_q = M(:,:,ii,jj);
                    pq = q - p;

                    eta = (1 + sqrt(pq' * M_p * pq) * log(params.beta))^(-2);
                    M_xp = eta * M_p;
                    M_q_new = real(metric_intersection(M_xp, M_q));
                    if any(abs(M_q_new(:) - M_q(:)) > tol)
                        M(:,:,ii,jj) = M_q_new;
                        did_update = true;
                    end

                    eta = (1 + sqrt(pq' * M_q * pq) * log(params.beta))^(-2);
                    M_xq = eta * M_q;
                    M_p_new = real(metric_intersection(M_xq, M_p));
                    if any(abs(M_p_new(:) - M_p(:)) > tol)
                        M(:,:,i,j) = M_p_new;
                        did_update = true;
                    end
                end
            end
        end
        update_history(iter) = did_update;
        if ~did_update
            update_history = update_history(1:iter);
            break;
        elseif iter == max_iter
            update_history = update_history(1:iter);
        end
    end

    metric(:,:,1) = reshape(M(1,1,:,:), Nx1, Nx2);
    metric(:,:,2) = reshape(M(1,2,:,:), Nx1, Nx2);
    metric(:,:,3) = reshape(M(2,2,:,:), Nx1, Nx2);
end

function M_inter = metric_intersection(M_xp, M_q)
    N = M_xp \ M_q;
    [V, ~] = eig(N);
    V(:,1) = V(:,1) / norm(V(:,1));
    V(:,2) = V(:,2) / norm(V(:,2));
    lam_xp = V' * M_xp * V;
    lam_q = V' * M_q * V;
    lam1_inter = max(lam_xp(1,1), lam_q(1,1));
    lam2_inter = max(lam_xp(2,2), lam_q(2,2));

    M_inter = inv(V') * diag([lam1_inter, lam2_inter]) * inv(V);
end
