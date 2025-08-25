% x1 and x2: the boundary points of the domain in 4*N array
% 4 ---- 3
% |      |
% 1 ---- 2
function [x1_new, x2_new] = UpdateCorrection(x1, x2, boundary_points)
    x1_new = x1;
    x2_new = x2;
    % For bottom boundary
    bx1 = x1(1,:);
    bx2 = x2(1,:);
    bx1_exact = boundary_points.b(1,:);
    bx2_exact = boundary_points.b(2,:);
    [bx1_new, bx2_new] = UpdateEachSide(bx1, bx2, bx1_exact, bx2_exact);
    x1_new(1,:) = bx1_new;
    x2_new(1,:) = bx2_new;
   
    % For top boundary
    bx1 = x1(end,:);
    bx2 = x2(end,:);
    bx1_exact = boundary_points.t(1,:);
    bx2_exact = boundary_points.t(2,:);
    [bx1_new, bx2_new] = UpdateEachSide(bx1, bx2, bx1_exact, bx2_exact);
    x1_new(end,:) = bx1_new;
    x2_new(end,:) = bx2_new;

    % For left boundary
    bx1 = x1(:,1);
    bx2 = x2(:,1);
    bx1_exact = boundary_points.l(:,1);
    bx2_exact = boundary_points.l(:,2);
    [bx1_new, bx2_new] = UpdateEachSide(bx1, bx2, bx1_exact, bx2_exact);
    x1_new(:,1) = bx1_new;
    x2_new(:,1) = bx2_new;

    % For right boundary
    bx1 = x1(:,end);
    bx2 = x2(:,end);
    bx1_exact = boundary_points.r(:,1);
    bx2_exact = boundary_points.r(:,2);
    [bx1_new, bx2_new] = UpdateEachSide(bx1, bx2, bx1_exact, bx2_exact);
    x1_new(:,end) = bx1_new;
    x2_new(:,end) = bx2_new;
end


function [bx1_new, bx2_new] = UpdateEachSide(bx1, bx2, bx1_exact, bx2_exact)
    N = length(bx1);
    N_exact = length(bx1_exact);
    bx1_new = zeros(size(bx1));
    bx2_new = zeros(size(bx2));
    
    for i = 1:N
        x1_star = bx1(i);
        x2_star = bx2(i);
        min_dist = Inf;
        best_idx = 0;
        best_r = 0;
        for j = 1:N_exact-1
            % Vector along the segment
            v = [bx1_exact(j+1) - bx1_exact(j); bx2_exact(j+1) - bx2_exact(j)];
            % Vector from the start of the segment to the target
            w = [x1_star - bx1_exact(j); x2_star - bx2_exact(j)];
            % Projection factor (can be <0 or >1)
            r = dot(w, v) / dot(v, v);
            % Clamp r to [0,1]
            r = max(0, min(1, r));
            % Compute the projection point on the segment
            proj = [bx1_exact(j); bx2_exact(j)] + r * v;
            % Distance from target to this projection
            dist = norm(proj - [x1_star; x2_star]);

            if dist < min_dist
                min_dist = dist;
                best_idx = j;
                best_r = r;
            end
        end
        bx1_new(i) = bx1_exact(best_idx) + best_r*(bx1_exact(best_idx+1) - bx1_exact(best_idx));
        bx2_new(i) = bx2_exact(best_idx) + best_r*(bx2_exact(best_idx+1) - bx2_exact(best_idx));
    end
end