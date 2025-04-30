nonlinear_list = [2, 3];
epsilon_list = 0.3:0.1:2;
C = 1;
converge_list = zeros(2,length(epsilon_list));

for i = 1:2
    for j = 1:length(epsilon_list)
        nonlinear = nonlinear_list(i);
        epsilon = epsilon_list(j);
        res_list = CurvedBv6(nonlinear, epsilon, C);
        if isnan(res_list(end))
            converge_list(i,j) = NaN;
        else
            converge_list(i,j) = length(res_list);
        end
        disp(['Finished epsilon=', num2str(epsilon)]);
    end
end
%converge_list(converge_list>=201) = NaN;
%converge_list(converge_list==0) = NaN;

plot(epsilon_list, converge_list, '-o', LineWidth=2);
legend('nonlinear 1', 'nonlinear 2');
title('Number of iteration to reach 1e-6 vs relaxation')
grid on