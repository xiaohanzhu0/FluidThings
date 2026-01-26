function [T1, T2, M11_fun,M12_fun,M22_fun] = InitProb9(cf)
        [data, ~, ~, ~, ~] = readTurtleFields(cf.metric_datapath);
        [x_cv, ~, ~, ~, ~, ~, ~, ~, ~, ~] = readTurtleGrid(cf.airfoil_datapath);

        M_samp.x_metric = x_cv{1,cf.block_idx}(:,:,1,1);
        M_samp.y_metric = x_cv{1,cf.block_idx}(:,:,1,2);

        M = data{1,1}{1,cf.block_idx}(:,:,1,1:2,1:2);
        dim = size(M);
        M = reshape(M, dim(2), dim(1), dim(4), dim(5));

        M_samp.metric(:,:,1) = M(:,:,1,1);
        M_samp.metric(:,:,2) = M(:,:,1,2);
        M_samp.metric(:,:,3) = M(:,:,2,2);

        M_samp.x_metric = M_samp.x_metric';
        M_samp.y_metric = M_samp.y_metric';
    
    %%
    if cf.grade
        [M_samp.metric, update_history] = problems.apply_metric_gradation(...
            M_samp.x_metric, M_samp.y_metric, M_samp.metric,cf);
    end

    M11_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
    M12_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
    M22_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));

    T1 = M_samp.x_metric;
    T2 = M_samp.y_metric;

end
