function [T1, T2, M_samp, Mfun] = InitProb9(cf)
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
    

    M_samp.F11 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
    M_samp.F12 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
    M_samp.F22 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));
    

    Mfun = @(x1,x2) Prob5Metric(x1, x2, M_samp);
    function M = Prob5Metric(x1,x2,M_samp)
        M.M11 = M_samp.F11(x1,x2);
        M.M12 = M_samp.F12(x1,x2);
        M.M22 = M_samp.F22(x1,x2);
    end

    T1 = M_samp.x_metric;
    T2 = M_samp.y_metric;

end

