function [T1, T2, M11_fun, M12_fun, M22_fun] = InitProb5(cf)

    [data, ~, ~, ~, ~] = readTurtleFields(cf.metric_datapath);
    [x_cv, ~, ~, ~, ~, ~, ~, ~, ~, ~] = readTurtleGrid(cf.airfoil_datapath);
    M_samp.x_metric = [x_cv{1,9}(:,:,1,1), x_cv{1,8}(:,:,1,1);
    x_cv{1,10}(:,:,1,1),x_cv{1,1}(:,:,1,1);
    x_cv{1,11}(:,:,1,1),x_cv{1,2}(:,:,1,1);
    x_cv{1,12}(:,:,1,1),x_cv{1,3}(:,:,1,1);
    x_cv{1,13}(:,:,1,1),x_cv{1,4}(:,:,1,1);
    x_cv{1,14}(:,:,1,1),x_cv{1,5}(:,:,1,1);
    x_cv{1,15}(:,:,1,1),x_cv{1,6}(:,:,1,1);
    x_cv{1,16}(:,:,1,1),x_cv{1,7}(:,:,1,1)];

    M_samp.y_metric = [x_cv{1,9}(:,:,1,2), x_cv{1,8}(:,:,1,2);
    x_cv{1,10}(:,:,1,2),x_cv{1,1}(:,:,1,2);
    x_cv{1,11}(:,:,1,2),x_cv{1,2}(:,:,1,2);
    x_cv{1,12}(:,:,1,2),x_cv{1,3}(:,:,1,2);
    x_cv{1,13}(:,:,1,2),x_cv{1,4}(:,:,1,2);
    x_cv{1,14}(:,:,1,2),x_cv{1,5}(:,:,1,2);
    x_cv{1,15}(:,:,1,2),x_cv{1,6}(:,:,1,2);
    x_cv{1,16}(:,:,1,2),x_cv{1,7}(:,:,1,2)];

    aux1 = cat(1,data{1,1}{1,9}(:,:,1,1:2,1:2), data{1,1}{1,10}(:,:,1,1:2,1:2), ...
        data{1,1}{1,11}(:,:,1,1:2,1:2), data{1,1}{1,12}(:,:,1,1:2,1:2), ...
        data{1,1}{1,13}(:,:,1,1:2,1:2), data{1,1}{1,14}(:,:,1,1:2,1:2), ...
        data{1,1}{1,15}(:,:,1,1:2,1:2), data{1,1}{1,16}(:,:,1,1:2,1:2));
    aux2 = cat(1,data{1,1}{1,8}(:,:,1,1:2,1:2), data{1,1}{1,1}(:,:,1,1:2,1:2), ...
        data{1,1}{1,2}(:,:,1,1:2,1:2), data{1,1}{1,3}(:,:,1,1:2,1:2), ...
        data{1,1}{1,4}(:,:,1,1:2,1:2), data{1,1}{1,5}(:,:,1,1:2,1:2), ...
        data{1,1}{1,6}(:,:,1,1:2,1:2), data{1,1}{1,7}(:,:,1,1:2,1:2));
    aux = cat(2,aux1,aux2);
    dim = size(aux);
    M = reshape(aux, dim(1), dim(2), dim(4), dim(5));
    M_samp.metric(:,:,1) = M(:,:,1,1);
    M_samp.metric(:,:,2) = M(:,:,1,2);
    M_samp.metric(:,:,3) = M(:,:,2,2);
    
    
    %% Gradation to regularize the metric field
    if cf.grade
        M_samp.metric = problems.apply_metric_gradation(...
            M_samp.x_metric, M_samp.y_metric, M_samp.metric,cf);
    end

%%
    M_samp.F11 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
    M_samp.F12 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
    M_samp.F22 = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));
    Mfun = @(x1,x2) Prob5Metric(x1, x2, M_samp);
    function M = Prob5Metric(x1,x2,M_samp)
        M.M11 = M_samp.F11(x1,x2);
        M.M12 = M_samp.F12(x1,x2);
        M.M22 = M_samp.F22(x1,x2);
    end
    
M11_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,1),[],1));
M12_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,2),[],1));
M22_fun = scatteredInterpolant(M_samp.x_metric(:),M_samp.y_metric(:),reshape(M_samp.metric(:,:,3),[],1));

    if cf.hyperbolic == 1
        gd = Hyperbolic(cf.Nx1, cf.Nx2, cf.alpha, cf.append_trail);
        T1 = gd.x';
        T2 = gd.y';
    else
        T1 = M_samp.x_metric';
        T2 = M_samp.y_metric';
    end
end
