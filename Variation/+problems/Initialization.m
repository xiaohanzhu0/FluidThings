function [x1,x2,M11_fun,M12_fun,M22_fun] = Initialization(problem,Nx1,Nx2)
    if problem == 1
        [x1,x2] = problems.InitProb1(Nx1,Nx2);
        M11_fun = @(x1,x2) 40000*(1+15*x1).^(-2);
        M22_fun = @(x1,x2) 40000*(1+15*x2).^(-2);
        M12_fun = @(x1,x2) 0;
    elseif problem == 2
        [x1,x2] = problems.InitProb2(Nx1,Nx2);
        M11_fun = @(x1,x2) 1000 + 600*sin(2*pi*x1).*sin(2*pi*x2);
        M22_fun = @(x1,x2) 1000 - 600*sin(2*pi*x1).*sin(2*pi*x2);
        M12_fun = @(x1,x2) 0;
    elseif problem == 3
        [x1,x2] = problems.InitProb3(Nx1,Nx2);
        M11_fun = @(x1,x2) 2000;
        M22_fun = @(x1,x2) 2000;
        M12_fun = @(x1,x2) 0;
    elseif problem == 4
        [x1,x2] = problems.InitProb4(Nx1,Nx2);
        M11_fun = @(x1,x2) 400;
        M22_fun = @(x1,x2) 400;
        M12_fun = @(x1,x2) 0;
    elseif problem == 5
        cf.grade = 1;
        cf.new_airfoil = 1;
        cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
        cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
        cf.Nx1 = Nx1; cf.Nx2 = Nx2;
        cf.alpha = 1.005;
        cf.append_trail = 1;
        [x1, x2, M11_fun, M12_fun, M22_fun] = problems.InitProb5(cf);
    elseif problem == 6
        [x1,x2] = problems.InitProb6(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob6Metric(x1, x2);
    elseif problem == 8
        [x1,x2] = problems.InitProb8(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob8Metric(x1, x2);
    elseif problem == 9
        cf.block_idx = 5;
        cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
        cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
        [x1, x2, M11_fun,M12_fun,M22_fun] = problems.InitProb9(cf);
    end
end