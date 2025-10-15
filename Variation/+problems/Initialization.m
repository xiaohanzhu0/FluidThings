function [x1,x2,Mfun] = Initialization(problem,Nx1,Nx2)
    if problem == 1
        [x1,x2] = problems.InitProb1(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob1Metric(x1, x2);
    elseif problem == 2
        [x1,x2] = problems.InitProb2(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob2Metric(x1, x2);
    elseif problem == 3
        [x1,x2] = problems.InitProb3(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob3Metric(x1, x2);
    elseif problem == 4
        [x1,x2] = problems.InitProb4(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob4Metric(x1, x2);
    elseif problem == 5
        cf.new_airfoil = 1;
        cf.metric_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/metricField.fields';
        cf.airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';
        cf.Nx1 = Nx1; cf.Nx2 = Nx2;
        cf.alpha = 1.005;
        cf.append_trail = 0;
        [x1, x2, M_samp, Mfun] = problems.InitProb5(cf);
    elseif problem == 6
        [x1,x2] = problems.InitProb6(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob6Metric(x1, x2);
    elseif problem == 8
        [x1,x2] = problems.InitProb8(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob8Metric(x1, x2);
    end
end