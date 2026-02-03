function [x1,x2,M11_fun,M12_fun,M22_fun] = Initialization(problem, Nx1, Nx2, params)
    if nargin < 4 || isempty(params)
        params = struct('problemId', problem, 'Nx1', Nx1, 'Nx2', Nx2);
    else
        if ~isfield(params, 'problemId') || isempty(params.problemId)
            params.problemId = problem;
        end
        if ~isfield(params, 'Nx1') || isempty(params.Nx1)
            params.Nx1 = Nx1;
        end
        if ~isfield(params, 'Nx2') || isempty(params.Nx2)
            params.Nx2 = Nx2;
        end
    end

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
        cf = meshgen.defaults_dual(params);
        [x1, x2, M11_fun, M12_fun, M22_fun] = problems.InitProb5(cf);
    elseif problem == 6
        [x1,x2] = problems.InitProb6(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob6Metric(x1, x2);
    elseif problem == 8
        [x1,x2] = problems.InitProb8(Nx1,Nx2);
        Mfun = @(x1,x2) problems.Prob8Metric(x1, x2);
    elseif problem == 9
        cf = meshgen.defaults_dual(params);
        [x1, x2, M11_fun,M12_fun,M22_fun] = problems.InitProb9(cf);
    end
end
