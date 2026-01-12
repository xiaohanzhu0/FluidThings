function Mfun = GetM(problem)
    if problem == 1
        Mfun = @(x1, x2) problems.Prob1Metric(x1,x2);
    elseif problem == 2
        Mfun = @(x1, x2) problems.Prob2Metric(x1,x2);
    elseif problem == 3
        Mfun = @(x1, x2) problems.Prob3Metric(x1,x2);
    elseif problem == 4
        Mfun = @(x1, x2) problems.Prob4Metric(x1,x2);
    elseif problem == 5
        Mfun = @(x1, x2) problems.Prob5Metric(x1,x2);
    elseif problem == 6
        Mfun = @(x1, x2) problems.Prob6Metric(x1,x2);
    elseif problem == 8
        Mfun = @(x1, x2) problems.Prob8Metric(x1,x2);
    end
end
