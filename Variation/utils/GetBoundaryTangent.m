function t = GetBoundaryTangent(x1, x2, normed)
    x1 = reshape(x1,1,[]);
    x2 = reshape(x2,1,[]);
    t = [x1(3:end)-x1(1:end-2); x2(3:end)-x2(1:end-2)];
    t0 = [x1(2)-x1(1); x2(2)-x2(1)];
    t1 = [x1(end)-x1(end-1); x2(end)-x2(end-1)];
    t = [t0, t, t1];
    if normed; t = normalize(t, 1, "norm"); end
end