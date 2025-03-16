function id = GetIndex(Nx1, Nx2)
    N = Nx1*Nx2;
    id.all = 1:N;
    id.l = 1:Nx2;
    id.r = N-Nx2+1:N;
    id.b = 1:Nx2:N;
    id.t = Nx2:Nx2:N;
    id.lb = 1;
    id.lt = Nx2;
    id.rb = N-Nx2+1;
    id.rt = N;
    id.corner = [id.lb, id.lt, id.rb, id.rt];

    id.inner = setdiff(id.all, [id.l, id.r, id.b, id.t]);
    id.l = setdiff(id.l, id.corner);
    id.r = setdiff(id.r, id.corner);
    id.b = setdiff(id.b, id.corner);
    id.t = setdiff(id.t, id.corner);
    id.boundary = unique([id.l, id.r, id.b, id.t]);
end