function id = GetIndex(Nx1, Nx2)
    N = Nx1*Nx2;
    id.l = 1:Nx2;
    id.r = N-Nx2+1:N;
    id.b = 1:Nx2:N;
    id.t = Nx2:Nx2:N;
    id.lb = 1;
    id.lt = Nx2;
    id.rb = N-Nx2+1;
    id.rt = N;
    id.corner = [id.lb, id.lt, id.rb, id.rt];

    id.inner = setdiff(1:N, [id.l, id.r, id.b, id.t]);
    id.l = setdiff(id.l, id.corner);
    id.r = setdiff(id.r, id.corner);
    id.b = setdiff(id.b, id.corner);
    id.t = setdiff(id.t, id.corner);
    id.boundary = unique([id.l, id.r, id.b, id.t]);

    %{
    id.l_sub1 = sub2ind([2*N,2*N], id.l, id.l);
    id.l_sub2 = sub2ind([2*N,2*N], id.l, N+id.l);
    id.l_sub3 = sub2ind([2*N,2*N], N+id.l, id.l);
    id.l_sub4 = sub2ind([2*N,2*N], N+id.l, N+id.l);
    id.l_sub5 = sub2ind([2*N,2*N], N+id.l, 1*Nx2+id.l);
    id.l_sub6 = sub2ind([2*N,2*N], N+id.l, N+1*Nx2+id.l);
    id.l_sub7 = sub2ind([2*N,2*N], N+id.l, 2*Nx2+id.l);
    id.l_sub8 = sub2ind([2*N,2*N], N+id.l, N+2*Nx2+id.l);

    id.r_sub1 = sub2ind([2*N,2*N], id.r, id.r);
    id.r_sub2 = sub2ind([2*N,2*N], id.r, N+id.r);
    id.r_sub3 = sub2ind([2*N,2*N], N+id.r, id.r);
    id.r_sub4 = sub2ind([2*N,2*N], N+id.r, N+id.r);
    id.r_sub5 = sub2ind([2*N,2*N], N+id.r, -1*Nx2+id.r);
    id.r_sub6 = sub2ind([2*N,2*N], N+id.r, N-1*Nx2+id.r);
    id.r_sub7 = sub2ind([2*N,2*N], N+id.r, -2*Nx2+id.r);
    id.r_sub8 = sub2ind([2*N,2*N], N+id.r, N-2*Nx2+id.r);

    id.b_sub1 = sub2ind([2*N,2*N], id.b, id.b);
    id.b_sub2 = sub2ind([2*N,2*N], id.b, N+id.b);
    id.b_sub3 = sub2ind([2*N,2*N], N+id.b, id.b);
    id.b_sub4 = sub2ind([2*N,2*N], N+id.b, N+id.b);
    id.b_sub5 = sub2ind([2*N,2*N], N+id.b, 1+id.b);
    id.b_sub6 = sub2ind([2*N,2*N], N+id.b, N+1+id.b);
    id.b_sub7 = sub2ind([2*N,2*N], N+id.b, 2+id.b);
    id.b_sub8 = sub2ind([2*N,2*N], N+id.b, N+2+id.b);

    id.t_sub1 = sub2ind([2*N,2*N], id.t, id.t);
    id.t_sub2 = sub2ind([2*N,2*N], id.t, N+id.t);
    id.t_sub3 = sub2ind([2*N,2*N], N+id.t, id.t);
    id.t_sub4 = sub2ind([2*N,2*N], N+id.t, N+id.t);
    id.t_sub5 = sub2ind([2*N,2*N], N+id.t, -1+id.t);
    id.t_sub6 = sub2ind([2*N,2*N], N+id.t, N-1+id.t);
    id.t_sub7 = sub2ind([2*N,2*N], N+id.t, -2+id.t);
    id.t_sub8 = sub2ind([2*N,2*N], N+id.t, N-2+id.t);
    %}
end