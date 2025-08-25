clear
Nx1 = 5;
Nx2 = 5;
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);
[x1_temp, x2_temp] = meshgrid(s1,s2);
x1 = x1_temp*cos(pi/4) - x2_temp*sin(pi/4);
x2 = x1_temp*sin(pi/4) + x2_temp*cos(pi/4);
x1 = x1+1/sqrt(2);


M = GetM(x1, x2, 1, 0.1);
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
%%
dM11dx1_ex = M.dM11dx1;
dM11dx2_ex = M.dM11dx2;

[dM11dx1, dM11dx2] = metric_grad(M.M11, x1, x2);
[dM11dx1_notgood, dM11dx2_notgood] = DCentralUneven(M.M11, M.M11, x1, x2);

norm(dM11dx2(2:end-1,2:end-1)-dM11dx2_ex(2:end-1,2:end-1), 'fro') / sqrt(Nx1*Nx2)
norm(dM11dx2_notgood(2:end-1,2:end-1)-dM11dx2_ex(2:end-1,2:end-1), 'fro') / sqrt(Nx1*Nx2)