clear
Nx1_samp = 160;
Nx2_samp = 160;
s1 = linspace(0, 1, Nx1_samp);
s2 = linspace(0, 1, Nx2_samp);
[x1_samp, x2_samp] = meshgrid(s1,s2);

M_samp = GetM(x1_samp, x2_samp, 2, 1);
plot(x1_samp, x2_samp, 'k'); hold on; plot(x1_samp', x2_samp', 'k');

x1 = linspace(0, 1, Nx1_samp/4);
x2 = linspace(0, 1, Nx2_samp/4);
[x1, x2] = meshgrid(x1,x2);

h1 = 1/Nx1_samp/1000000;
h2 = 1/Nx2_samp/1000000;

x1_p = x1+h1; x1_n = x1-h1;
x2_p = x2+h2; x2_n = x2-h2;

dM11dx1_p = griddata(x1_samp,x2_samp,M_samp.M11,x1_p,x2,"cubic");
dM11dx1_n = griddata(x1_samp,x2_samp,M_samp.M11,x1_n,x2,"cubic");

dM11dx2_p = griddata(x1_samp,x2_samp,M_samp.M11,x1,x2_p,"cubic");
dM11dx2_n = griddata(x1_samp,x2_samp,M_samp.M11,x1,x2_n,"cubic");

dM11dx1 = (dM11dx1_p - dM11dx1_n) ./ (x1_p - x1_n);
dM11dx2 = (dM11dx2_p - dM11dx2_n) ./ (x2_p - x2_n);