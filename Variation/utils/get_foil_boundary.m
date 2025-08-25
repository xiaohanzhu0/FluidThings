function [b1,b2,b3,b4,aux] = get_foil_boundary(append_trail)
addpath('~/Files/data/Mesh_generation/Airfoil');

m = load('metric_airfoil_bodyfitted.mat');
foil = readmatrix('A-airfoil.txt', 'NumHeaderLines', 1);
foil = foil(:,1:2);

%% Building interior boundary
trail_factor = 1.05;

N = 50;
j = 1:N;

%x = 0.5 - 0.5*cos(pi*(N-j)/(N-1));
x = 1 - cos(pi*(N-j)/(N-1)/2);
x = [x(1), x(3:end)];

[o,i_head] = min(foil(:,1));
foil_up = foil(i_head:end,:);
foil_down = foil(1:i_head,:);

y_up = interp1(foil_up(:,1),foil_up(:,2),x,'spline');
y_down = interp1(foil_down(:,1),foil_down(:,2),x,'spline');

counter = 0;

if append_trail
    while x(1) <= 30
        x = [x(1)+(x(1)-x(2))*trail_factor, x];
        y_up = [y_up(1), y_up];
        y_down = [y_down(1), y_down];
        counter = counter + 1;
    end
end

x = [x(1:end-1), flip(x)];
y = [y_up(1:end-1), flip(y_down)];
b1 = [x', y'];
b3 = b1;

%%
side_init_step = -0.1;
side_factor = 1.0;

b2_y = exp_stretch(b1(end,2), -27, 100, 1.05);
b2_x = b1(end,1)*ones(size(b2_y));
b2 = [b2_x; b2_y]';

b4_y = exp_stretch(b1(end,2), 27, 100, 1.05);
b4_x = b1(end,1)*ones(size(b4_y));
b4 = [b4_x; b4_y]';

%%

b3(1:counter, 2) = b4(end, 2);
b3(end-counter+1:end, 2) = b2(end, 2);
r0 = (b3(counter+1, :) + b3(end-counter, :))/2;
r = abs(b4(end, 2) - r0(2));
N = length(b3(counter+1:end-counter,:));

%t = GetBoundaryTangent(b1(:,1), b1(:,2), 1);
%n = [t(2,:); -t(1,:)];
%n(1,:) = min(n(1,:), b1(:,1)');
%b3 = b1 + n'*b4(end, 2);

theta = 0.5 - 0.5* cos(pi*((N-1)-((N-1):-1:0))/(N-1));
theta = 0.5 - 0.5* cos(pi*theta);
b3(counter+1:end-counter,1) = r*cos(pi/2 + pi*theta) + r0(1);
b3(counter+1:end-counter,2) = r*sin(pi/2 + pi*theta) + r0(2);

b3(counter+1:end-counter,1) = r*cos(pi/2 + pi*(0:(N-1))'/(N-1)) + r0(1);
b3(counter+1:end-counter,2) = r*sin(pi/2 + pi*(0:(N-1))'/(N-1)) + r0(2);

auxN = b3(ceil(size(b3,1)/2),:);
aux1 = b1(ceil(size(b3,1)/2),:);

aux_x = exp_stretch(aux1(1), auxN(1), 100, 1.05);
aux_y = linspace(aux1(2),auxN(2),100);
aux = [aux_x; aux_y]';

%%
scatter(b1(:,1),b1(:,2)); hold on
scatter(b2(:,1),b2(:,2)); hold on
scatter(b3(:,1),b3(:,2)); hold on
scatter(b4(:,1),b4(:,2)); hold on
scatter(aux(:,1),aux(:,2)); hold on
end

function a = exp_stretch(a1, aN, N, alpha)
    n = 1:N;
    a = ((alpha^N-alpha.^n)*a1 + (alpha.^n-alpha)*aN) / (alpha^N - alpha);
end