clear
Nx1 = 80;
Nx2 = 80;

e1 = ones(1,Nx1);
e2 = ones(1,Nx2);

N_line = floor(Nx1/3);
N_circ = Nx1 - 2*N_line;
bottom1 = [linspace(0,1,Nx1); 0*e1];
bottom2 = 1/6*[0.5*6+cos(linspace(pi,0,N_circ)); sin(linspace(pi,0,N_circ))];
bottom = [bottom1(:,1:N_line), bottom2, bottom1(:,N_line+N_circ+1:end)];
right = [e2; linspace(0,1,Nx2)];
top = [linspace(0,1,Nx1); e1];
left = [0*e2; linspace(0,2,Nx2)];

dx = diff(bottom(1,:));
dy = diff(bottom(2,:));
ds = sqrt(dx.^2 + dy.^2);
s = [0, cumsum(ds)];
s = s / s(end);

x = bottom(1,:);
y = bottom(2,:);
% Given target point (x_star, y_star)
x_star = 0.49; % your x* value
y_star = 0.2; % your y* value

N = length(x);
minDist = inf;
bestIdx = 0;
best_r = 0;


%%
clear
N = 41;
x1 = @(s) s;
x2 = @(s) -10*sin(2*pi*s);
F = @(phi) sqrt(1+(2*pi*10*cos(2*pi*phi)).^2);

phi_solution = BoundaryMesh(F, N);
plot(x1(phi_solution), x2(phi_solution),'o'); hold on
plot(linspace(0,1,N), x2(linspace(0,1,N)));
xlim([0,1]);
ylim([-10,10]);

%%
x1 = linspace(0,1,N);
x2 = linspace(0,1,N);
surf()
%%
clear
% Define the function F(phi)
F = @(phi) sqrt(1+(2*pi*4*cos(2*pi*phi)).^2); % Example function; replace with your F(phi)

% Compute the normalization constant c
c = integral(F, 0, 1);

% Create a fine grid for phi
phi_grid = linspace(0, 1, 1000);

% Evaluate F on the grid (ensure F is vectorized)
F_values = F(phi_grid);

% Compute cumulative integral of F(phi) using trapezoidal rule
integral_F = cumtrapz(phi_grid, F_values);

% Compute corresponding xi values
xi_grid = integral_F / c;

% Interpolate to get phi as a function of xi
xi_desired = linspace(0, 1, 40); % Desired xi points
phi_solution = interp1(xi_grid, phi_grid, xi_desired, 'pchip');

% Plot the result
plot(xi_desired, phi_solution);
xlabel('\xi');
ylabel('\phi(\xi)');
title('Solution of Equation 9.10');

plot(phi_grid, -4*sin(2*pi*phi_grid)); hold on
scatter(phi_solution, -4*sin(2*pi*phi_solution));


