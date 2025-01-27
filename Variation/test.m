clear
N = 41;
x1 = @(s) s;
x2 = @(s) -10*sin(2*pi*s);
F = @(phi) sqrt(1+(2*pi*10*cos(2*pi*phi)).^2);

phi_solution = SymBoundaryMesh(F, N);
plot(x1(phi_solution), x2(phi_solution),'o'); hold on
plot(linspace(0,1,N), x2(linspace(0,1,N)));
xlim([0,20]);
ylim([-10,10]);

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


