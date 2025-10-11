#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
% ---------------------------------------------------------------------------
% NOTE: Solving Numerically using finite difference method (Gauss-Seidel)
% ---------------------------------------------------------------------------

% Poisson equation in 1D: d²φ/dx² = -ρ(x)/ε₀
clear; clc;

% Parameters
Nx = 200;          % grid points
Lx = 1.0;          % length (m)
dx = Lx / (Nx - 1);
epsilon0 = 8.854e-12;

% Grid
x = linspace(0, Lx, Nx);

% Charge density: cosine wave
rho = sin(2 * pi * x);

% Potential initialization (Dirichlet BC: φ=0 at both ends)
phi = zeros(1, Nx);

% Iteration parameters
max_iter = 5000;
tol = 1e-6;

% Gauss–Seidel iteration
for iter = 1:max_iter
    phi_old = phi;
    for i = 2:Nx-1
        phi(i) = 0.5 * (phi(i+1) + phi(i-1) + dx^2 * rho(i) / epsilon0);
    end

% Print iteration info
    if mod(iter, 500) == 0
        printf("Iteration %d: max change = %.6e\n", iter, max(abs(phi - phi_old)));
    end

    if max(abs(phi - phi_old)) < tol
        printf("Converged in %d iterations\n", iter);
        break;
    end
end

% Analytical solution
phi_exact = sin(2*pi*x) ./ (4*pi^2*epsilon0);

% Plot results
subplot(3,1,1)
plot(x, rho, 'r', 'LineWidth', 1.5)
xlabel('x'); ylabel('\rho(x)')
title('Charge density \rho(x) = sin(2\pi x)')
grid on

subplot(3,1,2)
plot(x, phi, 'b', 'LineWidth', 1.5)
xlabel('x'); ylabel('\phi_{num}(x)')
title('Numerical potential (Gauss–Seidel)')
grid on

subplot(3,1,3)
plot(x, phi_exact, 'k--', 'LineWidth', 1.5)
xlabel('x'); ylabel('\phi_{ana}(x)')
title('Analytical potential \phi(x) = sin(2\pi x)/(4\pi^2 \epsilon_0)')
grid on

print("Numerical.png", "-dpng");


% ---------------------------------------------------------------------------
% NOTE: Solving Numerically using finite difference method (SOR)
% ---------------------------------------------------------------------------

% Poisson equation in 1D using SOR: d²φ/dx² = -ρ(x)/ε₀
clear; clc;

% Parameters
Nx = 200;
Lx = 1.0;
dx = Lx / (Nx - 1);
epsilon0 = 8.854e-12;
omega = 1.9;       % relaxation factor (1 < ω < 2 for acceleration)

% Grid
x = linspace(0, Lx, Nx);

% Charge density
rho = cos(2 * pi * x);

% Initialize potential (Dirichlet BC)
phi = zeros(1, Nx);

% Iteration parameters
max_iter = 10000;
tol = 1e-9;

for iter = 1:max_iter
    maxdiff = 0;
    for i = 2:Nx-1
        phi_new = 0.5 * (phi(i+1) + phi(i-1) + dx^2 * rho(i) / epsilon0);
        phi(i) = (1 - omega) * phi(i) + omega * phi_new;
        maxdiff = max(maxdiff, abs(phi_new - phi(i)));
    end
% Print iteration info
    if mod(iter, 500) == 0
        printf("Iteration %d: max change = %.6e\n", iter, maxdiff);
    end
    if maxdiff < tol
        printf("Converged in %d iterations\n", iter);
        break;
    end
end

% Plot results
subplot(2,1,1)
plot(x, rho, 'r', 'LineWidth', 1.5)
xlabel('x'); ylabel('\rho(x)')
title('Charge density \rho(x) = cos(2πx)')
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
grid on

subplot(2,1,2)
plot(x, phi, 'b', 'LineWidth', 1.5)
xlabel('x'); ylabel('\phi(x)')
title('Potential \phi(x) using SOR (\omega = 1.9)')
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
grid on

print ("Numerical_SOR.png", "-dpng");  % saves to PNG


