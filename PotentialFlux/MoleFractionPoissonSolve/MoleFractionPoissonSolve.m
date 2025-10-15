#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path


clear; clc;

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
Nx = 200;
Lx = 1.0;              % domain length (m)
dx = Lx / (Nx - 1);
x  = linspace(0, Lx, Nx);

epsilon0 = 8.854e-12;  % vacuum permittivity (F/m)
F = 96485;             % Faraday constant (C/mol)

mu_Fe = 5e-10;         % mobility of Fe2+
mu_Cl = 5e-10;         % mobility of Cl-
dt = 0.001;            % time step (s)
nsteps = 1;          % number of time steps

% -------------------------------------------------------------------------
% Initial mole fractions → concentrations (mol/m³)
% -------------------------------------------------------------------------
X_Fe = 0.5 * exp(-((x - 0.3).^2) / 0.002);
X_Cl = 0.5 * exp(-((x - 0.7).^2) / 0.002);

X_Fe = X_Fe / max(X_Fe);
X_Cl = X_Cl / max(X_Cl);
c_total = 1000;         % mol/m³ (≈1 M electrolyte)
c_Fe = X_Fe * c_total;
c_Cl = X_Cl * c_total;

z_Fe =  2;  % Fe2+
z_Cl = -1;  % Cl-

% -------------------------------------------------------------------------
% Time evolution loop
% -------------------------------------------------------------------------
for step = 1:nsteps
    % --- Total charge density from concentrations ---
    rho = F * (z_Fe * c_Fe + z_Cl * c_Cl);

    % --- Solve Poisson equation using FFT ---
    rho_hat = fft(rho);
    k = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];
    k2 = k.^2;  k2(1) = 1;
    phi_hat = rho_hat ./ (epsilon0 * k2);
    phi_hat(1) = 0;
    phi = real(ifft(phi_hat));

    % --- Electric field ---
    E = -gradient(phi, dx);

    % --- Ionic fluxes (drift only) ---
    J_Fe =  mu_Fe * z_Fe * F * c_Fe .* E;
    J_Cl =  mu_Cl * z_Cl * F * c_Cl .* E;

    % --- Time evolution of concentrations (ċ = -∇·J/Fz) ---
    c_Fe = c_Fe - dt * gradient(J_Fe, dx) / (F * z_Fe);
    c_Cl = c_Cl - dt * gradient(J_Cl, dx) / (F * z_Cl);

    % --- Keep within physical range ---
    c_Fe(c_Fe < 0) = 0;
    c_Cl(c_Cl < 0) = 0;

    if mod(step, 100) == 0
        fprintf('Step %d completed\n', step);
    end
end

% -------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------
figure;

subplot(3,1,1)
plot(x, c_Fe, 'r', 'LineWidth', 1.5); hold on;
plot(x, c_Cl, 'b', 'LineWidth', 1.5);
xlabel('x'); ylabel('c_i (mol/m^3)');
legend('Fe^{2+}','Cl^-');
title('Final Concentrations');
grid on;

subplot(3,1,2)
plot(x, phi, 'k', 'LineWidth', 1.5);
xlabel('x'); ylabel('\phi (V)');
title('Electric Potential');
grid on;

subplot(3,1,3)
plot(x, E, 'm', 'LineWidth', 1.5);
xlabel('x'); ylabel('E (V/m)');
title('Electric Field');
grid on;

print("ElectricField_Vectors.png", "-dpng");

figure;


plot(x, rho, 'r', 'LineWidth', 1.5);
xlabel('x'); ylabel('c_i (mol/m^3)');
title('Final Concentrations');
grid on;

print("Density.png", "-dpng");

