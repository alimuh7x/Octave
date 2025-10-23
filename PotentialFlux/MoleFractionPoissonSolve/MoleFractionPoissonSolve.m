#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path

clear; clc; close all;

% -------------------------------------------------------------------------
% Domain setup
% -------------------------------------------------------------------------
Nx = 200;
Lx = 1.0e-7;              % 100 nm domain
dx = Lx / (Nx - 1);
x  = linspace(0, Lx, Nx);

% -------------------------------------------------------------------------
% Physical constants (Cui et al. 2022)
% -------------------------------------------------------------------------
epsilon0 = 8.854e-12;     % vacuum permittivity (F/m)
F = 96485;                 % Faraday constant (C/mol)
R = 8.314;                 % gas constant (J/mol/K)
T = 298;                   % temperature (K)

% Diffusion coefficients (m^2/s)
D_Fe = 5e-10;
D_Cl = 5e-10;

% Mobilities (Nernst–Einstein relation)
mu_Fe = D_Fe * F / (R * T);
mu_Cl = D_Cl * F / (R * T);

% Atomic and saturation concentrations (mol/m³)
c_solid = 1.43e5;          % atomic density of Fe
c_sat   = 5.1e3;           % solubility limit (Fe²⁺ in electrolyte)

% Normalized mole fractions (dimensionless)
X_Fe_sat = c_sat / c_solid;    % = 0.036
X_Cl_sat = X_Fe_sat;           % same for neutrality

% Time settings
dt = 1e-3;
nsteps = 1;

% -------------------------------------------------------------------------
% Initial mole fraction profiles (dimensionless)
% -------------------------------------------------------------------------
X_Fe = X_Fe_sat * exp(-((x - 0.3*Lx).^2) / (0.002*Lx^2));
X_Cl = X_Cl_sat * exp(-((x - 0.7*Lx).^2) / (0.002*Lx^2));

% Normalize (max = X_sat)
X_Fe = X_Fe / max(X_Fe) * X_Fe_sat;
X_Cl = X_Cl / max(X_Cl) * X_Cl_sat;

% Convert to actual concentration (mol/m³)
c_Fe = X_Fe * c_solid;
c_Cl = X_Cl * c_solid;

z_Fe =  2;  % Fe²⁺
z_Cl = -1;  % Cl⁻

% -------------------------------------------------------------------------
% Time evolution (single iteration)
% -------------------------------------------------------------------------
for step = 1:nsteps

    % --- Charge density (C/m³)
    rho = F * (z_Fe * c_Fe + z_Cl * c_Cl);

    % --- Solve Poisson equation via FFT
    rho_hat = fft(rho);
    k = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];
    k2 = k.^2;  k2(1) = 1;
    phi_hat = -rho_hat ./ (epsilon0 * k2);
    phi_hat(1) = 0;
    phi = real(ifft(phi_hat));

    % --- Electric field
    E = -gradient(phi, dx);

    % --- Ionic fluxes (drift only)
    J_Fe =  mu_Fe * z_Fe * F * c_Fe .* E;
    J_Cl =  mu_Cl * z_Cl * F * c_Cl .* E;

    % --- Update concentrations (conservation)
    c_Fe = c_Fe - dt * gradient(J_Fe, dx) / (F * z_Fe);
    c_Cl = c_Cl - dt * gradient(J_Cl, dx) / (F * z_Cl);

    % --- Clamp to physical range
    c_Fe(c_Fe < 0) = 0;
    c_Cl(c_Cl < 0) = 0;
endfor

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------
figure(1);
subplot(3,1,1)
plot(x*1e9, c_Fe, 'r', x*1e9, c_Cl, 'b', 'LineWidth', 1.5);
xlabel('x (nm)'); ylabel('c_i (mol/m^3)');
legend('Fe^{2+}', 'Cl^-');
title('Final Concentrations (Cui et al. 2022 values)');
grid on;

subplot(3,1,2)
plot(x*1e9, phi, 'k', 'LineWidth', 1.5);
xlabel('x (nm)'); ylabel('\phi (V)');
title('Electric Potential');
grid on;

subplot(3,1,3)
plot(x*1e9, E, 'm', 'LineWidth', 1.5);
xlabel('x (nm)'); ylabel('E (V/m)');
title('Electric Field');
grid on;
print("ElectricField_Vectors.png", "-dpng");

figure(2);
plot(x*1e9, rho, 'r', 'LineWidth', 1.5);
xlabel('x (nm)'); ylabel('\rho (C/m^3)');
title('Charge Density');
grid on;
print("Density.png", "-dpng");

disp('Simulation completed using mole fractions from Cui et al. (2022).');

