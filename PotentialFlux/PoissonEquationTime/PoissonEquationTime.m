#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path

clear; clc; close all

Nx = 512;   % or 1024
dx = 1e-4;
Lx = Nx * dx;              % adjust Lx to match dx
x = linspace(0, Lx - dx, Nx);

epsilon0 = 8.854e-12;         % vacuum permittivity
mu = 1e-12;                   % mobility
D  = 5e-6;                    % diffusion coefficient
dt = 1e-4;
%dt = 1e-6;
Nt = 10000;                    % total steps
plot_interval = 100;           % plot frequency
Faraday = 96485.3329;          % C/mol

% -------------------------------------------------------------------------
% Positive carrier densities (Fe2+ and Cl-)
% -------------------------------------------------------------------------
x0_Fe = 0.3 * Lx;     % Fe2+ center
x0_Cl = 0.7 * Lx;     % Cl- center
sigma = 0.05 * Lx;    % width

c_Fe = exp(-((x - x0_Fe).^2) / (2*sigma^2));   % Fe2+ concentration (positive)
c_Cl = exp(-((x - x0_Cl).^2) / (2*sigma^2));   % Cl- concentration (positive)

z_Fe =  2;       % charge number for Fe2+
z_Cl = -2;       % charge number for Cl-


% -------------------------------------------------------------------------
% FFT wavenumbers
% -------------------------------------------------------------------------
k  = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];
k2 = k.^2;  
k2(1) = 1;

% -------------------------------------------------------------------------
% Storage
% -------------------------------------------------------------------------

rho_snapshots = [];
phi_snapshots = [];
time_labels   = [];
convergence   = zeros(1, Nt);
energy        = zeros(1, Nt);

dt_cfl = dx^2 / (2*D);
myprint("CFL : ", dt_cfl, " < then dt : ", dt);

rhoPlot = SnapshotPlotter('Charge_Evolution.png');

% -------------------------------------------------------------------------
% TIME LOOP
% -------------------------------------------------------------------------
for t = 1:Nt
    % ---- Compute total charge density ----
    rho = Faraday * (z_Fe * c_Fe + z_Cl * c_Cl);
    rho = rho - mean(rho);     % enforce neutrality

    rho(1)   = rho(end-1);
    rho(end) = rho(2);

    % ---- Solve Poisson: φ_xx = -ρ/ε0 ----
    rho_hat = fft(rho);
    phi_hat = rho_hat ./ (epsilon0 * k2 * 1e7);
    phi_hat(1) = 0;
    phi = real(ifft(phi_hat));

    % ---- Electric field: E = -φ_x ----
    E_hat = 1i * k .* phi_hat;
    E = -real(ifft(E_hat));

    % ---- Drift + diffusion for each species ----
    grad_Fe = real(ifft(1i * k .* fft(c_Fe)));
    grad_Cl = real(ifft(1i * k .* fft(c_Cl)));

    J_Fe = -D * grad_Fe + z_Fe * mu * c_Fe .* E;
    J_Cl = -D * grad_Cl + z_Cl * mu * c_Cl .* E;

    % ---- Continuity: c_t = -dJ/dx ----
    c_Fe_dot = -real(ifft(1i * k .* fft(J_Fe)));
    c_Cl_dot = -real(ifft(1i * k .* fft(J_Cl)));

    % ---- Ensure positive evolution ----
    %c_Fe_dot = abs(c_Fe_dot);
    %c_Cl_dot = abs(c_Cl_dot);

    % ---- Update concentrations ----
    c_Fe = c_Fe + dt * c_Fe_dot;
    c_Cl = c_Cl + dt * c_Cl_dot;

    % ---- Keep densities positive ----
    c_Fe = max(c_Fe, 0);
    c_Cl = max(c_Cl, 0);

    
    % ---- Periodic boundary wrapping (ghost cells) ----
    c_Fe(1)   = c_Fe(end-1);
    c_Fe(end) = c_Fe(2);
    c_Cl(1)   = c_Cl(end-1);
    c_Cl(end) = c_Cl(2);

    % ---- Diagnostics ----
    convergence(t) = sqrt(sum((rho).^2) / Nx);
    energy(t) = 0.5 * epsilon0 * sum(E.^2) * dx;

    if mod(t, plot_interval) == 0
        fprintf("Step %d:  <rho>=%.3e  Energy=%.3e\n", t, mean(rho), energy(t));


         rhoPlot.add_and_plot(x, rho, phi,  'x', '\rho(x)', '\phi(x)');
        % ---- Store snapshots ----
    end

    % ---- Stop if steady ----
    if convergence(t) < 1e-10
        fprintf('Equilibrium reached at step %d (t=%.3g s)\n', t, t*dt);
        break;
    end
end

fprintf('✓ Simulation complete. Final plots saved.\n');

