#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path

clear; clc;

% --- Diffusion rate and timestep estimator ---

% Input parameters
D  = 1e-10;       % diffusivity (m^2/s)
dx = 1e-6;        % grid spacing (m)

% --- Recommended c_dot magnitude (1/s)
c_dot_est = D / dx^2;

% --- Recommended stable timestep (s)
dt_recommended = dx^2 / (10 * D);
fprintf('=================================================================\n');
fprintf('  c_dot * dt     <  0.1 \n');
fprintf('  D * dt         <  1/6 dx^2  \n');

% --- Display results ---
fprintf('=================================================================\n');
fprintf('  Given Parameters:\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('  D  = %.1e m^2/s \n  dx = %.1e m\n', D, dx);
fprintf('=================================================================\n');
fprintf('  Calculated Parameters:\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('  Recommended c_dot  ≈   %.3f 1/s\n', c_dot_est);
fprintf('  Recommended dt     <   %.3f s\n', dt_recommended);
fprintf('  c_dot * dt         =   %.3f \n', c_dot_est * dt_recommended);
fprintf('=================================================================\n');

% --- Parameters ---
Nx = 400;
L  = 4.0;
dx = L / (Nx - 1);
x  = linspace(-L/2, L/2, Nx);

dt = 1e-4;
nsteps = 5000;

% Phase-field and material parameters
w = 0.2;                 % interface width
A = 1;

Ceq_alpha = 0.2;
Ceq_beta  = 1.0;

C_alpha = 0.0;
C_beta  = 1.0;

M = 1e-14;                 % constant mobility

% --- Phase-field profile ---
phi = 0.5 * (1 + tanh(x / w));

phi_alpha = 1 - phi;
phi_beta  = phi;

% --- Total composition ---
c   = phi_alpha .* C_alpha + phi_beta .* C_beta;
Ceq   = phi_alpha .* Ceq_alpha + phi_beta .* Ceq_beta;

mu  = A * (c - Ceq);

mu2 = A .* phi_alpha .* (c - Ceq_alpha) + A .* phi_beta .* (c - Ceq_beta);


SinglePlot(x, c, 'x', 'c(x)', 'Initial_Concentration.png');
SinglePlot(x, Ceq, 'x', 'c(x)', 'Initial_Equilibrium_Concentration.png');
SinglePlot(x, phi, 'x', 'phi', 'Phase_Field.png');
SinglePlot(x, mu, 'x', 'mu', 'InitialPotential.png');
SinglePlot(x, mu2, 'x', 'mu', 'InitialPotential2.png');

rhoPlot = SnapshotPlotter('Charge_Evolution.png');

for step = 1:nsteps
    % Chemical potential (single field)
    mu = A * (c - ( phi_alpha * Ceq_alpha + phi_alpha * Ceq_beta));

    % Laplacian with periodic BC
    mu_xx = (circshift(mu,-1) - 2*mu + circshift(mu,1)) / dx^2;

    % Update concentration
    dc = M * mu_xx;
    c = c + dt * dc;

    % Visualization
    if mod(step, 500) == 0
        printf('Step %d / %d\n', step, nsteps);
        %rhoPlot.add_and_plot(x, mu, c,  'x', '\mu(x)', 'c(x)');
    end
end

