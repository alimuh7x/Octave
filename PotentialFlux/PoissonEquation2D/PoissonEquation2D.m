#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
set(0, "defaultfigureposition", [100, 100, 800, 650]);

addpath(genpath("./../../src"));  % add src folder to path



% -------------------------------------------------------------------------
% 2D Poisson equation using FFT: ∇²φ = -ρ/ε₀
% -------------------------------------------------------------------------
clear; clc;

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
mu = 5e-10;
Nx = 512; Ny = 512;
Lx = 1.0; Ly = 1.0;

dx = Lx / Nx; 
dy = Ly / Ny;
epsilon0 = 8.854e-12;

x = linspace(0, Lx - dx, Nx);
y = linspace(0, Ly - dy, Ny);
[X,Y] = meshgrid(x, y);

% -------------------------------------------------------------------------
% Charge density (initial)
% -------------------------------------------------------------------------
rho = sin(2*pi*X) .* cos(2*pi*Y);

% -------------------------------------------------------------------------
% Solve Poisson equation in Fourier space
% -------------------------------------------------------------------------
rho_hat = fft2(rho);

kx = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];
ky = (2*pi/Ly) * [0:(Ny/2-1) -Ny/2:-1];
[KX, KY] = meshgrid(kx, ky);

k2 = KX.^2 + KY.^2;
k2(1,1) = 1;       % avoid division by zero

phi_hat = rho_hat ./ (epsilon0 * k2);
phi_hat(1,1) = 0;  % zero-mean potential
phi = real(ifft2(phi_hat));

% -------------------------------------------------------------------------
% Electric field
% -------------------------------------------------------------------------
[phi_x, phi_y] = gradient(phi, dx, dy);
Ex = -phi_x;
Ey = -phi_y;
E_mag = sqrt(Ex.^2 + Ey.^2);

% -------------------------------------------------------------------------
% Separate positive and negative charge densities
% -------------------------------------------------------------------------
rho_pos = max(rho, 0);   % positive charge (e.g., Fe2+)
rho_neg = min(rho, 0);   % negative charge (e.g., Cl-)

z_pos = +1;
z_neg = -1;
mu_pos = mu;
mu_neg = mu;

% -------------------------------------------------------------------------
% Fluxes for each charge type
% -------------------------------------------------------------------------
Jx_pos = mu_pos * rho_pos .* Ex * z_pos;
Jy_pos = mu_pos * rho_pos .* Ey * z_pos;

Jx_neg = mu_neg * rho_neg .* Ex * z_neg;
Jy_neg = mu_neg * rho_neg .* Ey * z_neg;

% Total flux
Jx_total = Jx_pos + Jx_neg;
Jy_total = Jy_pos + Jy_neg;

% -------------------------------------------------------------------------
% Rate of charge change (divergence of flux)
% -------------------------------------------------------------------------
[divJx, divJy] = gradient(Jx_total, dx, dy);
rho_dot = -(divJx + divJy);
rho_dot = rho_dot - mean(rho_dot(:));

% Update charge density
dt = 0.01;
new_rho = rho + dt * rho_dot;

% -------------------------------------------------------------------------
%  NOTE: Visualization using FieldPlot class (with electric field arrows)
% -------------------------------------------------------------------------

% Default arrow downsampling for clarity
skip = 15;

% --- Initial charge density ---
p = FieldPlot("rho_initial.png");
p.add_Image(x, y, rho, "Initial \\rho");
p.save();

% --- Positive charge density ---
p = FieldPlot("rho_positive.png");
p.add_Image(x, y, rho_pos, "Positive charge density \\rho^+", ColorMaps.jet);
p.save();

% --- Negative charge density ---
p = FieldPlot("rho_negative.png");
p.add_Image(x, y, rho_neg, "Negative charge density \\rho^-", ColorMaps.jet);
p.save();

% --- Electric potential with field arrows ---
p = FieldPlot("potential_phi_with_E.png");
p.add_Image(x, y, phi, "Electric potential \\phi with field arrows", ColorMaps.jet);
p.add_Arrows(X, Y, Ex, Ey, skip, 0.03, colors.black);
p.save();

% --- Electric field magnitude with arrows ---
p = FieldPlot("E_magnitude_with_arrows.png");
p.add_Image(x, y, E_mag, "Electric field magnitude |E|", ColorMaps.jet);
p.add_Arrows(X, Y, Ex, Ey, skip, 0.03, colors.black);
p.save();

% --- Updated charge density ---
p = FieldPlot("rho_updated.png");
p.add_Image(x, y, new_rho, "Updated charge density", ColorMaps.jet);
p.save();


