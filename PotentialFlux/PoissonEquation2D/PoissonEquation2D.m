#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
set(0, "defaultfigureposition", [100, 100, 800, 650]);


addpath(genpath("./../../src"));  % add src folder to path


% -------------------------------------------------------------------------
% 2D Poisson equation using FFT: ∇²φ = -ρ/ε₀
% -------------------------------------------------------------------------
clear; clc;

mu = 5e-10;
% Parameters
Nx = 128; Ny = 128;
Lx = 1.0; Ly = 1.0;

dx = Lx / Nx; 
dy = Ly / Ny;

epsilon0 = 8.854e-12;

x = linspace(0, Lx - dx, Nx);
y = linspace(0, Ly - dy, Ny);

[X, Y] = meshgrid(x, y);

% Charge density
rho = sin(2*pi*X) .* cos(2*pi*Y);

% FFT of charge
rho_hat = fft2(rho);

% Wavenumbers
kx = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];
ky = (2*pi/Ly) * [0:(Ny/2-1) -Ny/2:-1];

[KX, KY] = meshgrid(kx, ky);

k2 = KX.^2 + KY.^2;
k2(1,1) = 1;       % avoid division by zero

% Solve in Fourier space
phi_hat = rho_hat ./ (epsilon0 * k2);
phi_hat(1,1) = 0;  % zero-mean potential
phi = real(ifft2(phi_hat));

% -------------------------------------------------------------------------
% Electric field and flux
% -------------------------------------------------------------------------

[phi_x, phi_y] = gradient(phi, dx, dy);

Ex = -phi_x; 
Ey = -phi_y;
E_mag = sqrt(Ex.^2 + Ey.^2);


Jx = mu * rho .* Ex;
Jy = mu * rho .* Ey;

% -------------------------------------------------------------------------
% Rate of charge change
% -------------------------------------------------------------------------

[divJx, divJy] = gradient(Jx, dx, dy);

rho_dot = -(divJx + divJy);
rho_dot = rho_dot - mean(rho_dot(:));   % remove drift

new_rho = rho + 0.01 * rho_dot;

% -------------------------------------------------------------------------
% Plot: Charge Density
% -------------------------------------------------------------------------

figure;
imagesc(x, y, rho);
axis image; colorbar;
title('\rho(x,y) = sin(2πx)cos(2πy)');
xlabel('x'); 
ylabel('y');
set(gca, 'FontSize', 16, 'YDir', 'normal');
print("2D_Rho.png", "-dpng");

% -------------------------------------------------------------------------
% Plot: Potential
% -------------------------------------------------------------------------

figure;
imagesc(x, y, phi);
axis image; colorbar;
title('\phi(x,y)');
xlabel('x'); ylabel('y');
set(gca, 'FontSize', 16, 'YDir', 'normal');
print("2D_Potential.png", "-dpng");


% -------------------------------------------------------------------------
% Plot: Electric Field Vectors (normalized)
% -------------------------------------------------------------------------

figure;
imagesc(x, y, E_mag); 
axis image; hold on;
colorbar;
title('Electric Field Magnitude and Direction');
xlabel('x'); ylabel('y');

step = 4; % Downsample for clarity
Xs = X(1:step:end, 1:step:end);
Ys = Y(1:step:end, 1:step:end);
Exs = Ex(1:step:end, 1:step:end);
Eys = Ey(1:step:end, 1:step:end);
Ems = E_mag(1:step:end, 1:step:end);
quiver(Xs, Ys, Exs./Ems, Eys./Ems, 0.5, 'k');

set(gca, 'FontSize', 16, 'YDir', 'normal');
print("2D_ElectricField_Vectors.png", "-dpng");

% -------------------------------------------------------------------------
% Plot: Charge Flux
% -------------------------------------------------------------------------

figure;
J_mag = sqrt(Jx.^2 + Jy.^2);     % Flux magnitude
imagesc(x, y, J_mag);            % Background color map
axis image; hold on;
colorbar;
title('Charge Flux J(x,y) = μρE');
xlabel('x'); ylabel('y');

step = 4;  % Downsample for clarity
Xs = X(1:step:end, 1:step:end);
Ys = Y(1:step:end, 1:step:end);
Jxs = Jx(1:step:end, 1:step:end);
Jys = Jy(1:step:end, 1:step:end);
Jms = J_mag(1:step:end, 1:step:end);

% Normalize arrows for consistent size
quiver(Xs, Ys, Jxs./Jms, Jys./Jms, 0.5, 'k');

set(gca, 'FontSize', 16, 'YDir', 'normal');
print("2D_ChargeFlux_Map.png", "-dpng");


% -------------------------------------------------------------------------
% Plot: Rate of Change of Charge
% -------------------------------------------------------------------------
figure;
imagesc(x, y, rho_dot);
axis image; colorbar;
title('Rate of Change \rhȯ(x,y) = -∇·J');
xlabel('x'); ylabel('y');
set(gca, 'FontSize', 16, 'YDir', 'normal');
print("2D_RhoDot.png", "-dpng");

% -------------------------------------------------------------------------
% Plot: Old vs New Charge Density
% -------------------------------------------------------------------------

figure;

subplot(2,2,1)
imagesc(x, y, rho);
axis image; colorbar;
title('Old \rho(x,y)');
set(gca, 'YDir', 'normal', 'FontSize', 16);

subplot(2,2,2)
imagesc(x, y, new_rho);
axis image; colorbar;
title('New \rho(x,y)');
set(gca, 'YDir', 'normal', 'FontSize', 16);

% Difference map: Δρ = new_rho - rho
subplot(2,2,3)
imagesc(x, y, new_rho - rho);
axis image; colorbar;
colormap(viridis);
title('\Delta\rho(x,y) = New - Old');
set(gca, 'YDir', 'normal', 'FontSize', 16);

print("2D_Rho_Comparison_Diff.png", "-dpng");

printf("All 2D FFT Poisson computations and plots completed.\n");

