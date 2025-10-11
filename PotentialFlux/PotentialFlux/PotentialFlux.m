#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window

Nx = 101;                % number of grid points
L  = 1e-6;               % domain length (m)
dx = L / (Nx - 1);
x  = linspace(0, L, Nx); % grid points

D  = 1e-9;               % diffusion coefficient (m^2/s)
u  = 5e-8;               % mobility (m^2/V·s)
z  = 1;                 % charge number (Cl-)
F  = 96485;              % Faraday constant (C/mol)

## --- Tanh-shaped potential (well or step) ---
phi0 = 0.05;                  % amplitude (V)
width = L/20;                 % interface width
center = L/2;                 % center position

## Potential profile
phi = phi0 * tanh((x - center) / width);  % electric potential (V)

c = 1;

dphi_dx = gradient(phi, dx);   % dφ/dx
dc_dx   = gradient(c, dx);     % dc/dx

% J = -D ∇c - z*u*c*∇φ
J = - z .* u .* c .* dphi_dx;

% dc/dt = -∇·J = -dJ/dx
dJ_dx = gradient(J, dx);
c_dot = -dJ_dx;

figure;
subplot(3,1,1);
plot(x*1e6, phi, 'LineWidth', 2);
xlabel('x (µm)'); ylabel('\phi (V)');
title('Electric Potential');

subplot(3,1,2);
plot(x*1e6, J, 'LineWidth', 2);
xlabel('x (µm)'); ylabel('J (mol/m^2·s)');
title('Ionic Flux (Fe^2+)');

subplot(3,1,3);
plot(x*1e6, c_dot, 'LineWidth', 2);
xlabel('x (µm)'); ylabel('dc/dt (mol/m^3·s)');
title('Rate of Concentration Change');

print ("cooling_curve.png", "-dpng");  % saves to PNG
