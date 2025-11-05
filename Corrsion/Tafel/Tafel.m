#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src/"));  % add src folder to path

% --- Tafel curve for iron corrosion (Fe → Fe2+ + 2e−) ---
clear; clf;

% Constants
F  = 96485;     % C/mol
R  = 8.314;     % J/mol·K
T  = 298;       % K
n  = 2;         % electrons transferred

% Kinetic parameters
i0 = 1e-6;      % exchange current density (A/cm^2)
alpha_a = 0.5;  % anodic transfer coefficient
alpha_c = 0.5;  % cathodic transfer coefficient
E_eq = -0.44;   % equilibrium potential (V vs SHE for Fe)

% Potential range around E_eq
E = linspace(E_eq - 0.3, E_eq + 0.3, 400);

% Current densities (A/cm^2)
i_a =  i0 * exp((alpha_a * n * F * (E - E_eq)) / (R * T));
i_c = -i0 * exp((-alpha_c * n * F * (E - E_eq)) / (R * T));

% Net current
i_net = i_a + i_c;
% --- Plot log|i| vs E (Tafel plot) ---
figure;
semilogy(E, abs(i_net), 'k', 'LineWidth', 2); hold on;
semilogy(E, abs(i_a), '--r', 'LineWidth', 1.5);
semilogy(E, abs(i_c), '--b', 'LineWidth', 1.5);
xlabel('Potential E (V vs SHE)');
ylabel('log_{10}|Current density| (A/cm^2)');
title('Theoretical Tafel Plot for Iron');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickLength', [0.01, 0.01]);
set(gca, 'Position', [0.15 0.15 0.75 0.75]);  % manual tightening if needed
legend('Total current', 'Anodic', 'Cathodic');
print("Tafel plot.png", '-dpng', '-r300');

% --- Plot |i| vs E (linear scale) ---
figure;
plot(E, abs(i_net), 'k', 'LineWidth', 2); hold on;
plot(E, abs(i_a), '--r', 'LineWidth', 1.5);
plot(E, abs(i_c), '--b', 'LineWidth', 1.5);
xlabel('Potential E (V vs SHE)');
ylabel('|Current density| (A/cm^2)');
title('Current Density vs Potential for Iron');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickLength', [0.01, 0.01]);
set(gca, 'Position', [0.15 0.15 0.75 0.75]);  % manual tightening if needed
legend('Total current', 'Anodic', 'Cathodic');
print("Potential plot.png", '-dpng', '-r300');

figure;
plot(E, i_net, 'k', 'LineWidth', 2); hold on;
plot(E, i_a, '--r', 'LineWidth', 1.5);
plot(E, i_c, '--b', 'LineWidth', 1.5);
xlabel('Potential E (V vs SHE)');
ylabel('Current density (A/cm^2)');
title('Current Density vs Potential for Iron');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickLength', [0.01, 0.01]);
set(gca, 'Position', [0.15 0.15 0.75 0.75]);  % manual tightening if needed
legend('Total current', 'Anodic', 'Cathodic');
print("2Potential plot.png", '-dpng', '-r300');

