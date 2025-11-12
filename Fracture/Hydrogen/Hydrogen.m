#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path

% ============================================================
% Hydrogen-induced reduction of surface energy (complete plot)
% ============================================================

clear; clc; close all;

% --- Constants ---
R = 8.314;           % J/mol·K
T = 300;             % K
delta_gb = 30e3;     % J/mol
chi = 0.41;

% --- Concentration range ---
C = logspace(-8, -1, 600);

% --- Langmuir–McLean isotherm ---
theta = C ./ (C + exp(-delta_gb ./ (R .* T)));

% --- Degraded fracture/surface energy ---
Gc_ratio = 1 - chi .* theta;

% --- Plot ---
figure('Position',[100 100 800 500]);  % larger figure window
semilogx(C, Gc_ratio, 'k--', 'LineWidth', 2); hold on;
semilogx(C, theta, 'b-', 'LineWidth', 2);
legend('G_c(C)/G_c(0)', '\theta(C)', 'Location', 'southeast');
xlabel('Hydrogen concentration  C  (mole fraction)', 'FontSize', 12);
ylabel('Normalized value', 'FontSize', 12);
title('Hydrogen-induced degradation of fracture energy', 'FontSize', 12);
grid on; box on;
set(gca, 'LineWidth', 1.2, 'FontSize', 12);

% --- Axis control ---
xlim([1e-8 1e-1]);
ylim([-0.02 1.05]);
pbaspect([1.6 1 1]);  % balanced aspect ratio
set(gcf, 'PaperPositionMode', 'auto'); % ensure complete saving

% --- Save ---
filename = 'Hydrogen_FractureEnergy_vs_Concentration.png';
print(filename, '-dpng', '-r300');
fprintf('✅ Plot saved completely as: %s\n', filename);
