#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../src"));  % add src folder to path


% ============================================================
% Griffith energy release rate demonstration
% ============================================================

clc; clear; close all;

% --- Material & geometry parameters ---
E  = 210e9;       % Young's modulus (Pa)
sigma = 100e6;    % Applied stress (Pa)
delta = 0.001;    % Applied displacement (m)
B = 0.01;         % Thickness (m)

% --- Crack length range ---
a = linspace(0.001, 0.05, 200);  % crack length (m)

% --- Energy release rate formulas ---
% Fixed load (stress-controlled): G ~ sigma^2 * pi * a / E
G_load = (sigma.^2 .* pi .* a) / E;

% Fixed displacement (strain-controlled): G ~ (delta^2 * E * pi) / (a)
% Decreases with a (since same displacement over larger crack)
G_disp = (delta.^2 * E * pi) ./ a;

% --- Critical fracture energy (material property) ---
Gc = 1000;  % J/m^2 (typical for brittle material)

Subplot(a, G_load, G_disp, 'a', 'Fixed Load', 'Fixed Displacement', 'EnergyReleaseRate.png');
% --- Text output ---
fprintf("At G = Gc = %.1f J/m²:\n", Gc);
fprintf("Fixed load: crack length = %.4f m\n", interp1(G_load, a, Gc));
fprintf("Fixed displacement: crack length = %.4f m\n", interp1(G_disp, a, Gc));

