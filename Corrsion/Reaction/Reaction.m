#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path
% === Corrosion hydrolysis reaction rates ===
clc; clear; format long;
fprintf("===============================================================\n");
fprintf("Fe hydrolysis kinetics (Hageman & Martínez-Pañeda, 2023)\n");
fprintf("===============================================================\n");

% -------------------- CONSTANTS --------------------

Kw      = 1.0e-8;        % mol^2/m^6  (water ionization constant)
k_eq    = 1.0e6;         % m^3/(mol*s) (water equilibrium penalty)
k_fe    = 0.1;           % m^3/(mol*s) (Fe2+ + H2O -> FeOH+ + H+)
k_fe_r  = 1e-3;          % m^3/(mol*s) (reverse)
k_feoh  = 1e-3;          % s^-1        (FeOH+ + H2O -> Fe(OH)2 + H+)

% -------------------- INITIAL VALUES --------------------

cFe2   = 100;            % mol/m^3   Fe2+ initial
cFeOH  = 0.0;             % mol/m^3   FeOH+ initial
cFeOH2 = 0.0;             % mol/m^3   Fe(OH)2 initial
cH     = 1e-4;            % mol/m^3   (pH ≈ 7)
cOH    = Kw / cH;         % mol/m^3
cH2O   = 5.55e4;          % mol/m^3   water

dt   = 1e-7;             % time step [s]
tEnd = 0.1;              % total time [s]
Nt   = round(tEnd/dt);

% -------------------- STORAGE ARRAYS --------------------
time      = zeros(1,Nt);
pH_t      = zeros(1,Nt);

cFe2_t    = zeros(1,Nt);
cFeOH_t   = zeros(1,Nt);
cFeOH2_t  = zeros(1,Nt);
cH_t      = zeros(1,Nt);
cOH_t     = zeros(1,Nt);

xFe2_t    = zeros(1,Nt);
xFeOH_t   = zeros(1,Nt);
xFeOH2_t  = zeros(1,Nt);
xH_t      = zeros(1,Nt);
xOH_t     = zeros(1,Nt);

Fe_total0 = cFe2 + cFeOH + cFeOH2;

% -------------------- TIME INTEGRATION --------------------
for n = 1:Nt
    % --- water autoionization equilibrium ---
    Rw = k_eq*(Kw - cH*cOH);

    % --- Fe hydrolysis reactions ---
    R_Fe2  = -k_fe*cFe2 + k_fe_r*cFeOH*cH;
    R_FeOH =  k_fe*cFe2 - cFeOH*(k_feoh + k_fe_r*cH);
    R_H_fe =  k_fe*cFe2 - cFeOH*(k_fe_r*cH - k_feoh);

    % --- Update concentrations ---
    cFe2   += R_Fe2  * dt;
    cFeOH  += R_FeOH * dt;
    cFeOH2 += k_feoh * cFeOH * dt;
    cH     += (R_H_fe + Rw) * dt;
    cH      = max(cH, 1e-12);
    cOH     = Kw / cH;

    % --- Mole fractions ---
    xFe2   = cFe2   / cH2O;
    xFeOH  = cFeOH  / cH2O;
    xFeOH2 = cFeOH2 / cH2O;
    xH     = cH     / cH2O;
    xOH    = cOH    / cH2O;

    % --- Store ---
    time(n)     = n*dt;
    pH_t(n)     = -log10(cH*1e-3);
    cFe2_t(n)   = cFe2;
    cFeOH_t(n)  = cFeOH;
    cFeOH2_t(n) = cFeOH2;
    cH_t(n)     = cH;
    cOH_t(n)    = cOH;
    xFe2_t(n)   = xFe2;
    xFeOH_t(n)  = xFeOH;
    xFeOH2_t(n) = xFeOH2;
    xH_t(n)     = xH;
    xOH_t(n)    = xOH;

    Fe_total = cFe2 + cFeOH + cFeOH2;

    % --- Occasionally print status ---
    if mod(n, Nt/20) == 0
        fprintf("t = %.2fs | pH = %.2f | cH = %.3e | Fe2+ = %.3e | FeOH+ = %.3e | FeOH2 = %.3e\n", ...
                time(n), pH_t(n), xH, xFe2, xFeOH, xFeOH2);
        fprintf("    Total Fe = %.3e (ΔFe = %.1e)\n", Fe_total, Fe_total - Fe_total0);
        fprintf("---------------------------------------------------------------\n");
    end
end

fprintf("===============================================================\n");
fprintf("Final pH = %.2f\n", pH_t(end));
fprintf("Final Fe species:\n  Fe2+ = %.3e\n  FeOH+ = %.3e\n  Fe(OH)2 = %.3e mol/m^3\n", ...
        cFe2, cFeOH, cFeOH2);
fprintf("===============================================================\n");


% -------------------- OPTIONAL PLOTS --------------------

Subplot(time, pH_t, cFeOH2_t, 'x', 'PH',  'FeOH2', '1PH_FeOH.png', true);
Subplot(time, cFe2_t, cFeOH_t,'x', 'Fe2', 'FeOH',  '2Fe_FeOH.png', true);
Subplot(time, cH_t, cOH_t,    'x', 'H',   'OH',    '3H_OH.png', true);

Subplot(time, pH_t, xFeOH2_t, 'x', 'PH',  'xFeOH2',  '1PH_FeOH_moleFraction.png', true);
Subplot(time, xFe2_t, xFeOH_t,'x', 'xFe2', 'xFeOH',  '2Fe_FeOH_moleFraction.png', true);
Subplot(time, xH_t, xOH_t,    'x', 'xH',   'xOH',    '3H_OH_moleFraction.png', true);

