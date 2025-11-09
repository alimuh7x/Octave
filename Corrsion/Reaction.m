#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path
fprintf("===============================================================\n");
% === Corrosion hydrolysis reaction rates ===
clc; clear; format long;

% ---- Constants ----
Kb1 = 5e2;       % m^3/(mol·s)
Kb2 = 5e3;       % m^3/(mol·s)
K1  = 1.625e-4;  % mol/m^3
K2  = 1.0e-8;    % mol^2/m^6
cH2O = 5.55e4;   % mol/m^3 (water concentration)
cH = 1e-4;       % mol/m^3  (pH = 7)
cOH = K2 / cH;   % mol/m^3 from equilibrium
cFeOH = 0.0;     % mol/m^3 initial Fe(OH)+
ctot = cH2O;     % total concentration (mol/m^3)

% ---- Choose mole fraction of Fe2+ ----
xM = 0.1;                      % metal ion mole fraction
cM = xM * ctot;                % convert to mol/m^3

% ---- Reaction rates (Eq. 46 set) ----
RM  = -Kb1 * (K1 * cM - cFeOH * cH);
RFeOH =  Kb1 * (K1 * cM - cFeOH * cH) + Kb2 * (K2 - cH * cOH);
RH  =  Kb1 * (K1 * cM - cFeOH * cH) + Kb2 * (K2 - cH * cOH);
ROH = -Kb2 * (K2 - cH * cOH);

% ---- Print results ----
fprintf("Fe2+ mole fraction xM = %.3f\n", xM);
fprintf("Fe2+ concentration  cM = %.3e mol/m^3\n", cM);
fprintf("---------------------------------------------------------------\n");
fprintf("R_M        = %.3e mol/m^3/s (Fe2+ consumption)\n", RM);
fprintf("R_H        = %.3e  mol/m^3/s (H+ generation)\n", RH);
fprintf("R_OH       = %.3e mol/m^3/s (OH- change)\n", ROH);
fprintf("R_FeOH     = %.3e  mol/m^3/s (Fe(OH)+ formation)\n", RFeOH);

clc; clear; format long;

% ---- Constants ----
Kb1 = 5e2;       % m^3/(mol·s)
Kb2 = 5e3;       % m^3/(mol·s)
K1  = 1.625e-4;  % mol/m^3
K2  = 1.0e-8;    % mol^2/m^6
cH2O = 5.55e4;   % mol/m^3 (solvent concentration)
cH = 1e-4;       % mol/m^3  (pH = 7)
cOH = K2 / cH;   % mol/m^3
cFeOH = 0.0;     % mol/m^3
ctot = cH2O;     % total concentration (mol/m^3)

% ---- Choose Fe2+ mole fraction ----
xM = 0.01;
cM = xM * ctot;

% ---- Reaction rates in concentration form ----
RM  = -Kb1 * (K1 * cM - cFeOH * cH);
RH  =  Kb1 * (K1 * cM - cFeOH * cH) + Kb2 * (K2 - cH * cOH);
ROH = -Kb2 * (K2 - cH * cOH);
RFeOH = -RM + ROH;

% ---- Convert to mole-fraction form ----
x_dot_M    = RM / cH2O;
x_dot_H    = RH / cH2O;
x_dot_OH   = ROH / cH2O;
x_dot_FeOH = RFeOH / cH2O;

fprintf("---------------------------------------------------------------\n");
% ---- Display results ----
fprintf("Fe2+ mole fraction xM = %.3f\n", xM);
fprintf("x_dot_M    = %.3e  (Fe2+)\n", x_dot_M);
fprintf("x_dot_H    = %.3e  (H+)\n", x_dot_H);
fprintf("x_dot_OH   = %.3e  (OH-)\n", x_dot_OH);
fprintf("x_dot_FeOH = %.3e  (Fe(OH)+)\n", x_dot_FeOH);

fprintf("===============================================================\n");

dt = 1e-6;               % time step (s)
time = 0:dt:1e-4;          % simulate 1 seconds
cH = 1e-4;              % initial H+ (mol/m^3)
cOH = K2 / cH;

for n = 1:length(time)
    
    RM  = -Kb1 * (K1*cM - cFeOH*cH);
    RFeOH = -RM;
    RH  =  Kb1*(K1*cM - cFeOH*cH) + Kb2*(K2 - cH*cOH);
    ROH = -Kb2*(K2 - cH*cOH);
    
    cM    += RM * dt;
    cFeOH += RFeOH * dt;
    cH    += RH * dt;
    cOH    = K2 / cH;

    % update
    pH = -log10(cH/1000);     % convert to mol/L
    % reaction rate
    if mod(n,10) == 0 
        fprintf("t=%.1fs cOH = %.3e cH = %.3e  pH= %.2f\n", time(n), cOH, cH, pH);
    end
end

fprintf("===============================================================\n");
