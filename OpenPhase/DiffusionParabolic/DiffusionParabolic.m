#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("./../../src"));  % add src folder to path
warning("off", "all");

clear; clc;

% --- Diffusion rate and timestep estimator ---

% Input parameters
D  = 1e-18;       % diffusivity (m^2/s)
dx = 0.5e-9;        % grid spacing (m)
M = 5e-18;
dt = 0.01;
c_increment = 1e-6
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

fprintf('\n');


Nx = 100;
dx = 0.5e-9;              % 0.5 nm grid
x  = (0:Nx-1) * dx;       % physical coordinate in meters
i = 1:Nx;               % index array

myprint("dx : ", dx);

dt = 0.001;
nsteps = 5000000;

% Phase-field and material parameters
w = 8 * dx;                 % interface width

A = 4e8;

    % --- Phase-field evolution (Allen–Cahn type) ---
L_phi   = 1e-16;          % phase-field mobility
sigma   = 0.5;            % gradient energy coefficient



R = 8.3145;               % J/(mol K)
D_alpha = 1e-19;
D_beta = 1e-19;
Vm = 7.4e-6;               % m^3/mol

% --- Phase-field profile ---
ProfilePosition =  x-mean(x);;
phi = 0.5 * (1 - tanh(3 * ProfilePosition / w));

phi_alpha = phi;
phi_beta  = 1 - phi;

C_alpha = 0.5;
C_beta  = 0.8;

Ceq_alpha = 0.6;
Ceq_beta  = 0.6;

c   = phi_alpha .* C_alpha + phi_beta .* C_beta;
c_initial   = phi_alpha .* C_alpha + phi_beta .* C_beta;
% --- Total composition ---
C_eq   = phi_alpha .* Ceq_alpha + phi_beta .* Ceq_beta;

M = (D_alpha .* phi_alpha + D_beta .* phi_beta )/ (R * 300);

%mu  = A * (c - Ceq) * Vm;



%SinglePlot(x, c, 'x', 'c(x)', 'Initial_Concentration.png');
%SinglePlot(x, Ceq, 'x', 'C_eq', 'Initial_Equilibrium_Concentration.png');
%SinglePlot(x, phi, 'x', 'phi', 'Phase_Field.png');
%SinglePlot(x, mu, 'x', 'mu', 'InitialPotential.png');
%SinglePlot(x, M, 'x', 'Mobility', 'Mobility.png');

Phi_plot = SnapshotPlotter('1_PhaseEvolution.png');
C_plot = SnapshotPlotter('1_CompositionEvolution.png');
dc_plot = SnapshotPlotter('2_dc.png');
PotentialPlot = SnapshotPlotter('3_PotentialGrad.png');

f_plot = SnapshotPlotter('1_Drivingforce.png');

P = CurvePlot("physical.png");

for step = 1:nsteps
    % Chemical potential (single field)
     mu = A * Vm * (c - C_eq);


     % Laplacian with Neumann BC (zero flux)
     mu_xx = zeros(size(mu));
     mu_xx(2:end-1) = (mu(3:end) - 2*mu(2:end-1) + mu(1:end-2)) / dx^2;
     
     % Neumann BCs: zero derivative at boundaries
     mu_xx(1) = (mu(2) - mu(1)) / dx^2;          % left boundary
     mu_xx(end) = (mu(end-1) - mu(end)) / dx^2;  % right boundary
     
     % Update concentration
     dc = M .* mu_xx;
     %c = c + dt * dc;

    f_alpha = 0.5 * A * (c - Ceq_alpha).^2;
    f_beta  = 0.5 * A * (c - Ceq_beta).^2;




    % Visualization
    Delta_f = f_beta - f_alpha;

        % Laplacian of phi
    phi_xx = zeros(size(phi));
    phi_xx(2:end-1) = (phi(3:end) - 2*phi(2:end-1) + phi(1:end-2)) / dx^2;
    phi_xx(1) = (phi(2) - phi(1)) / dx^2;
    phi_xx(end) = (phi(end-1) - phi(end)) / dx^2;

        % Double-well derivative
    dW_dphi = 72 * 2 * phi .* (1 - phi) .* (1 - 2*phi) / w^2;

    % --- Phase-field rate and update ---
    phi_dot = L_phi * ( sigma * (phi_xx - dW_dphi) - (6 / w) .* phi_alpha .* phi_beta .* Delta_f );
    
    phi = phi + dt * phi_dot;
    phi = max(0, min(1, phi));   % clamp to [0,1]



    % Clamp phi to [0,1]
    %phi = max(0, min(1, phi));

    % --- Update phase fractions and equilibrium composition ---
    phi_alpha = 1 - phi;
    phi_beta  = phi;
    C_eq = phi_alpha .* Ceq_alpha + phi_beta .* Ceq_beta;


    if mod(step, 100) == 0
        % Total conservation check
        total_C = trapz(x, c);
        total_C_initial = trapz(x, c_initial);
        Error = abs((total_C - total_C_initial) / total_C_initial) * 100;
        fprintf('=================================================================\n');
        printf('Step %d / %d\n', step, nsteps);
        myprint("Max Laplacian       : ", max(mu_xx));
        myprint("Max c_dot           : ", max(dc));
        myprint("Max c_increment     : ", max(dc) * dt);
        myprint("Initial C Total     : ", total_C_initial);
        myprint("Total C             : ", total_C);
        myprint("Error               : ", Error);
        myprint("phi dot             : ", max(abs(phi_dot)));
        myprint("max Driving Force   : ", max(abs(Delta_f)));
        fprintf('=================================================================\n');


        Subplot(i, phi_dot, Delta_f, 'x', 'phi_dot', 'Driving Force', '1PhaseDot_and_DrivingForce.png');
        phi_phi = phi_alpha .* phi_beta;
        Phi_plot.add_and_plot(i, phi,   'x', 'phi');
        Phi_plot.add_and_plot(i, phi_phi,   'x', 'phi_alpha * phi_beta');
        C_plot.add_and_plot(i, c, C_eq,  'x', 'C', 'C_eq');
        dc_plot.add_and_plot(i, c, c - C_eq,  'x', 'C', 'C - C_eq');
        PotentialPlot.add_and_plot(i, mu, dc,  'x', 'mu(x) ', 'Laplacian of mu * M');

        dg_alpha = c - Ceq_alpha;
        dg_beta  = c - Ceq_beta;
        f_plot.add_and_plot(i, dg_alpha, dg_beta,  'x', 'f_alpha', 'f_beta');
        %P.add_plot(i, c, "x", "C");
        %P.add_plot(i, C_eq, "x", "C_eq");
        %P.add_plot(i, c_initial, "x", "C_Initial");
        %P.save();

    end
end

