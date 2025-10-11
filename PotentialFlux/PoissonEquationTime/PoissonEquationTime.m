#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window



clear; clc; close all

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
Nx = 256;                   % grid points
Lx = 1.0;                   % domain length
dx = Lx / Nx;
x = linspace(0, Lx - dx, Nx);
epsilon0 = 8.854e-12;       % vacuum permittivity
mu = 1e-13;                 % mobility
D  = 1e-4;                  % diffusion coefficient
dt = 3e-2;
Nt = 5000;                  % total steps
plot_interval = 50;         % plot/save frequency

% -------------------------------------------------------------------------
% Initial condition (periodic sine)
% -------------------------------------------------------------------------
rho = sin(2*pi*x);          % perfectly periodic
%rho = rho - mean(rho);      % ensure neutrality

fprintf('ρ(0)=%.3e, ρ(L)=%.3e, diff=%.3e\n', rho(1), rho(end), rho(1)-rho(end));

phi = zeros(1, Nx);

% -------------------------------------------------------------------------
% FFT wavenumbers
% -------------------------------------------------------------------------
k  = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];
k2 = k.^2;  k2(1) = 1;      % avoid division by zero

% -------------------------------------------------------------------------
% Storage
% -------------------------------------------------------------------------
rho_snapshots = [];
phi_snapshots = [];
time_labels   = [];
convergence   = zeros(1, Nt);
energy        = zeros(1, Nt);
rho_old       = rho;

% -------------------------------------------------------------------------
% TIME LOOP
% -------------------------------------------------------------------------
for t = 1:Nt
    % ---- Neutrality before Poisson ----
    rho = rho - mean(rho);

    % ---- Solve Poisson: φ_xx = -ρ/ε0 ----
    rho_hat = fft(rho);
    phi_hat = rho_hat ./ (epsilon0 * k2);   % <-- correct sign
    phi_hat(1) = 0;                          % remove DC mode
    phi_hat(end) = 0;                        % remove Nyquist
    phi = real(ifft(phi_hat));

    % ---- Electric field: E = -φ_x ----
    E_hat = 1i * k .* phi_hat;
    E_hat(end) = 0;
    E = -real(ifft(E_hat));

    % ---- Drift & diffusion flux ----
    J_drift = mu * rho .* E;                 % μρE (since E already = -φx)

    grad_rho = real(ifft(1i * k .* fft(rho)));
    J_diff  = -D * grad_rho;

    % J_total = J_drift + J_diff;
    J_total = J_drift;
    % J_total = J_diff;

    % ---- Continuity equation: ρ_t = -∂J/∂x ----
    rho_dot = -real(ifft(1i * k .* fft(J_total)));
    rho = rho + dt * rho_dot;
    rho = rho - mean(rho);                   % keep charge neutrality

    % ---- Diagnostics ----
    convergence(t) = sqrt(sum((rho - rho_old).^2) / Nx);
    energy(t) = 0.5 * epsilon0 * sum(E.^2) * dx;
    rho_old = rho;

    % ---- CFL info ----
    if mod(t, plot_interval) == 0
        Emax = max(abs(E));
        driftCFL = mu * Emax * dt / dx;
        diffCFL  = D * dt / dx^2;
        fprintf("Step %d: driftCFL=%.3e, diffCFL=%.3e\n", t, driftCFL, diffCFL);
        Q = sum(rho)*dx;
        fprintf("t=%.3e s, Total charge = %.3e\n", t*dt, Q);
    end

    % ---- Plot + Save periodically ----
    if mod(t, plot_interval) == 0 || t == 1
        rho_snapshots = [rho_snapshots; rho];
        phi_snapshots = [phi_snapshots; phi];
        time_labels   = [time_labels; t*dt];
        cmap = jet(size(rho_snapshots,1));

        % -------------------------------------------------------------
        % Figure 1: Charge & Potential evolution
        % -------------------------------------------------------------
        fig1 = figure(1, "visible", "off"); clf;
        subplot(2,1,1); hold on;
        for k1 = 1:size(rho_snapshots,1)
            plot(x, rho_snapshots(k1,:), 'Color', cmap(k1,:), 'LineWidth', 1.3);
        end
        xlabel('x'); ylabel('\rho(x)');
        title('Charge density evolution  \rho(x,t)  (Periodic BC)');
        grid on; ylim auto;

        subplot(2,1,2); hold on;
        for k1 = 1:size(phi_snapshots,1)
            plot(x, phi_snapshots(k1,:), 'Color', cmap(k1,:), 'LineWidth', 1.3);
        end
        xlabel('x'); ylabel('\phi(x)');
        title('Potential evolution  \phi(x,t)');
        grid on; ylim auto;

        print(fig1, "Charge_Potential_Evolution.png", "-dpng");
        fprintf('Saved Charge_Potential_Evolution.png at step %d\n', t);

        % -------------------------------------------------------------
        % Figure 2: Convergence & Energy
        % -------------------------------------------------------------
        fig2 = figure(2, "visible", "off"); clf;
        subplot(2,1,1);
        semilogy((1:t)*dt, convergence(1:t), 'b', 'LineWidth', 1.3);
        xlabel('Time (s)'); ylabel('||Δρ||₂');
        title('Convergence of charge density'); grid on;

        subplot(2,1,2);
        plot((1:t)*dt, energy(1:t), 'r', 'LineWidth', 1.3);
        xlabel('Time (s)'); ylabel('Energy (J)');
        title('Total electrostatic energy'); grid on;

        print(fig2, "Convergence_Energy.png", "-dpng");
    end

    % ---- Stop if steady ----
    if convergence(t) < 1e-10
        fprintf('Equilibrium reached at step %d (t=%.3g s)\n', t, t*dt);
        break;
    end
end

fprintf('✅ Simulation complete. Final plots saved.\n');

