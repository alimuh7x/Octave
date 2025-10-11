#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window

clear; clc; close all

% Parameters
Nx = 256;
Lx = 1.0;
dx = Lx / Nx;
x  = linspace(0, Lx - dx, Nx);
D  = 1e-4;
dt = 3e-2;
Nt = 5000;
plot_interval = 100;

% Initial condition (periodic sine)
rho = sin(2*pi*x);

% FFT wavenumbers
k  = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];
k2 = k.^2;

rho_snapshots = [];
lap_snapshots = [];
time_labels   = [];

for t = 1:Nt
    % ---- Compute Laplacian ∂²ρ/∂x² ----
    rho_hat = fft(rho);
    lap_rho = real(ifft(-k2 .* rho_hat));

    % ---- Pure diffusion: ρ_t = D ∂²ρ/∂x² ----
    rho = rho + dt * D * lap_rho;

    % ---- Store and plot ----
    if mod(t, plot_interval) == 0 || t == 1
        rho_snapshots = [rho_snapshots; rho];
        lap_snapshots = [lap_snapshots; lap_rho];
        time_labels = [time_labels; t*dt];
        cmap = jet(size(rho_snapshots,1));

        diffCFL  = D * dt / dx^2;
        fprintf('Step %d, Time %.3f, DiffCFL = %.3f\n', t, t*dt, diffCFL);
        fprintf('----------------------------------------------------------------\n');

        % Plot ρ(x,t)
        fig1 = figure(1, "visible", "off"); clf;
        subplot(2,1,1); hold on
        for i = 1:size(rho_snapshots,1)
            plot(x, rho_snapshots(i,:), 'Color', cmap(i,:), 'LineWidth', 1.2);
        end
        xlabel('x'); ylabel('\rho(x)');
        title('Diffusion of charge density  \rho(x,t)');
        grid on;

        % Plot ∂²ρ/∂x²
        subplot(2,1,2); hold on
        for i = 1:size(lap_snapshots,1)
            plot(x, lap_snapshots(i,:), 'Color', cmap(i,:), 'LineWidth', 1.2);
        end
        xlabel('x'); ylabel('\nabla^2\rho(x)');
        title('Laplacian of charge density  ∂²ρ/∂x²');
        grid on;

        print(fig1, "Diffusion_and_Laplacian.png", "-dpng");
        fprintf('Saved Diffusion_and_Laplacian.png at step %d\n', t);
    end
end

fprintf('✅ Diffusion-only simulation complete.\n');

