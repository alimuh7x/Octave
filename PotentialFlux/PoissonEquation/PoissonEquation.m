#!/usr/bin/env -S octave --no-gui --quiet
set(0, "defaultfigurevisible", "off");  % hide figure window
addpath(genpath("~/.config/nvim/Octave/src"));

% ---------------------------------------------------------------------------
% NOTE: Solving Numerically using finite difference method (SOR)
% ---------------------------------------------------------------------------

clear; clc;

Solver.Jacobi = 1;
Solver.FFT = 2;
Solver.SOR = 3;

solver = Solver.FFT;


% Poisson equation in 1D using SOR: d²φ/dx² = -ρ(x)/ε₀

% Parameters
Nx = 200;
Lx = 1.0;
dx = Lx / (Nx - 1);

epsilon0 = 8.854e-12;

omega = 1.8;       % relaxation factor (1 < ω < 2 for acceleration)

% Grid
x = linspace(0, Lx, Nx);

% Charge density
rho = sin(2 * pi * x);

% Initialize potential (Dirichlet BC)
phi = zeros(1, Nx);

% Iteration parameters
max_iter = 10000;
tol = 1e-9;

switch solver
    case Solver.Jacobi
        % Jacobi iteration (using old values)
        for iter = 1:max_iter
            phi_old = phi;   % keep previous step
            for i = 2:Nx-1
                phi(i) = 0.5 * (phi_old(i+1) + phi_old(i-1) + dx^2 * rho(i) / epsilon0);

                %  if iter == 9999
                %     if i == 99
                %         myprint("Phi(i-1)", phi_old(i-1), "old Phi(i)", phi_old(i), "Phi(i+1)", phi_old(i+1));
                %         myprint("neibhourAvg", 0.5 * (phi_old(i+1) + phi_old(i-1)));
                %         myprint("dx^2 * rho(i) / epsilon0", 0.5 * dx^2 * rho(i) / epsilon0);
                %         myprint("phi_new_Calculated", phi(i));
                %         diff = phi(i) - 0.5 * (phi_old(i+1) + phi_old(i-1) + dx^2);
                %            myprint("diff", diff);
                %         myprint("rho(i)", rho(i));
                %         printf("------------------------------------------------------------------------\n");
                %     end
                %     if i == 100
                %         myprint("Phi(i-1)", phi_old(i-1), "OLD Phi(i)", phi_old(i), "Phi(i+1)", phi_old(i+1));
                %         myprint("neibhourAvg", 0.5 * (phi_old(i+1) + phi_old(i-1)));
                %         myprint("dx^2 * rho(i) / epsilon0", 0.5 * dx^2 * rho(i) / epsilon0);
                %         myprint("phi_new_Calculated", phi(i));
                %         diff = phi(i) - 0.5 * (phi_old(i+1) + phi_old(i-1) + dx^2);
                %            myprint("diff", diff);
                %         myprint("rho(i)", rho(i));
                %         printf("------------------------------------------------------------------------\n");
                %     end
                %     if i == 101
                %         myprint("Phi(i-1)", phi_old(i-1), "OLD Phi(i)", phi_old(i), "Phi(i+1)", phi_old(i+1));
                %         myprint("neibhourAvg", 0.5 * (phi_old(i+1) + phi_old(i-1)));
                %         myprint("dx^2 * rho(i) / epsilon0", 0.5 * dx^2 * rho(i) / epsilon0);
                %         myprint("phi_new_Calculated", phi(i));
                %         diff = phi(i) - 0.5 * (phi_old(i+1) + phi_old(i-1) + dx^2);
                %            myprint("diff", diff);
                %         myprint("rho(i)", rho(i));
                %         printf("=========================================================================\n");
                %         pause;
                %     end
                %  end
            end
            % Periodic BCs
            phi(1) = 0.5 * (phi_old(2) + phi_old(end) + dx^2 * rho(1) / epsilon0);
            phi(end) = 0.5 * (phi_old(1) + phi_old(end-1) + dx^2 * rho(end) / epsilon0);
        

            % Convergence check

            maxdiff = max(abs(phi - phi_old));
            if mod(iter, 500) == 0
                fprintf("Iter %d: maxdiff = %.3e\n", iter, maxdiff);
            end
            if maxdiff < tol
                fprintf("Jacobi converged in %d iterations\n", iter);
                break;
            end
        end
        phi = phi - mean(phi);  % remove DC offset

    case Solver.SOR
         for iter = 1:max_iter
             maxdiff = 0;
             for i = 2:Nx-1
                 phi_new = 0.5 * (phi(i+1) + phi(i-1) + dx^2 * rho(i) / epsilon0);
                 phi(i) = (1 - omega) * phi(i) + omega * phi_new;
                 maxdiff = max(maxdiff, abs(phi_new - phi(i)));
             end

             % Print iteration info
             if mod(iter, 500) == 0
                 printf("Iteration %d: max change = %.6e\n", iter, maxdiff);
             end

             if maxdiff < tol
                 printf("Converged in %d iterations\n", iter);
                 break;
             end
         end
    case Solver.FFT
        rho_hat = fft(rho);
        k  = (2*pi/Lx) * [0:(Nx/2-1) -Nx/2:-1];  % periodic wavenumbers
        k2 = k.^2;
        k2(1) = 1;   % avoid division by zero at k = 0
        phi_hat = rho_hat ./ (epsilon0 * k2);
        phi_hat(1) = 0;  % set mean potential to zero
        phi = real(ifft(phi_hat));
end




% Plot results
subplot(2,1,1)
plot(x, rho, 'r', 'LineWidth', 1.5)
xlabel('x'); ylabel('\rho(x)')

title('Charge density \rho(x) = cos(2πx)')
set(gca, 'XTick', [0, 0.25, 0.5, 0.75, 1]);
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
grid on

subplot(2,1,2)
plot(x, phi, 'b', 'LineWidth', 1.5)
xlabel('x'); ylabel('\phi(x)')
title('Potential \phi(x)')
set(gca, 'XTick', [0, 0.25, 0.5, 0.75, 1]);
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
grid on

print ("Numerical_SOR.png", "-dpng");  % saves to PNG


% -------------------------------------------------------------------------
%  NOTE: Compute electric field and movement from potential
% -------------------------------------------------------------------------

E = -gradient(phi, dx);       % Electric field (negative gradient)
mu = 5e-10;                    % Mobility (m^2/V·s) — example value

J = mu * rho .* E;


% -------------------------------------------------------------------------
%  NOTE: Plot field, velocity, and flux
% -------------------------------------------------------------------------

figure

subplot(2,1,1)
plot(x, E, 'm', 'LineWidth', 1.5)
xlabel('x'); ylabel('E(x)')
title('Electric Field E(x) = -d\phi/dx')
set(gca, 'XTick', [0, 0.25, 0.5, 0.75, 1]);
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
grid on


subplot(2,1,2)
plot(x, J, 'k', 'LineWidth', 1.5)
xlabel('x'); ylabel('J(x)')
title('Charge Flux J(x) =  M ρ·E')
set(gca, 'XTick', [0, 0.25, 0.5, 0.75, 1]);
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
grid on

print("Charge_Movement_from_Potential.png", "-dpng");

% -------------------------------------------------------------------------
%  NOTE: Compute rho_dot (charge accumulation rate)
% -------------------------------------------------------------------------

rho_dot = -gradient(J, dx);
rho_dot = rho_dot - mean(rho_dot);  % remove drift
new_rho = rho + 0.01 * rho_dot;    % smaller time step

figure


subplot(2,1,1)
plot(x, rho_dot, 'k', 'LineWidth', 1.5)
title('Rate of Change \rhȯ(x) = -dJ/dx')
xlabel('x'); ylabel('\rhȯ(x)')
set(gca, 'XTick', [0, 0.25, 0.5, 0.75, 1]);
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
grid on

subplot(2,1,2)
plot(x, rho, 'r', 'LineWidth', 1.5); hold on
plot(x, new_rho, 'b--', 'LineWidth', 1.5);
title('Old vs New Charge Density \rho(x)')
set(gca, 'XTick', [0, 0.25, 0.5, 0.75, 1]);
set   (gca, 'FontSize', 14);
set   (gca, 'LineWidth', 1.5);
set   (gca, 'TickLength', [0.02, 0.02]);
xlabel('x'); ylabel('\rho(x)')
legend('Old \rho', 'New \rho')
grid on

print("Charge_Rho_Comparison.png", "-dpng");

printf("All computations and plots completed.\n");
